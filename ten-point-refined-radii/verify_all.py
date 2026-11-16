"""Exact marked-radius row verification for the refined ten-point kernel."""
import argparse
from fractions import Fraction as F
import importlib.util
import json
from hashlib import sha256
from math import isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REFINED = ROOT / 'ten-point-banded-kernel-refined'
BANDS = ((F(27, 500), F(217, 1000)), (F(37, 50), F(139, 100)))
EXCLUSIONS = {
    'outer-low': (1, F(37, 50), F(169, 200)),
    'outer-high': (1, F(164, 125), F(139, 100)),
    'inner-low': (0, F(27, 500), F(37, 500)),
    'inner-high': (0, F(43, 250), F(217, 1000)),
}
FLOOR = F(61, 200000)
ROW_FLOOR = 9 * FLOOR
GAP = F(1, 10**6)
C_OUTSIDE = F(3, 40000)
C_INNER = F(9, 100000)


def need(test, why):
    if not test:
        raise ArithmeticError(why)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build():
    data = json.loads((REFINED / 'kernel.json').read_bytes())
    report = json.loads((REFINED / 'verification_report.json').read_bytes())
    banded = load('banded', ROOT / 'ten-point-banded-kernel/verify_all.py')
    row = load('row', ROOT / 'ten-point-ward-radii/verify_all.py')
    ward = banded.upper
    need(report['passed'] is True and report['unresolved'] == 0
         and report['exact_inner_count'] == 2 and report['exact_outer_count'] == 8,
         'replayed upstream refined kernel and populations')
    inner = report['inner_pair']
    need(inner['applies_to_every_global_maximizer'] is True and inner['unresolved'] == 0
         and tuple(map(F, inner['inner_squared_radius_open_interval'])) == BANDS[0]
         and tuple(map(F, inner['inner_dot_product_open_interval'])) == (F(-27, 200), F(-1, 10))
         and F(inner['inner_midpoint_squared_strict_upper']) == F(1, 40)
         and F(inner['angle_cosine_strict_upper']) == F(-2, 3)
         and F(inner['inner_squared_distance_strict_upper']) == F(3, 5),
         'verified inherited inner-pair domain')
    banded.verify_dependencies()
    need(banded.BANDS[1] == BANDS[1], 'verified inherited outer band')
    need(F(report['c']) == C_OUTSIDE and F(report['c_inner']) == C_INNER,
         'matched refined slacks')
    need(ward.hx(data['B']) - ward.log.log_lower(ward.witness()) < GAP,
         'strict witness log gap')
    corrections, forces = ward.corrections(data)
    poly, shape = ward.polynomial(data, corrections, forces)
    need(bool(forces), 'defined rational force correction')
    params = banded.parsed(data)
    band_polys = {(i, j): banded.root_polynomial(poly, shape, params, i, j)
                  for i in range(2) for j in range(2)}
    for i, j in band_polys:
        need(banded.check_root(poly, band_polys[i, j], params, i, j) == 6,
             'exact independent band expansion crosscheck')
    return data, banded, row, shape, band_polys


def box_for(interval, neighbor, sign):
    band = BANDS[neighbor]
    bound_square = interval[1] * band[1] * 10000
    cap = isqrt(bound_square.numerator // bound_square.denominator)
    bound = F(cap + (cap * cap < bound_square), 100)
    return interval, band, (-bound, F(0)) if sign == -1 else (F(0), bound)


def expression(base, row, multipliers, mark, neighbor):
    a, b, force, count = multipliers
    balance = F(neighbor == 0) - F(2 - int(mark == 0), 9)
    poly = row.rowpoly(base, (a, b, force))
    difference = row.FLOOR - FLOOR - count * balance
    for key, d_coefficient in (((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), -2)):
        poly[key] = poly.get(key, F(0)) + d_coefficient * difference
    return poly


def can_prune(box, mark, neighbor, banded, row):
    ward = banded.upper
    (u0, u1), (v0, v1), (t0, t1) = box
    if mark == neighbor == 0:
        t_upper_squared = max(t0 * t0, t1 * t1)
        if (t0 >= F(-1, 10) or t1 <= F(-27, 200)
                or u0 + v0 + 2 * t0 >= F(1, 10)
                or 9 * t_upper_squared <= 4 * u0 * v0
                or u0 + v0 - 2 * t1 >= F(3, 5)):
            return True
    return (row.prunable(box, ward.U, ward.DMIN, ward.L)
            or banded.collapse.pair_box_excluded((u0,v0,t0), (u1,v1,t1)))


def replay(key, built, progress=False):
    _, banded, row, shape, expressions = built
    target = HERE / f'{key}.json'
    ward = banded.upper
    pack = ward.Packed(shape, True)
    ward.old.selftest()
    mark, begin, end = EXCLUSIONS[key]
    entries = json.loads(target.read_bytes())
    totals = dict.fromkeys('AP012',0)
    max_depth = 0
    cursor = begin
    roots = [(neighbor, sign) for neighbor in range(2) for sign in (-1,1)]
    for entry in entries:
        lo, hi = map(F, entry['interval'])
        multipliers = tuple(map(F, entry['multipliers']))
        need(lo == cursor and lo < hi <= end and len(multipliers)==4 and len(entry['trees'])==4,
             'all roots of a valid marked-radius slab')
        cursor = hi
        for (neighbor,sign),tree in zip(roots,entry['trees']):
            box = box_for((lo,hi),neighbor,sign)
            coeff = expression(expressions[mark,neighbor],row,multipliers,mark,neighbor)
            stack = [(box,ward.initial(coeff,shape,box),0)]
            counts = dict.fromkeys('AP012',0)
            for action in tree:
                need(stack and action in counts, 'valid tree action')
                box, packed, depth = stack.pop()
                max_depth = max(depth,max_depth)
                counts[action] += 1
                if action == 'P':
                    need(can_prune(box,mark,neighbor,banded,row), 'rigorous necessary-domain exclusion')
                elif action == 'A':
                    need(pack.dominates(packed,box,F(0)), 'exact Bernstein/log tangent lower bound')
                else:
                    stack.extend(banded.split(pack,packed,box,depth,int(action)))
            need(not stack and counts['A']+counts['P']==1+sum(counts[k] for k in '012'),
                 'complete signed root partition')
            for action, count in counts.items():
                totals[action] += count
        if progress:
            print(key, 'exact through', cursor, sum(totals.values()), file=sys.stderr, flush=True)
    need(cursor==end,'complete marked-radius interval')
    slack = (35*C_OUTSIDE+C_INNER) if mark==1 else 36*C_OUTSIDE
    need(ROW_FLOOR-slack>GAP,'row exceeds global remaining-pair budget')
    return {'passed':True, 'interval':list(map(str,(begin,end))), 'marked_band':mark,
            'row_lower':str(ROW_FLOOR), 'remaining_edge_slack':str(slack),
            'witness_log_gap_strict_upper':str(GAP), 'slabs':len(entries),
            'visited':sum(totals.values()), 'counts':totals,
            'maximum_depth':max_depth, 'unresolved':0}


def verify(progress=False):
    manifest = json.loads((HERE / 'source_manifest.json').read_bytes())
    for base, key in ((HERE, 'files'), (ROOT, 'dependencies')):
        for name, digest in manifest[key].items():
            need(sha256((base / name).read_bytes()).hexdigest() == digest, 'fixed input ' + name)
    built = build()
    reports = [replay(key, built, progress) for key in EXCLUSIONS]
    data, banded, _, _, _ = built
    midpoint = load('refined_midpoint', HERE/'midpoint.py')
    geometry = midpoint.verify(data, banded.upper,
                               ((F(37,500), F(43,250)), (F(169,200), F(164,125))),
                               2*(GAP + 44*C_OUTSIDE + C_INNER))
    return {
        'passed': True,
        'applies_to_every_global_maximizer': True,
        'normalization_mass': 9,
        'exact_inner_count': 2,
        'exact_outer_count': 8,
        'inner_squared_radius_open_interval': ['37/500', '43/250'],
        'outer_squared_radius_open_interval': ['169/200', '164/125'],
        'total_nodes': sum(r['visited'] for r in reports),
        'unresolved': 0,
        'covers': reports,
        'midpoint': geometry,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    parser.add_argument('--progress', action='store_true')
    args = parser.parse_args()
    result = verify(args.progress)
    raw = (json.dumps(result, indent=2) + '\n').encode()
    report = HERE / 'verification_report.json'
    if args.record_report:
        report.write_bytes(raw)
    else:
        need(result == json.loads(report.read_bytes()), 'exact replay matches saved report')
    print(raw.decode(), end='')


if __name__ == '__main__':
    main()
