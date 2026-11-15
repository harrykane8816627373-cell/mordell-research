"""Exact marked-row exclusions using the certified two-inner population."""
from argparse import ArgumentParser
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DOMAINS = ((F(0), F(3, 10)), (F(41, 60), F(71, 50)), (F(149, 100), F(1539, 1000)))
EXCLUSIONS = ((F(7, 1000), F(43, 1000), 0), (F(51, 200), F(3, 10), 0),
              (F(41, 60), F(37, 50), 1), (F(139, 100), F(71, 50), 1),
              (F(149, 100), F(61, 40), 2))
FLOOR = F(1, 800)
LOW = F(7, 1000)
INNER = (F(43, 1000), F(51, 200))
INNER_DOT = (F(-31, 200), F(-4, 125))
INNER_MIDPOINT_SQUARED = F(23, 200)


def load(module, path):
    spec = importlib.util.spec_from_file_location(module, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def need(ok, description):
    if not ok:
        raise ArithmeticError(description)


def counted_polynomial(row, polynomial, coefficients, marked, other):
    a, b, c, count = coefficients
    expression = row.rowpoly(polynomial, (a, b, c))
    need(row.FLOOR == FLOOR, 'same row floor')
    adjustment = count*(int(other == 0)-F(2-int(marked == 0), 9))
    for key, coefficient in (((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), -2)):
        expression[key] = expression.get(key, F(0))-coefficient*adjustment
    return expression


def prunable(row, ward, box):
    (ul, uh), (vl, vh), (tl, th) = box
    if uh+vh <= F(1007, 5000):
        return True
    if uh+vh-2*tl <= F(203, 500) and uh+vh+2*th <= F(1, 2500):
        return True
    return row.prunable(box, ward.U, ward.DMIN, ward.L)


def verify_tree(ward, row, pack, poly, shape, box, tree):
    stack = [(box, ward.initial(poly, shape, box), 0)]
    depth = 0
    counts = dict.fromkeys('AP012', 0)
    for node in tree:
        need(bool(stack) and node in counts, 'valid action on an existing interval')
        counts[node] += 1
        box, coefficients, level = stack.pop()
        depth = max(depth, level)
        if node == 'A':
            need(pack.dominates(coefficients, box, F(0)), 'exact Bernstein lower bound')
        elif node == 'P':
            need(prunable(row, ward, box), 'physical or deletion prune')
        else:
            axis = int(node)
            lower, upper = row.split(box, axis)
            first, second = pack.split(coefficients, axis)
            stack.extend(((upper, second, level+1), (lower, first, level+1)))
    need(not stack and counts['A']+counts['P']==1+sum(counts[k] for k in '012'), 'complete signed tree')
    return counts, depth


def inner_pair_contained(box):
    (_, u), (_, v), (t0, t1) = box
    return INNER_DOT[0] < t0 <= t1 < INNER_DOT[1] and u+v+2*t1 < 4*INNER_MIDPOINT_SQUARED


def verify_inner_pair(ward, row, pack, poly, shape):
    tree = (HERE/'inner-pair.txt').read_text(encoding='ascii').strip()
    gap = F(1, 10**6)
    budget = gap+44*ward.C
    box = (INNER, INNER, (-INNER[1], INNER[1]))
    stack = [(box, ward.initial(poly, shape, box), 0)]
    counts = dict.fromkeys('APG012', 0)
    depth = 0
    for node in tree:
        need(bool(stack) and node in counts, 'valid action for an inner pair box')
        counts[node] += 1
        box, coefficients, level = stack.pop()
        depth = max(depth, level)
        if node == 'A':
            need(pack.dominates(coefficients, box, -budget), 'pair deficit exceeds total budget')
        elif node == 'P':
            need(box[0][0]>box[1][1] or prunable(row, ward, box), 'exact symmetric pair prune')
        elif node == 'G':
            need(inner_pair_contained(box), 'both dot-product and midpoint restrictions')
        else:
            axis = int(node)
            lower, upper = row.split(box, axis)
            first, second = pack.split(coefficients, axis)
            stack.extend(((upper, second, level+1), (lower, first, level+1)))
    need(not stack and counts['A']+counts['P']+counts['G']==1+sum(counts[k] for k in '012'),
         'complete inner pair cover')
    return {'dot_interval': list(map(str, INNER_DOT)), 'midpoint_squared_strict_upper': str(INNER_MIDPOINT_SQUARED),
            'deficit_strict_upper': str(budget), 'nodes': sum(counts.values()), 'counts': counts,
            'maximum_depth': depth, 'unresolved': 0}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for base, field in ((HERE, 'files'), (ROOT, 'dependencies')):
        for path, digest in manifest[field].items():
            need(hashlib.sha256((base/path).read_bytes()).hexdigest()==digest, 'fixed '+path)
    upper = ROOT/'ten-point-ward-upper'
    ward = load('n10_ward_upper', upper/'verify_all.py')
    row = load('n10_row', ROOT/'ten-point-ward-radii/verify_all.py')
    needed = json.loads((ROOT/'ten-point-ward-inner-count/verification_report.json').read_bytes())
    need(needed['passed'] and needed['exact_inner_count']==2 and needed['outer_count']==8, 'exact prior population')
    data = json.loads((upper/'kernel.json').read_bytes())
    matrices, force = ward.corrections(data)
    need(bool(force), 'distance times the original two-point function')
    poly, shape = ward.polynomial(data, matrices, force)
    pack = ward.Packed(shape, True)
    ward.old.selftest()
    gap = 9*FLOOR-36*ward.C
    need(gap==F(9, 20000) and ward.log.log_lower(ward.witness())>ward.hx(data['B'])-F(1, 10**6),
         'marked row contradicts attainable witness')
    m = ward.DMIN-LOW-DOMAINS[0][1]
    need(m>0 and m*m>4*LOW*DOMAINS[0][1], 'inner pair too short if one radius is at most LOW')
    groups = []
    for index, (lo, hi, marked) in enumerate(EXCLUSIONS):
        slabs = json.loads((HERE/f'{index}.json').read_bytes())
        intervals = [tuple(map(F, slab['interval'])) for slab in slabs]
        need(bool(intervals) and intervals[0][0]==lo and intervals[-1][1]==hi, 'closed exclusion endpoints')
        need(all(a<b for a,b in intervals) and all(b==c for (_,b),(c,_) in zip(intervals, intervals[1:])),
             'consecutive marked slabs with no gaps')
        total = dict.fromkeys('AP012', 0)
        depth = 0
        for interval, slab in zip(intervals, slabs):
            coefficients = tuple(map(F, slab['coefficients']))
            need(len(slab['trees'])==len(DOMAINS), 'every neighbor radius band covered')
            for other, (bounds, tree) in enumerate(zip(DOMAINS, slab['trees'])):
                expression = counted_polynomial(row, poly, coefficients, marked, other)
                box = (interval, bounds, (-ward.U_NEW, ward.U_NEW))
                counts, top = verify_tree(ward, row, pack, expression, shape, box, tree)
                depth = max(depth, top)
                for k,v in counts.items():
                    total[k] += v
        groups.append({'excluded_closed_interval': [str(lo), str(hi)], 'slabs': len(slabs),
                       'nodes': sum(total.values()), 'counts': total, 'maximum_depth': depth, 'unresolved': 0})
    result = {'passed': True, 'mass': 9, 'applies_to_every_global_maximizer': True,
              'marked_row_lower': str(9*FLOOR), 'deficit_lower': str(gap),
              'small_radius_exclusion': ['0', str(LOW)], 'row_exclusions': groups,
              'inner_pair': verify_inner_pair(ward, row, pack, poly, shape)}
    output = json.dumps(result, indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_text(output, encoding='utf-8', newline='\n')
    print(output, end='')


if __name__ == '__main__':
    main()
