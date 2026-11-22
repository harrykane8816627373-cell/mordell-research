"""Replay the mass-nine upper bound using the cyclic separation of outer pairs."""

import argparse
import importlib.util
import json
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BANDS = ((F(37, 500), F(43, 250)), (F(169, 200), F(164, 125)))
CLASSES = ((0, 0, 0), (0, 1, 0), (1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 4))
WEIGHTS = (1, 16, 8, 8, 8, 4)
COSINES = ((F(-146, 1000), F(106, 125)), (F(-648, 1000), F(439, 1000)),
           (F(-954, 1000), F(-104, 1000)), (F(-1), F(-31, 32)))


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


curved = load('n10_curved', 'ten-point-curved-kernel/verify_all.py')
upper = curved.upper
need = upper.need


def read(path):
    return json.loads(path.read_bytes())


def check_domains():
    manifest = read(HERE / 'source_manifest.json')
    for base, group in ((HERE, 'files'), (ROOT, 'dependencies')):
        for name, expected in manifest[group].items():
            need(sha256((base / name).read_bytes()).hexdigest() == expected,
                 'fixed input ' + name)
    curved.verify_dependencies()
    radii = read(ROOT / 'ten-point-refined-radii/verification_report.json')
    need(radii['passed'] and radii['unresolved'] == 0 and radii['applies_to_every_global_maximizer']
         and radii['normalization_mass'] == 9 and (radii['exact_inner_count'], radii['exact_outer_count']) == (2, 8)
         and tuple(tuple(map(F, radii[k])) for k in ('inner_squared_radius_open_interval',
                                                     'outer_squared_radius_open_interval')) == BANDS
         and F(radii['midpoint']['inner_midpoint_squared_strict_upper']) == F(239, 2000)**2,
         'accepted radial and midpoint conditions')
    inner = read(ROOT / 'ten-point-conditioned-inner/verification_report.json')
    need(inner['passed'] and inner['unresolved'] == 0 and inner['applies_to_every_global_maximizer']
         and inner['normalization_mass'] == 9
         and tuple(map(F, inner['inner_squared_distance_open_interval'])) == (F(41, 100), F(11, 20))
         and F(inner['inner_angle_cosine_strict_upper']) == F(-77, 100), 'accepted inner-pair conditions')
    pairing = load('n10_tight_pairing', 'ten-point-tight-pairing/verify_all.py')
    kernel, bands, inner_midpoint = pairing.check_inputs()
    cheb = pairing.chebyshevs(max(pairing.HARMONICS))
    weights, _ = pairing.projection(kernel, bands, pairing.read(pairing.HERE / 'projection.json'), cheb)
    geometry = pairing.geometry(weights, cheb, bands[1], inner_midpoint)
    ca, cd = F(106, 125), F(31, 32)
    need(pairing.ALPHA_COS == ca and F(geometry['cyclic_rank_pairing_cosine_strict_upper']) == -cd
         and F(geometry['paired_squared_distance_strict_lower']) == F(10647, 3200),
         'replayed cyclic pair separation')

    scale = 10**15
    def sqrt_upper(x):
        return F(isqrt((scale**2 * x.numerator) // x.denominator) + 1, scale)

    sa, sd = sqrt_upper(1-ca**2), sqrt_upper(1-cd**2)
    sine_multipliers = (1, 2*ca, 4*ca**2-1)
    cosines = (F(1), ca, 2*ca**2-1, 4*ca**3-3*ca)
    need(F(1, 2)<ca<cd and all(v>0 for v in sine_multipliers),
         '0 < delta < alpha < pi/3 and positive sine multipliers')
    for k, (lo, hi) in enumerate(COSINES[:3], 1):
        m = 4-k
        upper_cosine = cosines[m]*cd + sa*sd*sine_multipliers[m-1]
        need(lo < -upper_cosine and cosines[k] <= hi, 'rational angular enclosure')


def outside(i, j, k, box):
    (ul, uh), (vl, vh), (tl, th) = box
    if curved.prune(box, i == j == 0):
        return True
    if i == j == 0:
        return (uh+vh-2*tl < F(41, 100) or ul+vl-2*th > F(11, 20)
                or tl**2 < F(77, 100)**2*ul*vl or ul+vl+2*tl > F(239, 1000)**2)
    if i != j:
        return False
    lo, hi = COSINES[k-1]
    if th < 0 and th**2 > lo**2*uh*vh:
        return True
    if hi < 0:
        if tl >= 0 or tl**2 < hi**2*ul*vl:
            return True
    elif tl > 0 and tl**2 > hi**2*uh*vh:
        return True
    return k == 4 and uh+vh-2*tl < F(10647, 3200)


def root_boxes(i, j, k):
    if i == j == 0:
        return [(BANDS[i], BANDS[j], (F(-27, 200), F(-1, 10)))]
    if i != j:
        end = F(12, 25)
        need(end**2 > BANDS[i][1]*BANDS[j][1], 'mixed dot product bound')
        return [(BANDS[i], BANDS[j], (-end, F(0))), (BANDS[i], BANDS[j], (F(0), end))]
    lo, hi = COSINES[k-1]
    left, right = BANDS[1][1]*lo, BANDS[1][int(hi >= 0)]*hi
    intervals = ((left, F(0)), (F(0), right)) if right > 0 else ((left, right),)
    return [(BANDS[i], BANDS[j], interval) for interval in intervals]


def check_tree(tree, box, poly, shape, packer, i, j, k, c):
    stack = [(box, upper.initial(poly, shape, box), 0)]
    counts = dict(P=0, A=0, Q=0, S=0)
    maximum = 0
    for tag in tree:
        need(bool(stack), 'a tree node has a remaining box')
        box, packed, depth = stack.pop()
        maximum = max(maximum, depth)
        if tag in '012':
            stack.extend(curved.split(packer, packed, box, depth, int(tag)))
            counts['S'] += 1
            continue
        accepted = (tag == 'P' and outside(i, j, k, box)
                    or tag == 'A' and packer.dominates(packed, box, c)
                    or tag == 'Q' and packer.curved_dominates(packed, box, c))
        need(accepted, f'valid {tag} leaf, class {(i, j, k)}, box {box}')
        counts[tag] += 1
    need(not stack and counts['P'] + counts['A'] + counts['Q'] == counts['S'] + 1,
         'complete binary coverage')
    return counts, maximum


def verify(progress=False):
    check_domains()
    data = read(curved.HERE / 'kernel.json')
    need((data['n'], data['mass'], len(data['Q_blocks'])) == (10, 9, 10), 'fixed kernel normalization')
    for matrix in upper.gram.decode(data['Q_blocks']):
        upper.gram.positive_ldl(matrix)
    ward, force = upper.corrections(data)
    need(bool(force), 'kernel is multiplied by the positive squared distance')
    base, shape = upper.polynomial(data, ward, force)
    upper.crosscheck(data, ward, force, base)
    params = curved.parsed(data)
    packer = curved.Packed(shape)
    upper.old.selftest()
    curved.curvature_selftest()

    certificates = read(HERE / 'cover.json')
    need(len(certificates) == len(CLASSES) and sum(WEIGHTS) == 45, 'all six pair classes')
    rows = []
    budget = upper.hx(data['B'])
    for (i,j,k), weight, cert in zip(CLASSES, WEIGHTS, certificates):
        c = F(cert['c'])
        need(c > 0, 'positive pair allowance')
        budget += weight*c
        poly = curved.root_polynomial(base, shape, params, i, j)
        curved.check_root(base, poly, params, i, j)
        boxes = root_boxes(i,j,k)
        need(len(boxes) == len(cert['trees']), 'one tree per root box')
        visited, maximum = 0, 0
        for box, tree in zip(boxes, cert['trees']):
            counts, depth = check_tree(tree, box, poly, shape, packer, i, j, k, c)
            visited += sum(counts.values())
            maximum = max(maximum, depth)
        rows.append({'class': [i,j,k], 'weight': weight, 'deficit_floor': str(-c),
                     'visited': visited, 'maximum_depth': maximum})
        if progress:
            print('accepted', (i,j,k), visited, flush=True)

    candidate_lower = upper.actual_root_lower()
    ratio = F(5011, 5000)
    margin = upper.log.log_lower(ratio*candidate_lower) - budget
    need(margin > F(1, 1_000_000), 'strict upper comparison against the true candidate')
    return {'passed': True, 'normalization_mass': 9, 'global_over_candidate_strict_upper': str(ratio),
            'log_comparison_margin_strict_lower': '1/1000000', 'log_upper_mass_nine': str(budget),
            'visited': sum(row['visited'] for row in rows), 'unresolved': 0, 'classes': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--progress', action='store_true')
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    report = verify(args.progress)
    path = HERE / 'verification_report.json'
    result = json.dumps(report, indent=2) + '\n'
    if args.record_report:
        path.write_text(result, encoding='utf-8', newline='\n')
    else:
        need(report == read(path), 'matches the saved report')
    print(result, end='')


if __name__ == '__main__':
    main()
