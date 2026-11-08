"""Accept fixed nonlinear supports, then independently exercise exact algebra."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json
import random
import sys

from core import (PI_LO, PI_HI, bernstein_integer_coefficients, check_certificate,
                  check_pi, cycle_geometry, edge_minorant, sine_minorant,
                  subdivide, validate_geometry)


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def power_value(poly, x):
    return sum((coefficient*x**i for i, coefficient in enumerate(poly)), F(0))


def bernstein_value(coeff, x):
    n = len(coeff)-1
    return sum((a*comb(n, k)*x**k*(1-x)**(n-k) for k, a in enumerate(coeff)), F(0))


def logarithm_lower(value, terms=100):
    need(value >= 1, 'positive log input')
    t = (value-1)/(value+1)
    return 2*sum((t**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))


def direct_square_value(radii, rotations):
    points = []
    for r, (x, y) in zip(radii, rotations):
        need(x*x+y*y == 1, 'rational unit direction')
        points.extend([(r*x, r*y), (-r*y, r*x), (-r*x, -r*y), (r*y, -r*x)])
    need(sum(x for x, y in points) == sum(y for x, y in points) == 0, 'direct point centroid')
    mass = sum((x*x+y*y for x, y in points), F(0))
    delta = F(1)
    for (x, y), (a, b) in combinations(points, 2):
        delta *= (x-a)**2+(y-b)**2
    n = len(points)
    return (F(n)/mass)**(n*(n-1)//2)*delta


def input_paths(base):
    return sorted(p for p in base.rglob('*') if p.is_file() and
                  (p.suffix in ('.py', '.md') or p.name == '.gitattributes' or
                   (p.parent.name == 'certificates' and p.suffix == '.json')))


def verify(base):
    manifest = json.loads((base/'source_manifest.json').read_text(encoding='utf-8'))
    need(set(manifest) == {'algorithm', 'files'} and manifest['algorithm'] == 'sha256', 'manifest schema')
    hashes = {p.relative_to(base).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths(base)}
    need(manifest['files'] == hashes, 'fixed input added, missing or modified')
    names = ['four_squares.json', 'four_squares_box.json', 'mixed_twentyfive_cycle.json']
    certs, results = {}, {}
    for name in names:
        raw = (base/'certificates'/name).read_bytes()
        need(len(raw) < 100000, 'bounded fixed certificate')
        certs[name] = json.loads(raw)
        results[name] = check_certificate(certs[name])
    sq, box, mix = [results[name] for name in names]
    need(sq['nonlinear_gap'] == F(119489373, 20000000), 'fixed square nonlinear exact gap')
    need(box['nonlinear_gap'] == F(119455433, 20000000), 'radius-box nonlinear exact gap')
    need(mix['nonlinear_gap'] == F(795833, 5000000), 'mixed nonlinear exact gap')
    for result in results.values():
        need(result['nonlinear_gap'] > result['jensen_gap'] > result['old_gap'], 'both reinforcement levels improve')
    need(sq['nonlinear_gap'] > F(5, 2)*sq['old_gap'], 'square gap greater than 2.5 times baseline')
    need(mix['nonlinear_gap'] > F(12, 5)*mix['old_gap'], 'mixed gap greater than 2.4 times baseline')
    need(sq['upper_ratio_to_old'] < F(483, 1000), 'square upper less than 48.3 percent of old')
    need(mix['upper_ratio_to_old'] < F(919, 1000), 'mixed upper less than 91.9 percent of old')

    # An independent polynomial/basis evaluator checks the exact identity
    # and both half-interval parameterizations, not merely coefficient signs.
    rng = random.Random(11082026)
    identities = 0
    for degree in range(1, 25):
        poly = [F(rng.randrange(-19, 20), rng.randrange(1, 10)) for _ in range(degree+1)]
        root, D = bernstein_integer_coefficients(poly)
        left, right = subdivide(root)
        scale = D*(1 << degree)
        for x in (F(0), F(1, 7), F(1, 3), F(3, 5), F(1)):
            need(power_value(poly, x) == bernstein_value(root, x)/D, 'independent Bernstein basis identity')
            need(power_value(poly, x/2) == bernstein_value(left, x)/scale, 'left subdivision parameterization')
            need(power_value(poly, (1+x)/2) == bernstein_value(right, x)/scale, 'right subdivision parameterization')
            identities += 3
    # Exact special angles have rational sine-squared values; log lower
    # bounds come from a separate positive atanh series.
    special_angles = [(F(0), F(0)), (F(1, 3), F(1, 4)), (F(1, 2), F(1, 2)),
                      (F(2, 3), F(3, 4)), (F(1), F(1))]
    trigonometric_checks = 0
    for cert in certs.values():
        ms, rs = validate_geometry(cert['populations'], cert['radius_intervals'], cert['center'])
        for cycle in cert['cycles']:
            for edge in cycle_geometry(ms, rs, cycle['vertices']):
                poly = edge_minorant(edge['g'], edge['eta'])
                for x, sine_squared in special_angles:
                    need(0 <= power_value(sine_minorant(), x) <= sine_squared, 'sine minorant at exact angle')
                    lower = edge['g']*logarithm_lower(1/(1-edge['eta']*sine_squared))
                    need(power_value(poly, x) <= lower, 'independent logarithm lower vs polynomial')
                    trigonometric_checks += 1
    check_pi()

    # Full-coordinate checks are implementation checks, not coverage proofs.
    dirs = [(F(1), F(0)), (F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)),
            (F(-8, 17), F(15, 17)), (F(7, 25), F(-24, 25))]
    endpoints = [[F(q) for q in pair] for pair in certs['four_squares_box.json']['radius_intervals']]
    direct_count = 0
    for radii in product(*endpoints):
        for shift in range(5):
            rotations = [dirs[(i+shift) % 5] for i in range(4)]
            actual = direct_square_value(radii, rotations)
            need(actual <= F(certs['four_squares_box.json']['all_phase_upper']), 'direct sixteen-point value exceeds box upper')
            direct_count += 1

    mutations = []
    original = certs['four_squares.json']
    def mutate(label, fn):
        c = deepcopy(original)
        fn(c)
        mutations.append((label, c))
    mutate('zeroed support intercept', lambda c: c['cycles'][0]['supports'][0].__setitem__('intercept', '0'))
    mutate('negative support intercept', lambda c: c['cycles'][0]['supports'][0].__setitem__('intercept', '-1'))
    mutate('forged dual slope', lambda c: c['cycles'][0].__setitem__('slope', '128'))
    mutate('negative dual slope', lambda c: c['cycles'][0].__setitem__('slope', '-1'))
    mutate('over-spent packing edge', lambda c: c['cycles'][0].__setitem__('weight', '1'))
    mutate('zero packing weight', lambda c: c['cycles'][0].__setitem__('weight', '0'))
    mutate('missing Bernstein leaf', lambda c: c['cycles'][0]['supports'][0]['cover'].pop())
    mutate('duplicate Bernstein leaf', lambda c: c['cycles'][0]['supports'][0]['cover'].append(c['cycles'][0]['supports'][0]['cover'][0]))
    mutate('hidden unverified subtree', lambda c: c['cycles'][0]['supports'][0]['cover'].append(''))
    mutate('invalid dyadic path', lambda c: c['cycles'][0]['supports'][0]['cover'].__setitem__(0, '02'))
    mutate('unbound edge index', lambda c: c['cycles'][0]['supports'][0].__setitem__('edge', [0, 3]))
    mutate('even-multiplier-sum cycle', lambda c: c['cycles'][0].__setitem__('vertices', [0, 1, 2, 3]))
    mutate('moved radius endpoint', lambda c: c['radius_intervals'][0].__setitem__(0, '8999/10000'))
    mutate('floating endpoint', lambda c: c['radius_intervals'][0].__setitem__(0, 0.9))
    mutate('altered reported gap', lambda c: c.__setitem__('nonlinear_gap', '100'))
    mutate('altered final upper', lambda c: c.__setitem__('all_phase_upper', '1'))
    mutate('scope escalation', lambda c: c.__setitem__('scope', 'unrestricted M_n'))
    mutate('unrecognized claim', lambda c: c.__setitem__('certified_global_optimum', True))
    mutate('undefined rational', lambda c: c['cycles'][0].__setitem__('slope', '1/0'))
    mutate('empty cycle list', lambda c: c.__setitem__('cycles', []))
    rejected = []
    for label, cert in mutations:
        try:
            check_certificate(cert)
        except ValueError:
            rejected.append(label)
        else:
            raise ArithmeticError('tampered certificate accepted: '+label)

    comparisons = {}
    previous = {'four_squares.json': 'four_squares_triangles.json',
                'four_squares_box.json': 'four_squares_radius_box.json',
                'mixed_twentyfive_cycle.json': 'mixed_twentyfive_cycles.json'}
    for name, old_name in previous.items():
        p = base.parent/'phase-packing'/'certificates'/old_name
        if p.exists():
            old = json.loads(p.read_text(encoding='utf-8'))
            need(F(old['optimal_log_gap']) == results[name]['old_gap'], 'old LP optimum comparison')
            need(F(old['pair_upper']) == F(certs[name]['pair_box_upper']), 'same normalization and family upper')
            comparisons[old_name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return {'passed': True, 'python': sys.version.split()[0], 'optimization_safe': True,
            'fixed_nonlinear_certificates': len(certs),
            'complete_bernstein_leaf_count': sum(r['bernstein_leaves'] for r in results.values()),
            'independent_polynomial_identity_checks': identities,
            'exact_special_angle_and_log_checks': trigonometric_checks,
            'exact_direct_sixteen_point_coordinate_values': direct_count,
            'distances_per_direct_value': 120,
            'tamper_cases_rejected': rejected,
            'results': {name: {k: str(v) if isinstance(v, F) else v for k, v in result.items()} for name, result in results.items()},
            'coverage': 'entire angular interval via Bernstein covers and entire stated radius boxes via analytic endpoints',
            'scope': 'fixed regular-ring populations and stated radius ranges; no unrestricted M_n or exact angular optimum claim',
            'previous_certificate_sha256_if_available': comparisons,
            'source_sha256': hashes}


if __name__ == '__main__':
    base = Path(__file__).resolve().parent
    report = verify(base)
    (base/'verification_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('passed', 'optimization_safe', 'fixed_nonlinear_certificates',
                                          'complete_bernstein_leaf_count', 'independent_polynomial_identity_checks',
                                          'exact_special_angle_and_log_checks', 'exact_direct_sixteen_point_coordinate_values')}, indent=2))
    print('Rejected', len(report['tamper_cases_rejected']), 'tampered inputs; exact results in verification_report.json')
