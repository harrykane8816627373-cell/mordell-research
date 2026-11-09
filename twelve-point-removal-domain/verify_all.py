"""Exact endpoint and witness checks for the mass-eleven twelve-point domain."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
from pathlib import Path
import argparse
import hashlib
import json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RADIUS = F(81, 50)
COMPLEMENT = F(7573, 1000)


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def exp_lower(x):
    term = total = F(1)
    for j in range(1, 101):
        term = term*x/j
        total += term
    return total


SQRT11_UPPER = F(33166248, 10**7)
C11 = 11**5*SQRT11_UPPER*F(1, 5)**55/exp_lower(F(55, 2))
C10 = 10**5*F(2, 9)**45/exp_lower(F(45, 2))


def old_witness():
    raw = json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())
    points = next(row['integer_points'] for row in raw['rows'] if row['n'] == 12)
    need(len(points) == 12, 'complete twelve-point witness')
    sx = sum(x for x, y in points)
    sy = sum(y for x, y in points)
    D = 12*sum(x*x+y*y for x, y in points)-sx*sx-sy*sy
    P = prod((x-a)**2+(y-b)**2 for (x, y), (a, b) in combinations(points, 2))
    need(D > 0 and P > 0, 'nondegenerate witness')
    centered = [(12*x-sx, 12*y-sy) for x, y in points]
    scale_squared = F(11, 12*D)
    need(scale_squared*sum(x*x+y*y for x, y in centered) == 11, 'exact witness mass')
    return F(132**66*P, D**66), centered, scale_squared


def schwarz_factor(delta):
    need(0 <= delta < 2, 'finite second-order Schwarz branch')
    return ((delta+2)/(2-delta))**2


def single_upper(u):
    q = 11-F(12, 11)*u
    return C11*q**55*(F(144, 121)*u)**11*schwarz_factor(6*q-55)


def uniform_pair_upper(q):
    return C10*q**45*(F(4, 7)*(11-q))**21*schwarz_factor(6*q-45)


def pair_terms(u, v, t):
    q = 11-F(11, 10)*(u+v)-t/5
    d = u+v-2*t
    b = (121*u+v+22*t)*(u+121*v+22*t)/10000
    return q, d, b


def pair_upper(u, v, t):
    q, d, b = pair_terms(u, v, t)
    return C10*q**45*d*b**10*schwarz_factor(6*q-45)


def pair_box_excluded(box, L):
    """A sufficient exact exclusion test for rational closed three-variable boxes."""
    (ul, uh), (vl, vh), (tl, th) = box
    need(0 <= ul <= uh and 0 <= vl <= vh and tl <= th, 'ordered physical radial box')
    if ul >= RADIUS or vl >= RADIUS:
        return True
    uh, vh = min(uh, RADIUS), min(vh, RADIUS)
    qhi = 11-F(11, 10)*(ul+vl)-tl/5
    if qhi <= COMPLEMENT:
        return True
    delta = 6*qhi-45
    if delta >= 2:
        return False
    dhi = uh+vh-2*tl
    rhi = (121*uh+vh+22*th)/100
    shi = (uh+121*vh+22*th)/100
    if dhi <= F(1, 3) or min(rhi, shi) < 0:
        return True
    upper = C10*qhi**45*dhi*(rhi*shi)**10*schwarz_factor(delta)
    return upper < L


def check_derivative_identity():
    # Homogeneity reduces d[s(6R-s)^20]/ds to R=1.
    left = [(-1)**i*(i+1)*comb(20, i)*6**(20-i) for i in range(21)]
    factor = [(-1)**i*comb(19, i)*6**(19-i) for i in range(20)]
    right = [0]*21
    for i, c in enumerate(factor):
        right[i] += 6*c
        right[i+1] -= 21*c
    need(left == right, 'complete derivative factorization')


def checked_manifest():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for name, sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == sha, 'fixed own input '+name)
    for name, sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha, 'fixed prior theorem/input '+name)
    return manifest


def single_radius_claim(L):
    need(SQRT11_UPPER**2 > 11, 'strict rational square-root upper bound')
    need(RADIUS < F(121, 72), 'radius strictly improved over the old endpoint')
    need(0 < 6*(11-F(12, 11)*RADIUS)-55 < 2, 'Schwarz range on the excluded radius interval')
    need(F(11)/RADIUS-F(144, 11) < 0, 'one-deletion function decreases')
    ratio = single_upper(RADIUS)/L
    need(ratio < F(997, 1000), 'strict one-deletion endpoint exclusion')
    return {'passed': True, 'strict_radius_upper': str(RADIUS),
            'old_radius_upper': '121/72',
            'endpoint_upper_over_witness_strict_upper': '997/1000',
            'endpoint_upper_over_witness_display': float(ratio)}


def verify_single_radius(V11):
    checked_manifest()
    L, _, _ = old_witness()
    need(V11 == L, 'caller witness equals the fixed feasible mass-eleven value')
    return single_radius_claim(L)


def verify():
    manifest = checked_manifest()
    L, points, scale_squared = old_witness()
    radius_claim = single_radius_claim(L)
    need(F(15, 2) < COMPLEMENT and 0 < 6*COMPLEMENT-45 < 2,
         'strictly improved complement endpoint in the finite Schwarz range')
    need(F(45)/COMPLEMENT-F(21)/(11-COMPLEMENT)+12 > 0,
         'two-deletion uniform function increases')
    check_derivative_identity()
    second_ratio = uniform_pair_upper(COMPLEMENT)/L
    need(second_ratio < F(992, 1000), 'strict two-deletion endpoint exclusion')
    qmin = F(11)
    max_radius = F(0)
    active = 0
    for (x, y), (a, b) in combinations(points, 2):
        u = scale_squared*(x*x+y*y)
        v = scale_squared*(a*a+b*b)
        t = scale_squared*(x*a+y*b)
        q, d, _ = pair_terms(u, v, t)
        need(max(u, v) < RADIUS and q > COMPLEMENT and d > F(1, 3) and t*t <= u*v,
             'every saved witness pair retained by the linear domain')
        max_radius = max(max_radius, u, v)
        qmin = min(qmin, q)
        if 6*q-45 < 2:
            active += 1
            need(pair_upper(u, v, t) >= L, 'active nonlinear cut retains saved witness pair')
    example = F(8, 5), F(8, 5), -F(1, 2)
    u, v, t = example
    q, d, _ = pair_terms(*example)
    need(u < RADIUS and v < RADIUS and q > COMPLEMENT and d > F(1, 3) and t*t < u*v,
         'example passes every stated simpler pair cut')
    need(pair_upper(*example)/L < F(3, 4), 'example violates nonlinear cut')
    eps = F(1, 1000)
    box = tuple((x-eps, x+eps) for x in example)
    (ul, uh), (vl, vh), (tl, th) = box
    need(max(uh, vh) < RADIUS and max(abs(tl), abs(th))**2 < ul*vl and ul+vl-2*th > F(1, 3),
         'entire example box meets radius, physicality and separation cuts')
    need(11-F(11, 10)*(uh+vh)-th/5 > COMPLEMENT, 'entire example box meets complement cut')
    need(pair_box_excluded(box, L), 'entire nonempty example box excluded by nonlinear cut')
    return {'passed': True, 'normalization': 'sum |z_i|^2=11', 'global_maximizers_only': True,
            'squared_radius_strict_upper': str(RADIUS), 'ten_point_centered_mass_strict_lower': str(COMPLEMENT),
            'linear_pair_cut': '11*(u+v)+2*t<3427/100',
            'endpoint_upper_over_witness_strict_upper': ['997/1000', '992/1000'],
            'endpoint_upper_over_witness_display': [radius_claim['endpoint_upper_over_witness_display'], float(second_ratio)],
            'single_radius': radius_claim,
            'witness_pairs_retained': 66, 'witness_pairs_in_active_nonlinear_branch': active,
            'witness_max_squared_radius_display': float(max_radius), 'witness_min_ten_mass_display': float(qmin),
            'nonlinear_excluded_pair': list(map(str, example)), 'whole_excluded_box_radius': str(eps),
            'sharp_global_value_proved': False,
            'input_sha256': manifest['files'], 'dependency_sha256': manifest['dependencies']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    report = verify()
    target = HERE/'verification_report.json'
    if args.record_report:
        target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    else:
        need(json.loads(target.read_bytes()) == report, 'fresh exact report equals the saved record')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
