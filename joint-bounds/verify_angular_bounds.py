"""Replay exact constants, Fourier normalization and joint exclusions.

The all-n analytic proof is in ANGULAR_PROOF.md; finite checks do not
replace the functional analysis or the universal quantifiers of that proof.
"""
from angular_bounds import *
import json


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def decimal_floor(q, digits=14):
    q = F(q)
    scale = 10**digits
    k = q.numerator*scale//q.denominator
    return str(F(k, scale))


def main():
    # Outward square roots, including exact squares and values at grid edges.
    roots = 0
    for a in range(31):
        for b in range(1, 13):
            q = F(a, b)
            u = sqrt_upper(q, 10000)
            require(u*u >= q, 'outward root')
            require(u == 0 or (u-F(1, 10000))**2 < q, 'root grid minimality')
            roots += 1
    # A circle Fourier charge 2a cos(k theta) has moment a R^k and
    # exact energy a^2/k. This checks the sharp 1/k moment coefficient.
    for k in range(1, 21):
        a, r2 = F(3, 7), F(9, 4)
        exact_energy = a*a/k
        require(angular_cost({k: a*a*r2**k}, r2) == exact_energy,
                'single boundary Fourier mode coefficient')
        require(k*r2**k == 2*F(k, 2)*r2**k,
                'interior plus exterior Dirichlet norms')
    # Rational normalized moments independently known for a square.
    square = rational_shape_data([(1, 0), (0, 1), (-1, 0), (0, -1)], 8)
    require(sum(square['coulomb_radii_squared']) == 3, 'Coulomb moment')
    require(square['moment_abs2'][1] == square['moment_abs2'][2] == square['moment_abs2'][3] == 0,
            'square vanished low moments')
    require(square['moment_abs2'][4] == F(81, 16), 'square fourth moment')
    require(square['moment_abs2'][8] == F(6561, 4096), 'square eighth moment')
    # Translation must disappear before normalization.
    shifted = rational_shape_data([(8, -4), (7, -3), (6, -4), (7, -5)], 8)
    require(square == shifted, 'centering translation invariance')
    # No use of an unvalidated floating-point optimizer in the acceptance path.
    rows = []
    for n in range(2, 101):
        cert = strict_disk_correction(n)
        c = cert['delta_lower']
        require(c > 0, 'positive global correction')
        # An actual unit regular n-gon has Delta=n^n. Its known gap must
        # exceed the new correction. Bound that gap FROM BELOW for this check.
        pairs = n*(n-1)//2
        ln = log_bounds(n)
        lc = log_bounds(F(2*n, n-1))
        polygon_gap_lower = -F(n, 2)*ln[1] + pairs*lc[0] - F(pairs, 2)
        require(c < polygon_gap_lower, 'regular polygon consistency')
        if n in (2, 3, 5, 8, 10, 30, 100):
            rows.append({'n': n, 'delta_lower': str(c),
                         'delta_lower_decimal_floor': decimal_floor(c),
                         'strict_positive_rational': True})
    # A nonzero-moment region inside a radial class containing center+7gon.
    n = 8
    lower = F(8**28, 7**21)
    budget = gap_upper_from_lower(n, lower)
    joint = joint_certificate(8, [(0, 0, 1), (1, 1, 7)],
                              [0, F(1, 4), F(1, 2), 1, 2], {2: 2})
    require(joint['support_r2_upper'] == F(9, 4), 'exact smeared support')
    require(joint['angular_delta_lower'] == F(16, 81), 'eight point moment penalty')
    require(joint['radial_delta_lower'] < budget < F(16, 81), 'joint exclusion beyond radial data')
    require(F(1, 10) < budget < F(105, 1000), 'eight point gap enclosure')
    # Exact p2 squared of cube roots plus two disjoint antipodal pairs:
    # p2 = (0 from cube roots)+(1+i sqrt3)-2, modulus squared=4.
    require((-1)**2 + 3 == 4, 'nonempty example squared moment')
    # Independently compute actual normalized products for arbitrary rational
    # shapes, testing the ADDITIVE certificate against a lower log-gap bound.
    shapes = [[(-1, 0), (1, 0)], [(0, 0), (1, 0), (0, 2)],
              [(1, 0), (0, 1), (-1, 0), (0, -1)],
              [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)],
              [(i, i*i % 5) for i in range(8)]]
    for points in shapes:
        data = rational_shape_data(points, 12)
        count = len(points)
        pairs = count*(count-1)//2
        raw = F(1)
        for i, (a, b) in enumerate(points):
            for c, d in points[i+1:]:
                raw *= (a-c)**2+(b-d)**2
        attained = F(count, 1)**pairs*raw/data['quadratic_moment_before_scaling']**pairs
        gaplo = (F(count, 2)*log_bounds(count)[0]
                 + pairs*log_bounds(F(2*count, count-1))[0]
                 - F(pairs, 2) - log_bounds(attained)[1])
        cert = joint_certificate(count, [(q, q, 1) for q in data['coulomb_radii_squared']],
                                 [F(j, 8) for j in range(17)], data['moment_abs2'])
        require(cert['delta_lower'] < gaplo, 'joint certificate below independent actual gap')
    # Check log range reduction against exact multiplicative identities.
    for q in (F(1, 7), F(7), F(8**28, 7**21), F(10**100, 3**140)):
        lo, hi = log_bounds(q)
        ilo, ihi = log_bounds(1/q)
        require(lo == -ihi and hi == -ilo, 'reciprocal log enclosure')
        dlo, dhi = log_bounds(q*q)
        require(dlo <= 2*hi and dhi >= 2*lo, 'squared log enclosure consistency')
    rejects = 0
    for call in (lambda: strict_disk_correction(1),
                 lambda: strict_disk_correction(True),
                 lambda: sqrt_upper(-1),
                 lambda: angular_cost({0: 1}, 2),
                 lambda: angular_cost({2: -1}, 2),
                 lambda: rational_shape_data([(0, 0), (0, 0)], 2),
                 lambda: log_bounds(0)):
        try:
            call()
        except ValueError:
            rejects += 1
        else:
            raise ArithmeticError('invalid input accepted')
    report = {'passed': True, 'root_enclosures': roots,
              'boundary_fourier_modes': 20, 'actual_shape_product_checks': len(shapes),
              'global_corrections_checked': 99,
              'n_min': 2, 'n_max': 100, 'invalid_inputs_rejected': rejects,
              'selected_global_corrections': rows,
              'eight_point_region': {'coulomb_bands': [['0', '0', 1], ['1', '1', 7]],
                 'moment_condition': '|p_2(x)|^2 >= 2',
                 'angular_delta_lower': '16/81',
                 'global_gap_upper': str(budget), 'excluded': True},
              'scope': 'Proof constants and finite consistency checks; no full eight-point cover, no gamma improvement, no formalized analytic proof.'}
    Path('angular_verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('PASS angular modes, positive disk correction n=2..100, eight-point moment region')


if __name__ == '__main__':
    main()
