"""Independent exact checks for the phase-frustration results; -O safe."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm, prod
from pathlib import Path
import hashlib
import json
import sys

from phase_frustration import (bad_triangles, encode, examples, family_bound,
                               maximizing_signs, pairwise_relaxation,
                               sign_product, simultaneous_pair_maxima,
                               three_equal_ring_value, triangle_gap)


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def determinant(matrix):
    """Independent dense rational elimination, not the binomial formula."""
    a = [list(map(F, row)) for row in matrix]
    d = F(1)
    for k in range(len(a)):
        pivot = next((j for j in range(k, len(a)) if a[j][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            d = -d
        p = a[k][k]
        d *= p
        for j in range(k+1, len(a)):
            q = a[j][k]/p
            for t in range(k+1, len(a)):
                a[j][t] -= q*a[k][t]
            a[j][k] = F(0)
    return d


def binomial_resultant(m, A, k, B):
    f = [F(1)] + [F(0)]*(m-1) + [-F(A)]
    g = [F(1)] + [F(0)]*(k-1) + [-F(B)]
    size = m+k
    matrix = [[F(0)]*size for _ in range(size)]
    for row in range(k):
        matrix[row][row:row+m+1] = f
    for row in range(m):
        matrix[k+row][row:row+k+1] = g
    return determinant(matrix)


def pair_conditions(ms, signs):
    for i, j in combinations(range(len(ms)), 2):
        g = gcd(ms[i], ms[j])
        if signs[i]**(ms[j]//g) == signs[j]**(ms[i]//g):
            return False
    return True


def qa(x, y):
    return (x[0]+y[0], x[1]+y[1])


def qneg(x):
    return (-x[0], -x[1])


def qm(x, y):
    """Q(sqrt(3)), represented by a+b*sqrt(3), exactly."""
    return (x[0]*y[0]+3*x[1]*y[1], x[0]*y[1]+x[1]*y[0])


def scale(x, a):
    return (a*x[0], a*x[1])


def direct_twelve_value():
    radii = (F(4, 5), F(1), F(6, 5))
    directions = (((F(1), F(0)), (F(0), F(0))),
                  ((F(0), F(1, 2)), (F(1, 2), F(0))),
                  ((F(1, 2), F(0)), (F(0), F(1, 2))))
    points = []
    for r, (x, y) in zip(radii, directions):
        x, y = scale(x, r), scale(y, r)
        for _ in range(4):
            points.append((x, y))
            x, y = qneg(y), x
    Q = (F(0), F(0))
    center_x = (F(0), F(0))
    center_y = (F(0), F(0))
    value = (F(1), F(0))
    for x, y in points:
        Q = qa(Q, qa(qm(x, x), qm(y, y)))
        center_x, center_y = qa(center_x, x), qa(center_y, y)
    need(center_x == center_y == (0, 0), 'direct twelve-point centroid')
    need(Q[1] == 0 and Q[0] == F(308, 25), 'direct twelve-point moment')
    for (x, y), (u, v) in combinations(points, 2):
        dx, dy = qa(x, qneg(u)), qa(y, qneg(v))
        d2 = qa(qm(dx, dx), qm(dy, dy))
        value = qm(value, d2)
    need(value[1] == 0 and value[0] > 0, 'direct twelve-point rational product')
    return value[0]*(F(12)/Q[0])**66


def verify():
    resultant_count = 0
    for m in range(2, 10):
        for k in range(2, 10):
            for sa, sb in product((-1, 1), repeat=2):
                A, B = sa*F(2, 3), sb*F(7, 5)
                g = gcd(m, k)
                actual = binomial_resultant(m, A, k, B)
                expected_squared = (A**(k//g)-B**(m//g))**(2*g)
                need(actual*actual == expected_squared, 'Sylvester vs binomial identity')
                resultant_count += 1
    need(binomial_resultant(4, 1, 6, 1) == 0, 'resultant collision case')

    family_count = 0
    sign_assignment_count = 0
    compatible_count = 0
    smallest_totals = {}
    for K in range(1, 6):
        for ms in combinations_with_replacement(range(2, 13), K):
            explicit = [s for s in product((-1, 1), repeat=K) if pair_conditions(ms, s)]
            sign_assignment_count += 2**K
            wanted = simultaneous_pair_maxima(ms)
            need(bool(explicit) == wanted, 'criterion vs all sign assignments')
            # Check the integer cycle independently of valuations.
            odd_cycles = []
            for tri in combinations(range(K), 3):
                L = lcm(*(ms[i] for i in tri))
                a = tuple(L//lcm(ms[i], ms[j]) for i, j in combinations(tri, 2))
                if sum(a) % 2:
                    odd_cycles.append(tri)
            need(bool(odd_cycles) == (not wanted), 'triangle coverage of every incompatible family')
            need([t[0] for t in bad_triangles(ms)] == odd_cycles, 'odd cycle records')
            witness = maximizing_signs(ms)
            need((witness is not None) == wanted, 'constructive equality witness availability')
            if wanted:
                need(pair_conditions(ms, witness), 'constructive signs satisfy every pair')
                compatible_count += 1
                smallest_totals[K] = min(sum(ms), smallest_totals.get(K, sum(ms)))
                if K >= 3:
                    need(sum(ms) >= 2**(K-1)+4, 'logarithmic cardinality bound')
                # Distinct radii; actual exact normalized equality for this witness.
                rings = tuple((m, F(i+2, 3)) for i, m in enumerate(ms))
                need(sign_product(rings, witness) == pairwise_relaxation(rings),
                     'all-pairs product equality for constructive signs')
            family_count += 1
    need(smallest_totals == {1: 2, 2: 4, 3: 8, 4: 12, 5: 20}, 'minimum total populations')

    ex = examples()
    triple = ex['square_triple']
    need(triple['log_gap_lower'] == F(207360000, 216002897), 'square triple exact gap')
    cv = ex['continuous_vs_sign']
    need(cv['continuous_value'] == direct_twelve_value(), 'independent Q(sqrt(3)) twelve-point product')
    need(cv['continuous_value'] > 2*cv['best_sign_value'], 'strict improvement over all sign phases')
    need(cv['continuous_value'] < triple['all_continuous_phase_upper'], 'continuous witness respects strict family upper')
    need(cv['best_sign_value'] < triple['all_continuous_phase_upper'], 'best sign respects family upper')
    # Independent exact upper comparisons for sign witnesses on mixed families.
    fixed_phase_checks = 0
    for ms in ((2, 4, 4), (4, 4, 4), (2, 2, 2, 2, 2, 2), (3, 3, 2), (2, 3, 4, 6)):
        rings = tuple((m, F(i+4, 5)) for i, m in enumerate(ms))
        bound = family_bound(rings)
        spent = [edge for t in bound['selected_triangles'] for edge in t['edges']]
        need(len(set(spent)) == len(spent), 'no double-counted triangle edges')
        need(sum((t['log_gap'] for t in bound['selected_triangles']), F(0)) == bound['log_gap_lower'],
             'triangle packing accounting')
        for s in product((-1, 1), repeat=len(ms)):
            need(sign_product(rings, s) <= bound['all_continuous_phase_upper'], 'all-phase upper vs exact sign product')
            fixed_phase_checks += 1
    # Resource-independent input checks; malformed data must not be silently accepted.
    rejected = 0
    for call in (lambda: simultaneous_pair_maxima([1, 4]),
                 lambda: family_bound([(2, 1), (4, 1)]),
                 lambda: family_bound([(2, -1)]),
                 lambda: family_bound([(2, 1)], center=True),
                 lambda: sign_product([(2, 1)], [0]),
                 lambda: triangle_gap([(2, 1), (4, 2), (8, 3)], [0, 1, 2])):
        try:
            call()
        except ValueError:
            rejected += 1
        else:
            raise ArithmeticError('invalid input was accepted')

    report = {'passed': True, 'python': sys.version.split()[0],
              'optimization_safe': True,
              'binomial_resultants_via_independent_sylvester_determinants': resultant_count,
              'population_multisets_checked': family_count,
              'sign_assignments_exhausted': sign_assignment_count,
              'compatible_multisets': compatible_count,
              'minimum_compatible_total_population_by_ring_count': smallest_totals,
              'twelve_point_product_independently_expanded_in_Q_sqrt3': True,
              'fixed_family_sign_values_checked_against_all_phase_upper': fixed_phase_checks,
              'malformed_inputs_rejected': rejected,
              'continuum_claim': 'proved analytically in CONTINUOUS_PHASES.md, not inferred from sampling',
              'unrestricted_Mn_upper_computed': False,
              'unrestricted_n12_lower_bound_improvement_claimed': False}
    return report, ex


if __name__ == '__main__':
    report, ex = verify()
    base = Path(__file__).resolve().parent
    (base/'certificates.json').write_text(json.dumps(encode(ex), ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    sources = [p for p in base.iterdir() if p.is_file() and p.suffix in ('.py', '.md')]
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    report['source_sha256'] = manifest
    (base/'verification_report.json').write_text(json.dumps(encode(report), ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    (base/'source_manifest.json').write_text(json.dumps({'algorithm': 'sha256', 'files': manifest}, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(encode(report), ensure_ascii=False, indent=2))
