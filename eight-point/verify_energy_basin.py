"""Exact larger objective basin and Parseval-assisted p7-free moment tube.

Analytic bridges are in ENERGY_BASIN.md. No global coverage or complete
classification of stationary shapes is inferred from these checks.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'power-moment'))
from complex_interval import C, I, as_i, outward, root_bounds
from witnesses import normalized
from verify_isolation import (Z, ONE, add, scale, mul, omega, conjugate,
    inv_one_minus, sum_fields, mm, poly_mul, run as verify_legacy_spectrum)
from verify_tube import INDICES, coefficient_errors, abs2_upper, abs2_lower, encode, require

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def positive_laplacian():
    """Reconstruct each actual 1/|w_i-w_j|^2 in Q(zeta_7)."""
    inverse = [[Z]*8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            if i == j:
                continue
            if i == 0:
                inverse[i][j] = scale(omega(-(j-1)), F(-1, 2))
            elif j == 0:
                inverse[i][j] = scale(omega(-(i-1)), F(1, 2))
            else:
                inverse[i][j] = scale(mul(omega(-(i-1)), inv_one_minus(j-i)), F(1, 2))
    matrix = [[Z]*8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            if i != j:
                weight = mul(inverse[i][j], conjugate(inverse[i][j]))
                require(weight == conjugate(weight), 'real squared-distance weight')
                matrix[i][j] = scale(weight, -1)
        matrix[i][i] = scale(sum_fields(matrix[i]), -1)
    basis = [[ONE]+[Z]*7] + [[Z]+[omega(k*j) for j in range(7)] for k in range(7)]
    action = [[F(0)]*8 for _ in range(8)]
    action[0][0], action[1][0] = F(7, 4), F(-1, 4)
    action[0][1], action[1][1] = F(-7, 4), F(1, 4)
    for k in range(1, 7):
        action[k+1][k+1] = F(1, 4)+F(k*(7-k), 8)
    metric = [1]+[7]*7
    for i in range(8):
        for j in range(8):
            require(metric[i]*action[i][j] == metric[j]*action[j][i], 'Hermitian Fourier action')
    for k, vector in enumerate(basis):
        actual = [sum_fields([mul(matrix[i][j], vector[j]) for j in range(8)]) for i in range(8)]
        expected = [sum_fields([scale(basis[j][i], action[j][k]) for j in range(8)]) for i in range(8)]
        require(actual == expected, 'actual distance Laplacian Fourier action')
    power = [[F(int(i == j)) for j in range(8)] for i in range(8)]
    traces = [F(8)]
    for _ in range(8):
        power = mm(power, action)
        traces.append(sum((power[i][i] for i in range(8)), F(0)))
    polynomial = [F(1)]
    for k in range(1, 9):
        polynomial.append(-sum((polynomial[j]*traces[k-j] for j in range(k)), F(0))/k)
    spectrum = [0, 2, 1, F(3, 2), F(7, 4), F(7, 4), F(3, 2), 1]
    expected = [F(1)]
    for value in spectrum:
        expected = poly_mul(expected, [1, -F(value)])
    require(polynomial == expected and max(spectrum) == 2, 'complete positive Laplacian spectrum')
    for k in INDICES:
        for ell in INDICES:
            actual = sum_fields([omega(((8-k)-(8-ell))*j) for j in range(7)])
            require(actual == scale(ONE, 7 if k == ell else 0), 'exact Parseval columns')
    return {'spectrum': spectrum, 'characteristic_coefficients': polynomial,
            'operator_norm': 2, 'Parseval_column_pair_checks': len(INDICES)**2}


def certify(epsilon=F(1, 8), root_radius=F(3, 50), basin_radius=F(1, 5)):
    epsilon, r, rho = F(epsilon), F(root_radius), F(basin_radius)
    require(0 < r < F(1, 5) and 0 < rho < F(1, 2), 'local radii')
    bounds = coefficient_errors(epsilon)
    initial_root_lower, center_radius = F(187, 200), F(1, 20)
    require(initial_root_lower**2 < F(7, 8), 'first root lower bound from total mass')
    initial_relative_error = sum(bounds[k]/initial_root_lower**k for k in INDICES)
    require(initial_relative_error < 1, 'initial free coefficient nonzero')
    initial_a_lower = initial_root_lower**7*(1-initial_relative_error)
    center_error = center_radius**8+sum(bounds[k]*center_radius**(8-k) for k in INDICES)
    require(center_error < initial_a_lower*center_radius, 'one small actual q-root')
    improved_root_lower = F(999, 1000)
    require(improved_root_lower**2 < 1-center_radius**2/7, 'remaining seven roots force a larger maximum')
    improved_relative_error = sum(bounds[k]/improved_root_lower**k for k in INDICES)
    require(improved_relative_error < 1, 'positive improved free coefficient bound')
    a_lower = improved_root_lower**7*(1-improved_relative_error)
    t_lower, t_upper = F(24, 25), F(11, 10)
    require(a_lower > t_lower**7, 'improved radial scale t>24/25')
    d = [bounds[k]/t_lower**k for k in range(9)]
    inner_error = sum(d[k]*r**(8-k) for k in INDICES)
    inner_boundary = r*(1-r**7)
    outer_error = sum(d[k]*(1+r)**(8-k) for k in INDICES)
    outer_boundary = (1-r)*(7*r-sum(F(comb(7, j))*r**j for j in range(2, 8)))
    require(inner_error < inner_boundary and outer_error < outer_boundary, 'eight Rouche disks')
    require(2*r < F(6, 7), 'disjoint disks')
    slope = F(15, 14)
    require(slope*slope > F(8, 7) and slope*r < 1, 'normalization slope')
    require(1/(1-slope*r) < t_upper, 'radial scale upper bound')
    central_remainder = r**7+sum(d[k]*r**(7-k) for k in INDICES if k < 7)
    outer_remainder = (sum(F(comb(8, j))*r**(j-1) for j in range(2, 9)) +
        sum(d[k]*sum(F(comb(8-k, j))*r**(j-1) for j in range(1, 9-k)) for k in INDICES))/7
    remainder = max(central_remainder, outer_remainder)
    require(remainder < F(213, 500), 'root-equation remainder factor <0.426<1')
    linear_root_norm_square = d[8]**2+sum(d[k]**2 for k in INDICES)/7
    root_norm_square = linear_root_norm_square/(1-remainder)**2
    distance_square = 4*t_upper*root_norm_square
    require(distance_square < F(3, 100) < rho*rho, 'Parseval root distance squared <3/100<basin radius squared')
    sqrt2_upper, separation_lower = F(99, 70), F(12, 7)
    require(sqrt2_upper**2 > 2, 'sqrt2 upper bound')
    u = sqrt2_upper*rho/separation_lower
    require(u < 1, 'every logarithmic edge series converges')
    tail = u/3+u*u/4+u**3/(5*(1-u))
    gap = F(2, 7)-4*tail
    require(gap > F(1, 30), 'direct log-product gap >1/30')
    return {'epsilon_in_s_k': epsilon, 'restricted_indices': INDICES, 'p7_restricted': False,
        'force_one_moment_caps': {k: 2**k*epsilon for k in INDICES},
        'initial_c7_modulus_lower': initial_a_lower, 'central_q_disk_radius': center_radius,
        'central_q_disk_error_upper': center_error,
        'improved_c7_modulus_lower': a_lower, 'scale_lower': t_lower, 'scale_upper': t_upper,
        'root_matching_radius': r, 'Rouche_inner_error_upper': inner_error,
        'Rouche_inner_boundary_lower': inner_boundary, 'Rouche_outer_error_upper': outer_error,
        'Rouche_outer_boundary_lower': outer_boundary,
        'central_root_remainder_factor': central_remainder, 'outer_root_remainder_factor': outer_remainder,
        'Parseval_linear_root_norm_squared_upper': linear_root_norm_square,
        'matched_root_norm_squared_upper': root_norm_square,
        'orbit_distance_squared_upper': distance_square, 'objective_basin_radius': rho,
        'relative_edge_change_upper': u, 'logarithmic_tail_factor_upper': tail,
        'exact_quadratic_gap_minus_tail': gap, 'claimed_log_product_gap': F(1, 30),
        'log_product_ratio_upper': '-distance_to_candidate_orbit_squared/30'}


def reference_intervals():
    """Machin + real Taylor remainder; only used to check a geometric witness."""
    def atan_bounds(x, n):
        a = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)), F(0))
        b = a+(-1)**n*x**(2*n+1)/F(2*n+1)
        return min(a, b), max(a, b)
    a, b = atan_bounds(F(1, 5), 12), atan_bounds(F(1, 239), 6)
    lo, hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    mesh = 10**18
    lo = F((lo*mesh).__floor__(), mesh)
    hi = F((hi*mesh).__ceil__(), mesh)
    require(F(3) < lo < hi < F(22, 7), 'rational pi enclosure')
    points = [as_i(0)]
    for j in range(7):
        xlo, xhi = 2*j*lo/7, 2*j*hi/7
        x, delta = (xlo+xhi)/2, (xhi-xlo)/2
        sin_value = sum(((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(28)), F(0))
        cos_value = sum(((-1)**k*x**(2*k)/factorial(2*k) for k in range(28)), F(0))
        error = x**56/factorial(56)+delta
        points.append(I(outward(2*(cos_value-error), 2*(cos_value+error)),
                        outward(2*(sin_value-error), 2*(sin_value+error))))
    return points


def run():
    require(verify_legacy_spectrum()['passed'] is True, 'legacy force and complex stability spectrum')
    laplacian, cert = positive_laplacian(), certify()
    rejected = []
    for name, kwargs in [('oversized_moment_cap', {'epsilon': F(1, 2)}),
                         ('too_small_root_disks', {'root_radius': F(1, 100)}),
                         ('too_large_logarithmic_basin', {'basin_radius': F(49, 100)})]:
        try:
            certify(**kwargs)
        except ArithmeticError:
            rejected.append(name)
        else:
            raise ArithmeticError('unsafe mutation accepted: '+name)
    source = json.loads((HERE/'energy_witnesses.json').read_text(encoding='utf-8'))
    target = 4**28*7**7
    raw = [C(*v) for v in source['tube_raw_points']]
    centered, ratio, powers, product = normalized(raw, 28)
    for k in INDICES:
        require(abs2_upper(powers[k]) <= F(2**k, 8)**2, 'actual normalized witness in fivefold tube')
    require(abs2_lower(powers[2]) > F(1, 10)**2, 'outside every rotation of the preceding six-modulus tube')
    require(abs2_lower(powers[2])*abs2_lower(powers[3]) > 0, 'nonzero p2 and p3')
    residual = sum((1/(centered[0]-centered[j]) for j in range(1, 8)), C(0))-ratio*centered[0].conjugate()
    require(residual.norm() > 0 and 0 < product < target, 'independent force residual and exact pair product')
    subset_checks, smallest_slack = 0, None
    for size in range(2, 8):
        for indices in combinations(range(8), size):
            points = [centered[i] for i in indices]
            average = sum(points, C(0))/size
            slack = ratio*sum(((v-average).norm() for v in points), F(0))-F(size*(size-1), 2)
            require(slack > 0, 'tube witness passes inherited strict subset condition')
            smallest_slack = slack if smallest_slack is None else min(smallest_slack, slack)
            subset_checks += 1
    ellipse = F(2, 7)**2*abs2_upper(powers[2])+F(27, 100)*F(2, 7)**3*abs2_upper(powers[3])
    require(subset_checks == 246 and ellipse < F(64, 49), 'inherited geometric and low-moment conditions')
    raw_geometry = [C(*v) for v in source['basin_raw_points']]
    geo_points, geo_ratio, _, geo_product = normalized(raw_geometry, 28)
    root_scale = I(root_bounds(geo_ratio))
    actual = [as_i(v)*root_scale for v in geo_points]
    reference = reference_intervals()
    matching_distance_upper = sum(abs2_upper(a-b) for a, b in zip(actual, reference))
    central_norm_square = geo_ratio*geo_points[0].norm()
    require(F(1, 20)**2 < central_norm_square < 1, 'one point is too far from every old-orbit reference radius')
    require(matching_distance_upper < F(1, 5)**2, 'geometric witness inside the enlarged objective basin')
    require(0 < geo_product < target, 'independent exact geometric-witness product')
    return {'passed': True, 'date': '2026-11-01', 'n': 8, 'force_constant': 1,
        'positive_distance_laplacian': laplacian, 'larger_tube_certificate': cert,
        'all_six_moment_caps_multiplied_by': 5,
        'tube_witness_outside_previous_tube_all_rotations': True,
        'tube_witness_exact_product_over_candidate': product/target,
        'tube_witness_inherited_strict_subsets_passed': subset_checks,
        'tube_witness_minimum_subset_slack': smallest_slack,
        'tube_witness_old_low_moment_ellipse_upper': ellipse,
        'geometric_witness_orbit_distance_squared_lower': central_norm_square,
        'geometric_witness_matching_distance_squared_upper': matching_distance_upper,
        'geometric_witness_outside_old_radius_1_over_20_basin': True,
        'invalid_parameter_cases_rejected': rejected,
        'full_stationary_uniqueness_in_larger_basin_proved': False,
        'full_eight_point_exclusion_proved': False,
        'remaining_gap': 'The complement of the enlarged tube and previously excluded regions has not been exhaustively excluded.'}


if __name__ == '__main__':
    (HERE/'energy_basin_verification.json').write_text(
        json.dumps(encode(run()), ensure_ascii=False, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('PASS eight-point/verify_energy_basin.py')
