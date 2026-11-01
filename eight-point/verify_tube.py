"""Exact constants for a p7-free, rotation-invariant eight-point exclusion tube.

The analytic bridges (Rouche and integration of the real Hessian) are proved
in MOMENT_TUBE.md. Rational checks do not replace those bridges.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'power-moment'))
from complex_interval import C
from witnesses import normalized
from verify_isolation import run as verify_reference_spectrum

INDICES = (2, 3, 4, 5, 6, 8)


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def coefficient_errors(epsilon):
    """c7 is FREE. Neither c1*s7 nor c7*s1 occurs because s1=c1=0."""
    epsilon = F(epsilon)
    require(epsilon > 0, 'positive moment cap')
    bounds = [F(0)] * 9
    for k in INDICES:
        bounds[k] = (epsilon + epsilon * sum(bounds[2:k-1])) / k
    return bounds


def rouche_bounds(bounds, scale_lower, radius):
    scaled = [bounds[k] / scale_lower**k for k in range(9)]
    inner_boundary = radius * (1 - radius**7)
    outer_boundary = (1-radius) * (
        7*radius - sum(F(comb(7, j))*radius**j for j in range(2, 8)))
    inner_error = sum(scaled[k]*radius**(8-k) for k in INDICES)
    outer_error = sum(scaled[k]*(1+radius)**(8-k) for k in INDICES)
    require(inner_error < inner_boundary, 'central Rouche disk')
    require(outer_error < outer_boundary, 'seven outer Rouche disks')
    return {
        'inner_boundary_lower': inner_boundary, 'inner_error_upper': inner_error,
        'outer_boundary_lower': outer_boundary, 'outer_error_upper': outer_error,
    }


def certify(epsilon=F(1, 40), radius=F(7, 800), rho=F(1, 20)):
    epsilon, radius, rho = F(epsilon), F(radius), F(rho)
    require(0 < radius < F(1, 10) and 0 < rho < F(1, 10), 'local radii')
    bounds = coefficient_errors(epsilon)
    max_root_lower, scale_lower = F(187, 200), F(23, 25)
    require(max_root_lower**2 < F(7, 8), 'mass gives a sufficiently large root')
    relative_error = sum(bounds[k]/max_root_lower**k for k in INDICES)
    require(relative_error < 1, 'free c7 cannot vanish')
    coefficient_modulus_lower = max_root_lower**7 * (1-relative_error)
    require(coefficient_modulus_lower > scale_lower**7, 't=abs(c7)^(1/7)>23/25')
    rouche = rouche_bounds(bounds, scale_lower, radius)
    # Minimum reference separation in q=w/2 is >6/7 (legacy exact check).
    require(2*radius < F(6, 7), 'eight disjoint root disks')
    slope = F(15, 14)
    require(slope**2 > F(8, 7) and slope*radius < 1, 'reverse triangle inequality')
    scale_upper = F(100, 99)
    require(1/(1-slope*radius) < scale_upper, 'normalization fixes the free radial scale')
    distance_squared_upper = 32*scale_upper*radius**2
    require(distance_squared_upper < rho**2, 'normalized roots enter the convexity ball')
    # For h in C^8, sum_(i<j)|h_i-h_j|^2 <= 8*||h||^2.
    sqrt2_upper = F(99, 70)
    require(sqrt2_upper**2 > 2, 'sqrt(2) rational upper enclosure')
    separation = F(12, 7)
    require(separation > sqrt2_upper*rho, 'no collisions in the ball')
    perturbation = 16*sqrt2_upper*rho/(separation-sqrt2_upper*rho)**3
    require(F(2, 7)-perturbation > F(1, 35), 'uniform Hessian gap greater than 1/35')
    return {
        'epsilon_in_s_k_equals_p_k_over_2_to_k': epsilon,
        'restricted_indices': INDICES, 'p7_restricted': False,
        'force_one_moment_caps': {k: 2**k*epsilon for k in INDICES},
        'Newton_coefficient_error_bounds_c7_free': {k: bounds[k] for k in INDICES},
        'max_root_lower': max_root_lower, 'free_c7_modulus_lower': coefficient_modulus_lower,
        'scale_lower': scale_lower, 'scale_upper': scale_upper,
        'root_disk_radius_in_q_over_t': radius, **rouche,
        'Euclidean_distance_squared_upper_in_w': distance_squared_upper,
        'Euclidean_ball_radius_in_w': rho,
        'Hessian_operator_perturbation_upper': perturbation,
        'uniform_slice_Hessian_gap_lower': F(1, 35),
        'log_product_ratio_upper': '-distance_to_candidate_orbit_squared/35',
    }


def abs2_upper(value):
    return max(map(abs, value.re))**2 + max(map(abs, value.im))**2


def abs2_lower(value):
    def distance(interval):
        return F(0) if interval[0] <= 0 <= interval[1] else min(map(abs, interval))
    return distance(value.re)**2 + distance(value.im)**2


def run():
    legacy = verify_reference_spectrum()
    require(legacy['passed'] is True and legacy['singular_squares'] ==
            [0, 1, F(1, 2), F(1, 2), F(1, 4), F(1, 4), F(1, 16), F(1, 16)],
            'independently regenerate actual reference forces and complete spectrum')
    result = certify()
    # Deliberately wrong constants must fail; no assert can disappear under -O.
    rejected = []
    for name, kwargs in [('oversized_moment_cap', {'epsilon': F(1)}),
                         ('too_small_root_disks', {'radius': F(1, 10000)}),
                         ('oversized_convexity_ball', {'rho': F(9, 100)})]:
        try:
            certify(**kwargs)
        except ArithmeticError:
            rejected.append(name)
        else:
            raise ArithmeticError('unsafe parameter mutation accepted: '+name)
    source = json.loads((HERE/'witness.json').read_text(encoding='utf-8'))
    raw = [C(*v) for v in source['raw_points']]
    centered, ratio, powers, product = normalized(raw, 28)
    for k in INDICES:
        require(abs2_upper(powers[k]) <= (F(2**k, 40))**2, 'actual witness in the new tube')
    require(abs2_lower(powers[6]) > F(1, 10)**2, 'witness outside every rotation of the old local box')
    require(abs2_lower(powers[2])*abs2_lower(powers[3]) > 0, 'nonzero p2 and p3')
    target = 4**28*7**7
    require(0 < product < target, 'independent direct pair-product comparison')
    residual = sum((1/(centered[0]-centered[j]) for j in range(1, 8)), C(0))
    residual -= ratio*centered[0].conjugate()
    require(residual.norm() > 0, 'normalized witness is explicitly nonstationary')
    return {
        'passed': True, 'date': '2026-11-01', 'n': 8, 'force_constant': 1,
        'inherited_exact_reference_spectrum_recomputed': True,
        'tube_certificate': result,
        'old_uniform_radius_1_over_10_box_strictly_contained': True,
        'nonempty_nonzero_moment_witness_outside_old_rotation_orbits': True,
        'witness_exact_product_over_candidate': product/target,
        'invalid_parameter_cases_rejected': rejected,
        'full_eight_point_exclusion_proved': False,
        'remaining_gap': 'No complete cover of the complement of this tube and the legacy excluded regions.',
    }


if __name__ == '__main__':
    report = encode(run())
    (HERE/'tube_verification.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('PASS eight-point/verify_tube.py')
