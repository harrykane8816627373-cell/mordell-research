"""Replay the first-harmonic midpoint bound on the refined radial bands."""
from fractions import Fraction as F
from math import comb
from hashlib import sha256
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
COEFFICIENTS = tuple(map(F, ('-.0080963', '-.0817488', '.6778129',
                            '-.9620101', '-.1909644', '.1782666')))
ERROR_BOUNDS = tuple(map(F, ('.00534', '.02000')))
SLOPE = F('-.0432535')
SUM_NORM = F(239, 1000)


def need(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def rational(obj):
    return F(int(obj['num']), int(obj['den']))


def verify(data, ward, bands, total_budget):
    spec = importlib.util.spec_from_file_location(
        'radial_polynomial_tools', ROOT/'ten-point-ward-inner-count/verify_radial.py')
    radial = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(radial)
    need(F(0) < total_budget < 1, 'ordered Gram energy budget')
    transform = [[ward.hx(x) for x in row] for row in data['transforms'][0]]
    need(transform == [[F(i == j) for j in range(6)] for i in range(6)],
         'first harmonic has the untransformed power basis')
    Q = [[ward.hx(x) for x in row[2:]] for row in data['Q_blocks'][1][2:]]
    ward.gram.positive_ldl(Q)
    saved = (HERE/'jacobian_midpoint_solve.json').read_bytes()
    native = json.loads(saved)
    need(native['operation_id'] == 'matrix.rational_linear_system.solve', 'rational solve operation')
    result = native['output']
    need(result['outcome'] == 'UNIQUE'
         and [[rational(x) for x in row] for row in result['matrix']['entries']] == Q
         and tuple(map(rational,result['rhs'])) == COEFFICIENTS[2:],
         'same restricted Gram system')
    y = list(map(rational, result['solution']))
    need(len(y) == 4 and
         tuple(sum(a*b for a,b in zip(row,y)) for row in Q) == COEFFICIENTS[2:],
         'exact solution substitution')
    inverse_norm_square = sum(a*b for a,b in zip(COEFFICIENTS[2:],y))
    need(inverse_norm_square > 0, 'positive squared dual norm')

    degree = len(COEFFICIENTS) - 1
    polynomial = [sum(COEFFICIENTS[j] * comb(j,i) * (-1)**(j-i)
                      for j in range(i, degree+1)) for i in range(degree+1)]
    covers = []
    for is_inner, band, error in zip((True, False), bands, ERROR_BOUNDS):
        difference = polynomial[:]
        difference[0] -= int(is_inner)
        affine = [error - F(9,10)*SLOPE, SLOPE]
        need(all(affine[0] + SLOPE * u > 0 for u in band),
             'strictly positive affine error bound')
        square = radial.mul(affine, affine)
        margin = radial.mul([F(0), F(-1)], radial.mul(difference, difference))
        for index, value in enumerate(square):
            margin[index] += value
        covers.append(radial.strictly_above(margin, band, F(0)))

    summed_error = 2*ERROR_BOUNDS[0] + 8*ERROR_BOUNDS[1]
    remainder = SUM_NORM - summed_error
    need(remainder > 0 and remainder**2 > total_budget*inverse_norm_square,
         'error plus Gram norm strictly below the stated sum cap')
    return {
        'passed': True,
        'jacobian_solve_sha256': sha256(saved).hexdigest(),
        'first_harmonic_gram_budget_strict_upper': str(total_budget),
        'indicator_approximation_coefficients': list(map(str,COEFFICIENTS)),
        'affine_error_intercepts_at_nine_tenths': list(map(str,ERROR_BOUNDS)),
        'affine_error_slope': str(SLOPE),
        'summed_indicator_error_upper': str(summed_error),
        'squared_inverse_gram_norm_decimal_display': float(inverse_norm_square),
        'two_inner_point_sum_norm_strict_upper': str(SUM_NORM),
        'inner_midpoint_norm_strict_upper': str(SUM_NORM/2),
        'inner_midpoint_squared_strict_upper': str((SUM_NORM/2)**2),
        'outer_centroid_norm_strict_upper': str(SUM_NORM/8),
        'indicator_error_covers': covers,
    }

