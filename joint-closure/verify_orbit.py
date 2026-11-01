"""Integer interval existence/uniqueness and full-Hessian exclusion of one orbit box."""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import json
import argparse
import interval_engine as E


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def interval_det(matrix):
    n = len(matrix)
    out = E.ci(0)
    for perm in permutations(range(n)):
        inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n))
        value = E.ci(-1 if inversions % 2 else 1)
        for i, j in enumerate(perm):
            value = E.cm(value, matrix[i][j])
        out = E.ca(out, value)
    return out


def polynomial_K(q, modes):
    powers, action, _ = E.polynomial_columns(q)
    matrix = []
    for j in modes:
        row = []
        for k in modes:
            value = E.ci(0)
            for i in range(8):
                value = E.ca(value, E.cs(E.cscale(E.cm(E.cc(powers[i][j]), powers[i][k]), 64),
                                           E.cm(E.cc(action[i][j]), action[i][k])))
            row.append(value)
        matrix.append(row)
    return matrix


def fmt(interval):
    return [str(F(v, E.S)) for v in interval]


def check(path):
    data = json.loads(path.read_text())
    E.setup(data['precision_bits'])
    midpoint, radius, R = data['center'], data['root_radius'], data['preconditioner']
    need(len(midpoint) == 13 and len(R) == 13 and all(len(row) == 13 for row in R), '13 actual-coordinate variables')
    box = [(x-radius, x+radius) for x in midpoint]
    need(box[0][0] > 0, 'rotation gauge anchor is positive')
    q = E.positions(box)
    inverse, distances = E.invert_edges(q)
    need(all(d[0] > 0 for d in distances), 'the entire root box is collision-free')
    J = E.polynomial_jacobian(q)
    g = E.polynomial_equations(E.positions([E.point(x) for x in midpoint]))
    qbound, selfmap = 0, 0
    for i in range(13):
        correction = E.ZERO
        for j in range(13):
            correction = E.add(correction, E.mul(E.point(R[i][j]), g[j]))
        row_bound = 0
        for k in range(13):
            value = E.rat(int(i == k))
            for j in range(13):
                value = E.sub(value, E.mul(E.point(R[i][j]), J[j][k]))
            row_bound += max(abs(value[0]), abs(value[1]))
        qbound = max(qbound, row_bound)
        # This bounds x-RG(x)-center on the entire closed coordinate cube.
        error = max(abs(correction[0]), abs(correction[1])) + -(-row_bound*radius//E.S)
        selfmap = max(selfmap, error)
    need(qbound < E.S and selfmap < radius, 'strict contraction and strict self-map')
    # ||I-RJ||<1 implies R nonsingular, so the fixed point is a zero of all 13 equations.
    _, sums = E.power_columns(q)
    a2, b2 = E.cnorm(sums[2]), E.cnorm(sums[3])
    need(a2[0] > 0 and b2[0] > 0, 'both low moments provably nonzero at the exact root')
    need(175*a2[1]+54*b2[1] < 175*E.S, 'exact root is inside the inherited ellipse')
    eps = E.ZERO
    for z in q:
        t = E.cnorm(z)
        eps = E.add(eps, E.mul(t, E.square(E.sub(t, E.ONE))))
    need(eps[0]*1200 > E.S, 'exact root lies beyond the certified radial exclusion')
    target_product = E.ONE
    for distance in distances:
        target_product = E.mul(target_product, distance)
    need(target_product[1] < 7**7*E.S, 'this physical stationary orbit has objective below candidate')
    K = polynomial_K(q, [0, 2, 3, 4])
    minors = []
    for n in range(1, 5):
        determinant = interval_det([row[:n] for row in K[:n]])
        need(determinant[0][0] > 0 and determinant[1][0] <= 0 <= determinant[1][1],
             'first four-mode Hermitian block is strictly positive, minor '+str(n))
        minors.append(fmt(determinant[0]))
    hbox = [(x-data['Hessian_box_radius'], x+data['Hessian_box_radius']) for x in midpoint]
    hq = E.positions(hbox)
    hi, _ = E.invert_edges(hq)
    Hbound = E.rayleigh(hq, hi, data['integer_displacement'])
    need(Hbound[1] < 0, 'entire larger actual-coordinate box violates full real Hessian PSD')
    pbox = [(x-data['polynomial_box_radius'], x+data['polynomial_box_radius']) for x in midpoint]
    Pbound = E.polynomial_rayleigh(E.positions(pbox), data['integer_polynomial_coefficients'])
    need(Pbound[1] < 0, 'entire polynomial box violates the denominator-free full stability matrix')
    return {'passed': True, 'precision_bits': E.BITS,
        'unique_exact_physical_stationary_root_in_box': True,
        'root_box_radius': str(F(radius, E.S)), 'contraction_factor_upper': str(F(qbound, E.S)),
        'self_map_radius_upper': str(F(selfmap, E.S)),
        's2_modulus_squared': fmt(a2), 's3_modulus_squared': fmt(b2),
        'ellipse_upper': str(F(175*a2[1]+54*b2[1], 175*E.S)),
        'E_interval': fmt(eps),
        'Delta_over_target_interval': [str(F(v, E.S*7**7)) for v in target_product],
        'first_four_modes': [0, 2, 3, 4], 'first_four_mode_block_positive_definite': True,
        'first_four_mode_leading_minors': minors,
        'larger_Hessian_box_radius': str(F(data['Hessian_box_radius'], E.S)),
        'full_Hessian_integer_direction_bound': fmt(Hbound),
        'polynomial_box_radius': str(F(data['polynomial_box_radius'], E.S)),
        'full_polynomial_integer_direction_bound': fmt(Pbound),
        'all_global_stationary_orbits_exhausted': False,
        'full_eight_point_exclusion_proved': False}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('certificate', type=Path)
    args = p.parse_args()
    print(json.dumps(check(args.certificate), indent=2))
