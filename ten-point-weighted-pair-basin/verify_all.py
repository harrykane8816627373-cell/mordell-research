"""Check weighted Taylor bounds for the ten-point pair basin using integers."""
import argparse
from fractions import Fraction as F
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module


contact = load('weighted_contact', ROOT / 'ten-point-contact-geometry' / 'verify_all.py')
basin = load('weighted_basin', ROOT / 'ten-point-local-basin' / 'verify_all.py')
Ball, SCALE = contact.Ball, contact.SCALE
ZERO, ONE, exact, need = contact.ZERO, contact.ONE, contact.exact, contact.need
matmul = contact.ref.prior.matmul
transpose = contact.ref.prior.transpose

AMPLITUDE = F(7, 250)
MIDPOINT_ENERGY = F(7, 40)
MIDPOINT_STATIONARY = F(13, 100)
GAP_MINUS, GAP_PLUS = F(17, 100), F(29, 4)
LENGTH_AMPLITUDE_4 = 99
LENGTH_MIDPOINT_RELATIVE_4 = F(13, 2)
CUBIC_PLUS, CUBIC_MIXED, CUBIC_MIDPOINT = 21, 19, 12


def sqrt(x):
    return contact.ref.sym.sqrt_ball(x, Ball, SCALE)


def eye(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def add(a, b, alpha=1, beta=1):
    return [[alpha*x+beta*y for x, y in zip(u, v)] for u, v in zip(a, b)]


def outer(a, b):
    return [[x*y for y in b] for x in a]


def subtract(a, b):
    return add(a, b, beta=-1)


def inverse(a):
    n = len(a)
    rows = [u + v for u, v in zip(a, eye(n))]
    for k in range(n):
        pivot = rows[k][k]
        need(pivot.lo > 0, 'strict midpoint inverse pivot')
        rows[k] = [v/pivot for v in rows[k]]
        for i, row in enumerate(rows):
            if i != k:
                factor = row[k]
                rows[i] = [u-factor*v for u, v in zip(row, rows[k])]
    return [row[n:] for row in rows]


def value(row, v):
    return tuple(sum((coeff*z[k] for coeff, z in zip(v, row)), ZERO) for k in (0, 1))


def normalized_functional(d, p):
    """Coefficients of p.x / d for consecutive real and imaginary coordinates."""
    x, y = d
    squared = x.square()+y.square()
    u, v = x/squared, y/squared
    return [(t*u, -t*v) for t in p]+[(t*v, t*u) for t in p]


def triple_real(a, b, c):
    x, y = a
    u, v = b
    s, t = c
    return (x*u-y*v)*s-(x*v+y*u)*t


def cubic_slices(pieces, first, second):
    """Matrices T_i with T(h,s,s)=sum_i h_i s^T T_i s in the chosen bases."""
    rows = [(coeff, [value(x, b) for b in first], [value(y, b) for b in second])
            for coeff, x, y in pieces]
    return [[[sum((coeff*triple_real(u[i], v[j], v[k]) for coeff, u, v in rows), ZERO)
              for k in range(len(second))] for j in range(len(second))] for i in range(len(first))]


def gram_from_slices(slices, right_inverse=None):
    if right_inverse is None:
        flat = [[x for row in a for x in row] for a in slices]
        return matmul(flat, transpose(flat))
    weighted = [matmul(a, right_inverse) for a in slices]
    return [[sum((x[i][j]*y[j][i] for i in range(len(x)) for j in range(len(x))), ZERO)
             for y in weighted] for x in weighted]


def add_quadratic(target, functional, weight):
    r, i = map(list, zip(*functional))
    return add(target, subtract(outer(r, r), outer(i, i)), beta=weight)


def reflection_bases():
    e = eye(10)
    half = ONE/sqrt(Ball.integer(2))
    first = [e[0]] + [[half*(x+y) for x, y in zip(e[i], e[j])] for i, j in ((1, 2), (3, 4))]
    first += [[half*(x-y) for x, y in zip(e[i], e[j])] for i, j in ((6, 7), (8, 9))]
    second = [[half*(x-y) for x, y in zip(e[i], e[j])] for i, j in ((1, 2), (3, 4))]
    second += [e[5]] + [[half*(x+y) for x, y in zip(e[i], e[j])] for i, j in ((6, 7), (8, 9))]
    return first, second


def algebraic_data():
    contact.manifest_check()
    nodes, xy, _, _, _ = contact.candidate()
    reduced_root, _, _ = contact.check_reduced_box()
    need(all(lo < F(v.lo, SCALE) <= F(v.hi, SCALE) < hi
             for (lo, hi), v in zip(reduced_root, nodes + xy)), 'common reduced and full force root')
    a, b, c, x, y = nodes + xy
    reps = [(sqrt(a), ZERO)]
    for r, t in ((b, x), (c, y)):
        u, v = sqrt((r+t)/2), sqrt((r-t)/2)
        reps.extend([(u, v), (u, -v)])
    return reps


def matrix_bounds(reps):
    coords = eye(5)
    full = eye(10)
    square_amp = add(full, full, alpha=10, beta=0)
    square_mid = square_amp
    length_amp = [[ZERO]*5 for _ in coords]
    length_mid = [[ZERO]*5 for _ in coords]
    cubic_amp, cubic_mid = [], []

    def accumulate_amp(weight, d, p):
        nonlocal square_amp, length_amp
        squared = sum((x.square() for x in d), ZERO)
        need(F(squared.lo, SCALE) > F(1, 9), 'candidate denominators for amplitude remainder')
        fun = normalized_functional(d, p)
        square_amp = add_quadratic(square_amp, fun, weight//2)
        length_amp = add(length_amp, outer(p, p), beta=(weight//2)/squared.square())
        cubic_amp.append((exact(F(weight, 3)), fun, fun))
        return fun, squared

    for i, d in enumerate(reps):
        accumulate_amp(2, d, coords[i])
    for i, j in combinations(range(5), 2):
        m = [a-b for a, b in zip(coords[i], coords[j])]
        for sign in (-1, 1):
            d = tuple(a+sign*b for a, b in zip(reps[i], reps[j]))
            p = [a+sign*b for a, b in zip(coords[i], coords[j])]
            fun, squared = accumulate_amp(4, d, p)
            need(F(squared.lo, SCALE) > F(73, 100)**2, 'strict cross-pair candidate separation')
            mid = normalized_functional(d, m)
            square_mid = add_quadratic(square_mid, mid, 2)
            length_mid = add(length_mid, outer(m, m), beta=ONE/squared.square())
            cubic_mid.append((4, fun, mid))

    plus, minus = reflection_bases()
    q = [-y for _, y in reps]+[x for x, _ in reps]
    q_projector = [[2*x/9 for x in row] for row in outer(q, q)]
    plus_projector = matmul(transpose(plus), plus)
    minus_projector = matmul(transpose(minus), minus)
    mixed = [(3*coeff, a, b) for coeff, a, b in cubic_amp]
    gram_amp_plus = gram_from_slices(cubic_slices(cubic_amp, plus, plus))
    gram_amp_mixed = gram_from_slices(cubic_slices(mixed, plus, minus))
    gram_mid = gram_from_slices(cubic_slices(cubic_mid, full, full), inverse(square_mid))
    plus_quadratic = matmul(matmul(plus, square_amp), transpose(plus))
    real_length_mid = ([row+[ZERO]*5 for row in length_mid] +
                       [[ZERO]*5+row for row in length_mid])

    amplitude_curvature = add(add(square_amp, q_projector), plus_projector, beta=-exact(GAP_PLUS))
    amplitude_curvature = add(amplitude_curvature, minus_projector, beta=-exact(GAP_MINUS))
    bounds = {
        'amplitude_curvature': amplitude_curvature,
        'amplitude_length_4': add(eye(5), length_amp, alpha=LENGTH_AMPLITUDE_4, beta=-1),
        'amplitude_even_cubic': add(eye(5), gram_amp_plus, alpha=CUBIC_PLUS**2, beta=-1),
        'amplitude_mixed_cubic': add(plus_quadratic, gram_amp_mixed, alpha=CUBIC_MIXED**2, beta=-1),
        'midpoint_curvature': add(square_mid, full, beta=-exact(F(7, 5))),
        'midpoint_relative_length_4': add(square_mid, real_length_mid,
                                           alpha=exact(LENGTH_MIDPOINT_RELATIVE_4), beta=-1),
        'midpoint_relative_cubic': add(full, gram_mid, alpha=CUBIC_MIDPOINT**2, beta=-1),
    }
    return {name: basin.rational_spd(matrix, SCALE, name) for name, matrix in bounds.items()}


def remainder_bounds():
    radius = AMPLITUDE
    separation = F(73, 100)
    pair_displacement = F(10, 7)*radius
    x = pair_displacement / separation
    theta = 1-x
    need(1-3*radius > 0 and theta > 0 and F(10, 7)**2 > 2, 'denominator margins')
    length_tail = LENGTH_AMPLITUDE_4*radius**2/(1-3*radius)
    s_diagonal = GAP_MINUS-length_tail
    even_diagonal = 1-(CUBIC_PLUS*radius+length_tail)/GAP_PLUS
    amplitude_det = s_diagonal*even_diagonal-(CUBIC_MIXED*radius/2)**2
    need(s_diagonal > 0 and amplitude_det > F(3, 500), 'positive two-by-two amplitude loss')

    def midpoint_budget(m, factor):
        rho = 2*m**2/(separation-pair_displacement)**2
        need(rho < 1, 'positive midpoint denominator margin')
        inverse_square_tail = 4*radius**2*(3-2*x)*LENGTH_MIDPOINT_RELATIVE_4/theta**2
        midpoint_tail = factor*2*m**2*LENGTH_MIDPOINT_RELATIVE_4/theta**4/(1-rho)
        return CUBIC_MIDPOINT*radius+inverse_square_tail+midpoint_tail

    pos_energy = 1-midpoint_budget(MIDPOINT_ENERGY, 1)
    pos_stationary = 1-midpoint_budget(MIDPOINT_STATIONARY, 2)
    need(pos_energy > F(1, 40) and pos_stationary > F(1, 200), 'positive midpoint loss margins')
    return amplitude_det, pos_energy, pos_stationary


def verify():
    checks = matrix_bounds(algebraic_data())
    det, energy, stationary = remainder_bounds()
    return ({
        'passed': True,
        'normalization_mass': 9,
        'candidate_pair_amplitude_radius': str(AMPLITUDE),
        'midpoint_radius_for_sharp_energy': str(MIDPOINT_ENERGY),
        'midpoint_radius_for_stationary_symmetry': str(MIDPOINT_STATIONARY),
        'global_complement_covered': False,
        'amplitude_loss_matrix_determinant': str(det),
        'midpoint_energy_loss_relative_to_quadratic': str(energy),
        'midpoint_stationary_loss_relative_to_quadratic': str(stationary),
        'matrix_checks': {name: report for name, (report, _) in checks.items()},
    }, {name: payload for name, (_, payload) in checks.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-jacobian', type=Path)
    args = parser.parse_args()
    report, payloads = verify()
    if args.export_jacobian is not None:
        args.export_jacobian.write_text(json.dumps(payloads, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
