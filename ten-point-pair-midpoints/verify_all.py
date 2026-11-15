"""Exact antipodal pairing and anisotropic local upper bound for the ten-point candidate."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SEPARATION = F(73, 100)
INVERSE_SQUARE_NORM = F(427, 100)
LENGTH_NORM_3 = F(11)
LENGTH_NORM_4 = F(13)
AMPLITUDE_RADIUS = F(1, 150)
MIDPOINT_STATIONARY_RADIUS = F(13, 100)
MIDPOINT_ENERGY_RADIUS = F(17, 100)


def need(condition, detail):
    if not condition:
        raise ArithmeticError(detail)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_manifest():
    raw = (HERE/'source_manifest.json').read_bytes()
    manifest = json.loads(raw)
    need(set(manifest['files']) == {'.gitattributes', 'README.md', 'verify_all.py'},
         'complete own proof file list')
    need(set(manifest['dependencies']) == {
        'candidate-global-screen/symmetry_certificates.json',
        'candidate-global-screen/verify_symmetry.py',
        'ten-point-local-basin/README.md',
        'ten-point-local-basin/source_manifest.json',
        'ten-point-local-basin/verify_all.py',
        'validated-maxima/certificates.json',
        'validated-maxima/dyadic.py',
        'validated-maxima/verify_all.py',
    }, 'complete immediate proof dependencies')
    for directory, section in ((HERE, 'files'), (ROOT, 'dependencies')):
        for name, sha in manifest[section].items():
            actual = hashlib.sha256((directory/name).read_bytes()).hexdigest()
            need(actual == sha, 'fixed input '+name)
    return hashlib.sha256(raw).hexdigest()


def candidate_representatives(symmetry):
    prior = symmetry.load_prior()
    common = json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    catalog = json.loads((ROOT/'candidate-global-screen/symmetry_certificates.json').read_bytes())
    row = next(item for item in common['rows'] if item['n'] == 10)
    entry = next(item for item in catalog['rows'] if item['n'] == 10)
    certified = symmetry.verify_row(entry, prior, common)
    need(certified['passed'] and certified['rotation_order_exact'] == 2,
         'replayed exact half-turn symmetry')
    permutation = entry['rotation_permutation']
    pairs = [(i, j) for i, j in enumerate(permutation) if i < j]
    need(len(pairs) == 5 and sorted(i for pair in pairs for i in pair) == list(range(10)),
         'five distinct antipodal pairs')
    need(pairs[-1] == tuple(certified['radial_layers'][0]['point_indices']),
         'last antipodal pair consists of the two inner points')
    Ball, scale = prior.Ball, prior.SCALE
    error = F(certified['root_refinement_error_upper'])*scale
    units = -((-error.numerator)//error.denominator)
    center = [Ball.dyadic(n, common['center_exponent']) for n in row['center_numerators']]
    coords = [Ball(v.lo-units, v.hi+units) for v in center]
    coords.insert(1, Ball.integer(0))
    mass_scale = symmetry.sqrt_ball(Ball.integer(9)/10, Ball, scale)
    reps = [(mass_scale*coords[2*i], mass_scale*coords[2*i+1]) for i, _ in pairs]
    need(all(F((x.square()+y.square()).lo, scale) > AMPLITUDE_RADIUS**2 for x, y in reps),
         'all amplitudes stay nonzero throughout the accepted region')
    return pairs, reps, Ball, scale


def matrix_bounds(symmetry, basin):
    pairs, reps, Ball, scale = candidate_representatives(symmetry)
    zero = Ball.integer(0)
    real, imag, len3, len4 = ([[zero for _ in pairs] for _ in pairs] for _ in range(4))
    for i, (x, y) in enumerate(reps):
        for j, (s, t) in enumerate(reps[:i]):
            for sign in (-1, 1):
                a, b = x+sign*s, y+sign*t
                a2, b2 = a.square(), b.square()
                squared = a2+b2
                need(F(squared.lo, scale) > SEPARATION**2, 'strict candidate cross-pair separation')
                length = symmetry.sqrt_ball(squared, Ball, scale)
                squared_squared = squared.square()
                weights = ((a2-b2)/squared_squared, -2*a*b/squared_squared,
                           1/(squared*length), 1/squared_squared)
                for matrix, weight in zip((real, imag, len3, len4), weights):
                    matrix[i][i] += weight
                    matrix[j][j] += weight
                    matrix[i][j] -= weight
                    matrix[j][i] -= weight

    # The real symmetric matrix represents the antilinear operator z -> L conjugate(z).
    real_form = ([row_a+row_b for row_a, row_b in zip(real, imag)] +
                 [row_b+[-v for v in row_a] for row_a, row_b in zip(real, imag)])
    norms = [
        ('inverse_square_positive', INVERSE_SQUARE_NORM, real_form),
        ('inverse_square_negative', INVERSE_SQUARE_NORM, [[-x for x in row] for row in real_form]),
        ('length_3', LENGTH_NORM_3, len3),
        ('length_4', LENGTH_NORM_4, len4),
    ]
    reports = {}
    for label, norm, matrix in norms:
        bound = Ball.integer(norm.numerator)/norm.denominator
        shifted = [[(bound if i == j else zero)-value for j, value in enumerate(row)]
                   for i, row in enumerate(matrix)]
        reports[label], _ = basin.rational_spd(shifted, scale, label)
    return pairs, reports


def rational_bounds():
    sqrt2_upper = F(10,7)
    need(sqrt2_upper**2 > 2, 'strict pair difference factor')
    shift = sqrt2_upper*AMPLITUDE_RADIUS
    separation = SEPARATION-shift
    relative_separation = separation/SEPARATION
    need(relative_separation > 0, 'positive displaced cross-pair separation')
    bound0 = INVERSE_SQUARE_NORM+2*shift*LENGTH_NORM_3/relative_separation**3
    bound4 = LENGTH_NORM_4/relative_separation**4

    def correction(radius):
        rho = 2*radius**2/separation**2
        need(rho < 1, 'positive denominator margin on entire midpoint region')
        return 2*radius**2*bound4/(1-rho)

    stationary_norm = bound0+correction(MIDPOINT_STATIONARY_RADIUS)
    energy_loss = 10-2*bound0-correction(MIDPOINT_ENERGY_RADIUS)
    need(stationary_norm < F(499,100), 'strict stationary contraction')
    need(energy_loss > F(1,8), 'strict positive midpoint energy loss factor')
    return stationary_norm, energy_loss


def verify():
    sha = check_manifest()
    symmetry = load('pair_midpoint_symmetry', ROOT/'candidate-global-screen/verify_symmetry.py')
    basin = load('pair_midpoint_basin', ROOT/'ten-point-local-basin/verify_all.py')
    previous, _ = basin.verify()
    need(previous['passed'] and previous['normalization_mass'] == 10
         and previous['candidate_value_cap_proved']
         and previous['equality_only_candidate_orbit_within_this_radius'],
         'replayed prior local L upper bound')
    basin_radius = F(previous['euclidean_orbit_radius'])
    need(2*F(10,9)*AMPLITUDE_RADIUS**2 < basin_radius**2,
         'symmetrized amplitude vector enters prior mass-ten L basin')
    pairs, checks = matrix_bounds(symmetry, basin)
    stationary_norm, energy_loss = rational_bounds()
    return {
        'passed': True,
        'source_manifest_sha256': sha,
        'normalization_mass': 9,
        'candidate_antipodal_pairs_zero_based': pairs,
        'candidate_interpair_separation_strict_lower': str(SEPARATION),
        'candidate_complex_pair_laplacian_norm_strict_upper': str(INVERSE_SQUARE_NORM),
        'candidate_pair_length_laplacian_norm_strict_upper': [str(LENGTH_NORM_3), str(LENGTH_NORM_4)],
        'matrix_checks': checks,
        'joint_amplitude_norm_upper': str(AMPLITUDE_RADIUS),
        'stationary_midpoint_norm_upper': str(MIDPOINT_STATIONARY_RADIUS),
        'stationary_pair_laplacian_norm_rational_upper': str(stationary_norm),
        'stationary_pair_laplacian_norm_strict_upper': '499/100',
        'unconstrained_midpoint_norm_for_energy_comparison': str(MIDPOINT_ENERGY_RADIUS),
        'midpoint_quadratic_loss_coefficient_rational_lower': str(energy_loss),
        'midpoint_quadratic_loss_coefficient_strict_lower': '1/8',
        'prior_mass_ten_L_basin_replayed': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    report = json.dumps(verify(), indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_bytes(report.encode('utf-8'))
    sys.stdout.buffer.write(report.encode('utf-8'))


if __name__ == '__main__':
    main()
