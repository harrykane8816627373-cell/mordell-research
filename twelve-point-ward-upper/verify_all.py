"""Full exact replay of the twelve-point Ward-corrected nonsharp upper bound."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SLACK = F(1, 25000)
RATIO = F(379, 378)
LOG_MARGIN = F(1, 1000000)
HARMONIC_TERMS = (
    {(0, 0, 0): 1},
    {(0, 0, 1): 1},
    {(0, 0, 2): 2, (1, 1, 0): -1},
    {(0, 0, 3): 4, (1, 1, 1): -3},
)


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hex_rational(value):
    return F.from_float(float.fromhex(value))


def ward_matrices(data):
    weights = {name: hex_rational(value) for name, value in data['ward_parameters'].items()}
    need(len(weights) == 24, 'nineteen Ward parameters and five mass parameters')
    matrices = [[[F(0) for j in range(5)] for i in range(5)] for k in range(4)]
    for i in (0, 1):
        for j in range(i, 5):
            value = weights.pop(f'A1_{i}{j}')
            matrices[1][i][j] = value
            matrices[1][j][i] = value
    for k, name, constant in ((0, 'mass', F(1, 12)), (2, 'a2', F(-3, 4)), (3, 'a3', F(-2, 3))):
        v = [constant, F(1), F(0), F(0), F(0)]
        a = [weights.pop(f'{name}_{i}') for i in range(5)]
        matrices[k] = [[v[i]*a[j]+a[i]*v[j] for j in range(5)] for i in range(5)]
    need(not weights, 'no unconstrained correction parameters')
    return matrices


def verify_restrictions(matrices, exact_kernel):
    unit_vectors = [[F(i == j) for j in range(5)] for i in range(5)]
    checks = 0
    for k, constant in enumerate((F(1, 12), None, F(-3, 4), F(-2, 3))):
        spanning_vectors = unit_vectors[2:]
        if constant is not None:
            spanning_vectors = [[F(1), -constant, F(0), F(0), F(0)]]+spanning_vectors
        for left in spanning_vectors:
            for right in spanning_vectors:
                need(exact_kernel.bilinear(left, matrices[k], right) == 0,
                     'Ward form vanishes on its whole stated linear subspace')
                checks += 1
    need(checks == 57, 'complete four-block restriction check')
    return checks


def combined_polynomial(data, matrices, old_cover):
    coefficients, shape = old_cover.polynomial(data)
    powers = [[F(comb(j, i)*(-1)**(j-i)) if i <= j else F(0)
               for j in range(5)] for i in range(5)]

    def add(power, value):
        coefficients[power] = coefficients.get(power, F(0))+value

    for k, matrix in enumerate(matrices):
        for i in range(5):
            for j in range(5):
                value = sum(powers[i][a]*matrix[a][b]*powers[j][b]
                            for a in range(5) for b in range(5))
                for (a, b, c), multiplier in HARMONIC_TERMS[k].items():
                    add((a+i, b+j, c), -multiplier*value)
                add((i+j+k, 0, 0), -value/22)
                add((0, i+j+k, 0), -value/22)
    coefficients = {power: value for power, value in coefficients.items() if value}
    need(shape == (20, 20, 12), 'fixed Bernstein tensor shape')
    need(all(all(0 <= exponent < length for exponent, length in zip(power, shape))
             for power in coefficients), 'every exponent fits the complete tensor')
    need(all(coefficients.get((j, i, k), F(0)) == value
             for (i, j, k), value in coefficients.items()), 'exact u/v symmetry of combined h')
    return coefficients, shape


def check_evaluations(data, matrices, coefficients, exact_kernel):
    b = (hex_rational(data['B'])-66)/12

    def pair_kernel(u, v, t):
        phi_u = [(u-1)**j for j in range(5)]
        phi_v = [(v-1)**j for j in range(5)]
        harmonics = (F(1), t, 2*t*t-u*v, 4*t**3-3*u*v*t)
        correction = sum(H*exact_kernel.bilinear(phi_u, A, phi_v)
                         for H, A in zip(harmonics, matrices))
        return exact_kernel.encoded_kernel(u, v, t, data)+correction

    U = F(121, 72)
    points = [(F(0), F(0), F(0)), (F(1), F(1), F(-1)), (U, U, -U),
              (F(2, 7), F(9, 8), F(-3, 13)), (F(1, 9), F(3, 2), F(1, 5)),
              (F(7, 8), F(1, 3), F(-3, 8)), (F(1, 4), U, F(0)),
              (F(2, 3), F(4, 3), F(1, 3))]
    for u, v, t in points:
        diagonal_sum = pair_kernel(u, u, u)+pair_kernel(v, v, v)
        direct = (6*(u+v)+2*b-diagonal_sum/2)/11-pair_kernel(u, v, t)
        expanded = sum(value*u**i*v**j*t**k for (i, j, k), value in coefficients.items())
        need(direct == expanded, 'independent PSD-plus-Ward two-point expression agrees exactly')
    return len(points)


def replay(data, old_cover, coefficients, shape, progress):
    old_cover.selftest()
    # This isolated module instance supplies the frozen integer conversion and
    # cover algorithms with the already expanded combined polynomial.
    old_cover.polynomial = lambda unused_data: (coefficients, shape)
    initial, actual_shape, initial_box, initial_report = old_cover.initial_coefficients(data)
    need(actual_shape == shape, 'initial tensor uses the checked shape')
    packed = old_cover.Packed(shape)
    raw = (HERE/'cover.txt').read_bytes()
    actions = raw.decode('ascii').strip()
    stack = [(initial_box, initial, 0)]
    visited = splits = accepted = pruned = maximum_depth = 0
    for action in actions:
        need(bool(stack), 'tree action has a remaining box')
        box, control, depth = stack.pop()
        visited += 1
        maximum_depth = max(maximum_depth, depth)
        if action == 'P':
            need(old_cover.prune(box), f'exact domain exclusion at node {visited}')
            pruned += 1
        elif action == 'A':
            if not packed.dominates(control, box, SLACK):
                raise ArithmeticError(f'continuous bound failed at node {visited}, box {box}')
            accepted += 1
        else:
            need(action in '012', 'valid split axis')
            axis = int(action)
            left, right = packed.split(control, axis)
            middle = sum(box[axis])/2
            left_box, right_box = list(box), list(box)
            left_box[axis] = (box[axis][0], middle)
            right_box[axis] = (middle, box[axis][1])
            stack.append((tuple(right_box), right, depth+1))
            stack.append((tuple(left_box), left, depth+1))
            splits += 1
        if progress and visited % 50000 == 0:
            print(f'Ward exact cover {visited}/{len(actions)} nodes', file=sys.stderr, flush=True)
    need(not stack and visited == 2*splits+1 and accepted+pruned == splits+1,
         'complete binary partition with zero pending regions')
    return {'passed': True, 'c': str(SLACK), 'visited': visited, 'split': splits,
            'accepted': accepted, 'domain_discarded': pruned, 'pending': 0,
            'maximum_depth': maximum_depth, 'initial_coefficients': initial_report,
            'tree_sha256': hashlib.sha256(raw).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for name, expected in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, 'fixed Ward input '+name)
    for name, expected in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, 'fixed prior theorem/helper '+name)
    source = (HERE/'kernel.json').read_bytes()
    data = json.loads(source)
    need((data['n'], data['mass'], data['feature_order'], data['harmonics_through']) == (12, 11, 5, 11),
         'fixed twelve-point kernel parameterization')
    gram = load('candidate-kernel/verify_psd.py', 'ward_gram')
    matrices = gram.decode(data['Q_blocks'])
    need(len(matrices) == 12, 'all twelve Gram blocks')
    for matrix in matrices:
        need(all(pivot > 0 for pivot in gram.positive_ldl(matrix)), 'positive exact LDL reconstruction')
    correction = ward_matrices(data)
    exact_kernel = load('candidate-kernel/verify_all.py', 'ward_direct_kernel')
    restrictions = verify_restrictions(correction, exact_kernel)
    old_cover = load('twelve-point-kernel-upper/verify_cover.py', 'ward_old_cover')
    coefficients, shape = combined_polynomial(data, correction, old_cover)
    evaluations = check_evaluations(data, correction, coefficients, exact_kernel)
    cover = replay(data, old_cover, coefficients, shape, args.record_report)
    old_bound = load('twelve-point-kernel-upper/verify_all.py', 'ward_old_bound')
    V11, V12 = old_bound.witness()
    need(old_bound.actual_root_lower() > V12, 'true candidate still strictly exceeds feasible witness')
    B = hex_rational(data['B'])
    log_upper = B+66*SLACK
    need(old_bound.log_lower(RATIO*V11)-log_upper > LOG_MARGIN,
         'strict logarithm comparison to the feasible witness')
    report = {
        'passed': True, 'date': '2026-11-09', 'n': 12,
        'arithmetic': 'Python standard-library exact rational and unbounded integer arithmetic',
        'full_continuous_cover_replayed': True, 'cover': cover,
        'exact_positive_definite_Gram_blocks': len(matrices), 'Ward_parameters': 19,
        'mass_parameters': 5, 'exact_zero_bilinear_restrictions': restrictions,
        'combined_kernel_sum_nonnegative_at_complete_stationary_configurations': True,
        'combined_polynomial_exact_swap_symmetry': True,
        'maximum_exponents': [max(power[i] for power in coefficients) for i in range(3)],
        'independent_exact_polynomial_evaluations': evaluations,
        'kernel_sha256': hashlib.sha256(source).hexdigest(),
        'mass_eleven_logarithmic_global_upper': str(log_upper),
        'global_maximum_over_old_feasible_value_strict_upper': str(RATIO),
        'true_candidate_over_global_maximum_strict_lower': str(1/RATIO),
        'log_comparison_margin_strict_lower': str(LOG_MARGIN),
        'actual_candidate_root_replayed_and_exceeds_witness': True,
        'nonsharp_global_upper_proved': True, 'sharp_global_optimality_proved': False,
        'all_maximizers_classified': False, 'Jacobian_or_SDP_required_for_acceptance': False,
        'input_sha256': manifest['files'], 'dependency_sha256': manifest['dependencies'],
    }
    if args.record_report:
        (HERE/'verification_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
