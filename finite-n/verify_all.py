"""Replay published exact certificates. Finite checks do not replace proofs."""
from collections import defaultdict
from fractions import Fraction
from math import factorial, comb, prod
from pathlib import Path
import hashlib
import json

from finite_lower_bounds import CoefficientMoments, subset_slack, gaussian_lower, stopping_check
from triangular_lattice import audit as lattice_audit

ROOT = Path(__file__).resolve().parent


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def independent_direct_norm(n, k):
    """Direct pair-factor expansion; independent of partition/sign-mask DP."""
    coefficients = {(0,) * n: 1}
    for i in range(n):
        for j in range(i + 1, n):
            nxt = defaultdict(int)
            for alpha, c in coefficients.items():
                for a in range(k + 1):
                    exponent = list(alpha)
                    exponent[i] += a
                    exponent[j] += k - a
                    nxt[tuple(exponent)] += c * comb(k, a) * (-1) ** a
            coefficients = {alpha: c for alpha, c in nxt.items() if c}
    return sum(c * c * prod(factorial(a) for a in alpha) for alpha, c in coefficients.items())


def main():
    manifest = json.loads((ROOT / 'source_manifest.json').read_text(encoding='utf-8-sig'))
    for row in manifest['files']:
        p = ROOT / row['file']
        check(p.parent == ROOT and p.is_file(), 'invalid source manifest member')
        b = p.read_bytes()
        check(len(b) == row['bytes'] and hashlib.sha256(b).hexdigest() == row['sha256'], 'source digest mismatch: ' + row['file'])
    engines = {n: CoefficientMoments(n) for n in range(2, 8)}
    direct_count = 0
    for n in range(2, 6):
        for k in range(1, 4):
            check(engines[n].coefficient_norm(k) == independent_direct_norm(n, k), 'independent coefficient norm mismatch')
            direct_count += 1
    stored = json.loads((ROOT / 'moment_values.json').read_text(encoding='utf-8-sig'))
    check([row['n'] for row in stored['rows']] == list(range(2, 8)), 'unexpected computed moment range')
    for row in stored['rows']:
        n = row['n']
        engine = engines[n]
        check(row['J'] == [str(engine.coefficient_norm(k)) for k in range(4)], 'stored J mismatch')
        moments = [engine.moment(k) for k in range(4)]
        check(row['moments'] == [str(x) for x in moments], 'stored moment mismatch')
        check(moments[1] == gaussian_lower(n), 'centered Gaussian normalization mismatch')
        rayleigh, shift = engine.rayleigh_lower()
        check(row['rayleigh_lower'] == str(rayleigh) and row['rational_shift'] == str(shift), 'stored Rayleigh mismatch')
        if n >= 3:
            check(rayleigh > moments[3] / moments[2] > moments[2] / moments[1] > moments[1], 'strict moment improvement failed')
        else:
            check(moments == [Fraction(4 ** k) for k in range(4)] and rayleigh == 4, 'two-point exact value failed')
    for k in range(9):
        direct = Fraction(sum(comb(k, j) ** 2 * 9 ** j * factorial(3 * k - 2 * j) * factorial(2 * j) for j in range(k + 1)), 2 ** k)
        check(engines[3].coefficient_norm(k) == direct, 'independent three-point closed formula mismatch')
    lattice_payload, comparison = lattice_audit(30)
    check(json.loads(json.dumps(lattice_payload)) == json.loads((ROOT / 'triangular_lattice_values.json').read_text(encoding='utf-8-sig')), 'lattice coordinate/product certificate mismatch')
    check(comparison.rstrip() == (ROOT / 'triangular_comparison.md').read_text(encoding='utf-8-sig').rstrip(), 'lattice comparison mismatch')
    subset_count = 0
    for n in range(3, 31):
        for m in range(2, n):
            B, epsilon = subset_slack(n, m)
            check(B > 1 and epsilon > 0, 'subset rational constants not positive')
            subset_count += 1
    B, epsilon = subset_slack(8, 7)
    check(B == Fraction(8, 3) and epsilon == Fraction(2187, 9217846), 'eight-point epsilon mismatch')
    check(2 * epsilon / 7 == Fraction(2187, 32262461), 'eight-point inertia increment mismatch')
    check(2 * epsilon / 8 == Fraction(2187, 36871384), 'eight-point radius increment mismatch')
    stopping = []
    for n in range(2, 15):
        k = 0
        while not stopping_check(n, k, Fraction(1, 10)):
            k += 1
        stopping.append({'n': n, 'relative_error': '1/10', 'sufficient_ratio_order': k})
    report = {
        'format': 'mordell-subset-finite-lower-replay-v1',
        'result': 'PASS',
        'source_files_sha256_matched': len(manifest['files']),
        'independent_direct_norm_cases': direct_count,
        'coefficient_moments_replayed_n': list(range(2, 8)),
        'independent_n3_closed_formula_orders': list(range(9)),
        'lattice_certificates_replayed_n': list(range(2, 31)),
        'rational_subset_constant_cases': subset_count,
        'eight_point_subset_constants': 'PASS',
        'integer_stopping_checks': stopping,
        'floating_point_acceptance': False,
        'global_eight_point_sharp_bound_proved': False,
        'finite_cases_do_not_replace_all_n_proofs': True,
    }
    (ROOT / 'verification_report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('PASS: source hashes, independent coefficient expansion, exact moments, lattice certificates, subset constants, integer stopping checks.')


if __name__ == '__main__':
    main()
