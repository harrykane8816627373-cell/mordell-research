"""Whole fixed-ring-family exclusions from nonlinear deficits and feasible points.

All arithmetic is exact. Pinned nonlinear core/data are accepted again;
the old integer witness is used only through its coordinates.
"""
from copy import deepcopy
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import factorial
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CORE_PATH = 'phase-deficit-sharp/core.py'
WITNESS_PATH = 'symmetry-breaking/integer_witnesses.json'
NONLINEAR_CERTS = (
    'phase-deficit-sharp/certificates/four_squares.json',
    'phase-deficit-sharp/certificates/four_squares_box.json',
    'phase-deficit-sharp/certificates/mixed_twentyfive_cycle.json',
)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def read_claims():
    claims = json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    need(type(claims) is dict and set(claims) == {'format', 'date', 'exponential_order', 'inputs', 'cases'}, 'claim schema')
    need(claims['format'] == 'mordell-family-exclusion-v1' and claims['date'] == '2026-11-08', 'format/date')
    need(claims['exponential_order'] == 8, 'fixed exponential order')
    need(set(claims['inputs']) == {CORE_PATH, WITNESS_PATH, *NONLINEAR_CERTS}, 'pinned input set')
    need(type(claims['cases']) is list and len(claims['cases']) == 3, 'three fixed family claims')
    need([c['nonlinear_certificate'] for c in claims['cases']] == list(NONLINEAR_CERTS), 'fixed family coverage')
    return claims


def verify_pins(pins):
    need(set(pins) == {CORE_PATH, WITNESS_PATH, *NONLINEAR_CERTS}, 'exact external input set')
    for name, expected in pins.items():
        need(type(expected) is str and len(expected) == 64, 'sha256 hexadecimal input')
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, 'changed pinned external input: '+name)


def load_core():
    # Called only after validating the exact source hash.
    spec = spec_from_file_location('mordell_nov8_nonlinear_core', ROOT/CORE_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def feasible_value(n, points):
    need(type(n) is int and 2 <= n <= 100, 'integer benchmark population')
    need(type(points) is list and len(points) == n, 'benchmark point count')
    need(all(type(p) is list and len(p) == 2 and all(type(q) is int for q in p) for p in points), 'integer coordinate pairs required')
    sx, sy = sum(x for x, y in points), sum(y for x, y in points)
    Dmass = n*sum(x*x+y*y for x, y in points)-sx*sx-sy*sy
    product = 1
    mass_from_pairs = 0
    for (x, y), (a, b) in combinations(points, 2):
        d2 = (x-a)**2+(y-b)**2
        need(d2 > 0, 'benchmark collision')
        mass_from_pairs += d2
        product *= d2
    need(Dmass == mass_from_pairs > 0, 'independent benchmark centered mass identity')
    N = n*(n-1)//2
    V = F(n**(2*N)*product, Dmass**N)
    # A second normalization uses explicit integer centered coordinates.
    centered = [(n*x-sx, n*y-sy) for x, y in points]
    Qc = sum(x*x+y*y for x, y in centered)
    need(sum(x for x, y in centered) == sum(y for x, y in centered) == 0, 'exact centered witness')
    need(Qc == n*Dmass, 'independent centered moment')
    direct_product = 1
    for (x, y), (a, b) in combinations(centered, 2):
        direct_product *= (x-a)**2+(y-b)**2
    need(F(n, Qc)**N*direct_product == V, 'independent exact benchmark value')
    return V, {'n': n, 'distance_count': N, 'centered_mass_integer': str(Dmass),
               'sum_integer_coordinates': [sx, sy], 'value_recomputed_exactly': True}


def exp_lower(D, order=8):
    need(type(order) is int and 0 <= order <= 32 and D >= 0, 'nonnegative exponential argument and bounded integer order')
    result = term = F(1)
    for k in range(1, order+1):
        term *= D/F(k)
        result += term
    return result


def decimal_interval(q, places=14):
    need(q >= 0 and type(places) is int and 0 <= places <= 20, 'nonnegative fixed-place interval')
    scale = 10**places
    integer = q.numerator*scale//q.denominator
    def format_integer(v):
        a, b = divmod(v, scale)
        return str(a)+'.'+str(b).zfill(places)
    return {'lower_rational': str(integer)+'/'+str(scale),
            'upper_rational': str(integer+1)+'/'+str(scale),
            'lower_decimal': format_integer(integer), 'upper_decimal': format_integer(integer+1),
            'lower_closed': True, 'upper_open': True}


def check_case(case, nonlinear, value, core):
    fields = {'label', 'n', 'nonlinear_certificate', 'ratio_lower', 'ratio_upper',
              'threshold', 'previous_first_order_excludes', 'previous_quadratic_excludes',
              'previous_quadratic_E8_excludes', 'qualitative_status'}
    need(type(case) is dict and set(case) == fields, 'case schema')
    accepted = core.check_certificate(nonlinear)
    n = sum(nonlinear['populations'])+nonlinear['center']
    need(type(case['n']) is int and case['n'] == n, 'family and benchmark population match')
    D, G = accepted['nonlinear_gap'], accepted['old_gap']
    U = F(nonlinear['pair_box_upper'])
    E = exp_lower(D)
    need(E == sum((D**k/F(factorial(k)) for k in range(9)), F(0)), 'independent E8 expansion')
    upper = U/E
    ratio = upper/value
    lo, hi, threshold = map(F, (case['ratio_lower'], case['ratio_upper'], case['threshold']))
    need(0 <= lo < ratio < hi < threshold <= 1, 'strict whole-family benchmark comparison')
    previous_first_order = U/(1+D)/value
    previous_quadratic = U/(1+G)/value
    previous_quadratic_E8 = U/exp_lower(G)/value
    expected_old = (previous_first_order < 1, previous_quadratic < 1, previous_quadratic_E8 < 1)
    flags = (case['previous_first_order_excludes'], case['previous_quadratic_excludes'], case['previous_quadratic_E8_excludes'])
    need(all(type(v) is bool for v in flags) and flags == expected_old, 'qualitative comparison to previous bounds')
    if n == 16:
        need(threshold == F(1, 20) and case['qualitative_status'] == 'new exclusion from the displayed certified bounds', 'sixteen-point consequence scope')
        need(expected_old == (False, False, False), 'do not call an existing displayed exclusion new')
    else:
        need(n == 25 and threshold == 1 and case['qualitative_status'] == 'already excluded by previous quadratic bound; stronger arithmetic only', 'twentyfive-point status')
        need(expected_old == (True, True, True), 'existing twentyfive-point exclusion retained')
    return {'label': case['label'], 'n': n, 'nonlinear_certificate_accepted_again': True,
            'nonlinear_log_gap': str(D), 'E8_exact': str(E), 'upper_to_feasible_benchmark': decimal_interval(ratio),
            'previous_first_order_ratio': decimal_interval(previous_first_order),
            'previous_quadratic_ratio': decimal_interval(previous_quadratic),
            'previous_quadratic_E8_ratio': decimal_interval(previous_quadratic_E8),
            'strict_exclusion_threshold': str(threshold), 'entire_family_excluded': True,
            'qualitative_status': case['qualitative_status'],
            'meaning': 'every value in this specified family is below an independently feasible n-point value; not an upper bound for M_n'}


def source_paths():
    return sorted(p for p in HERE.iterdir() if p.is_file() and
                  (p.suffix in ('.py', '.md') or p.name in ('.gitattributes', 'certificates.json')))


def run():
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    need(set(manifest) == {'algorithm', 'files'} and manifest['algorithm'] == 'sha256', 'source manifest schema')
    own_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths()}
    need(manifest['files'] == own_hashes, 'local proof or verifier input changed')
    claims = read_claims()
    verify_pins(claims['inputs'])
    core = load_core()
    baseline = json.loads((ROOT/WITNESS_PATH).read_text(encoding='utf-8'))
    rows = baseline['rows']
    benchmarks, metadata, integer_points = {}, {}, {}
    invariant_checks = 0
    for n in (16, 25):
        matching = [row for row in rows if row['n'] == n]
        need(len(matching) == 1, 'unique pinned benchmark row')
        points = matching[0]['integer_points']
        V, meta = feasible_value(n, points)
        benchmarks[n], metadata[n], integer_points[n] = V, meta, points
        # Translation, scaling and ordering must not change the normalized value.
        for altered in ([[x+1234567, y-9876543] for x, y in points],
                        [[2*x, 2*y] for x, y in points], list(reversed(points))):
            need(feasible_value(n, altered)[0] == V, 'benchmark similarity/order invariance')
            invariant_checks += 1
    results, input_certificates = [], {}
    for case in claims['cases']:
        nonlinear = json.loads((ROOT/case['nonlinear_certificate']).read_text(encoding='utf-8'))
        input_certificates[case['nonlinear_certificate']] = nonlinear
        results.append(check_case(case, nonlinear, benchmarks[case['n']], core))

    rejected = []
    def rejection(label, fn):
        try:
            fn()
        except (ValueError, ZeroDivisionError):
            rejected.append(label)
        else:
            raise ArithmeticError('bad comparison input accepted: '+label)
    altered = deepcopy(claims['inputs']); altered[WITNESS_PATH] = '0'*64
    rejection('wrong pinned benchmark hash', lambda: verify_pins(altered))
    altered_core_pin = deepcopy(claims['inputs']); altered_core_pin[CORE_PATH] = '0'*64
    rejection('wrong pinned nonlinear core hash', lambda: verify_pins(altered_core_pin))
    bad_points = deepcopy(integer_points[16]); bad_points[1] = bad_points[0]
    rejection('colliding feasible witness', lambda: feasible_value(16, bad_points))
    floating_points = deepcopy(integer_points[16]); floating_points[0][0] = float(floating_points[0][0])
    rejection('floating witness coordinate', lambda: feasible_value(16, floating_points))
    base_case = claims['cases'][0]
    original = input_certificates[base_case['nonlinear_certificate']]
    mutations = [
        ('incorrect population', 'n', 17),
        ('false ratio bracket', 'ratio_upper', '1/1000000'),
        ('false prior-exclusion status', 'previous_first_order_excludes', True),
        ('false qualitative claim', 'qualitative_status', 'unrestricted global sharpness'),
    ]
    for label, key, value in mutations:
        bad_case = deepcopy(base_case); bad_case[key] = value
        rejection(label, lambda bad_case=bad_case: check_case(bad_case, original, benchmarks[16], core))
    bad_nonlinear = deepcopy(original); bad_nonlinear['nonlinear_gap'] = '100'
    rejection('unaccepted nonlinear gap', lambda: check_case(base_case, bad_nonlinear, benchmarks[16], core))
    bad_cover = deepcopy(original); bad_cover['cycles'][0]['supports'][0]['cover'].pop()
    rejection('unaccepted angular continuum cover', lambda: check_case(base_case, bad_cover, benchmarks[16], core))

    report = {'passed': True, 'python': sys.version.split()[0], 'optimization_safe': True,
              'pinned_nonlinear_certificates_accepted_again': len(results),
              'feasible_integer_benchmarks_recomputed': metadata,
              'benchmark_similarity_invariance_checks': invariant_checks,
              'rejected_bad_inputs': rejected, 'cases': results,
              'inputs_sha256': claims['inputs'], 'source_sha256': own_hashes,
              'scope': 'excludes exactly the three stated regular-ring families; unrestricted M_n sharpness remains unproved',
              'unrestricted_Mn_upper_claimed': False}
    (HERE/'verification_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'passed': True, 'nonlinear_certificates_reaccepted': len(results),
                      'exact_benchmarks': 2, 'bad_inputs_rejected': len(rejected),
                      'ratio_intervals': {row['label']: row['upper_to_feasible_benchmark'] for row in results}}, indent=2))


if __name__ == '__main__':
    run()
