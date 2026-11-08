"""Standard-library exact acceptance of Nov5 integer witnesses; safe under -O.

Discovery logs, gradients and floating stopping flags never decide acceptance.
The only lower-bound inputs are the integer coordinates and the pinned baseline.
"""
from pathlib import Path
from fractions import Fraction
from math import prod
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASELINE_SHA256 = '73e5f42ef78db59dd0ff3e05227caf75eb3834a02486c509767c8534c3c954d5'
RELATIVE_THRESHOLD = Fraction(1, 10000)


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def raw_value(points):
    if not isinstance(points, list) or len(points) < 2:
        raise ValueError('at least two integer coordinate pairs required')
    if any(not isinstance(p, list) or len(p) != 2 or any(type(x) is not int for x in p) for p in points):
        raise ValueError('integer coordinate pairs required')
    n = len(points)
    need(len(set(map(tuple, points))) == n, 'distinct integer points')
    sx = sum(p[0] for p in points)
    sy = sum(p[1] for p in points)
    mass = n*sum(x*x+y*y for x, y in points)-sx*sx-sy*sy
    distances = [(x-u)**2+(y-v)**2 for i, (x, y) in enumerate(points) for u, v in points[i+1:]]
    need(mass > 0 and min(distances) > 0 and sum(distances) == mass, 'centering and pair-mass identity')
    return prod(distances), mass


def integer_value(points):
    product, mass = raw_value(points)
    n = len(points)
    exponent = n*(n-1)//2
    return Fraction(n**(2*exponent)*product, mass**exponent)


def independently_centered_value(points):
    n = len(points)
    sx, sy = sum(p[0] for p in points), sum(p[1] for p in points)
    centered = [(n*x-sx, n*y-sy) for x, y in points]
    moment = sum(x*x+y*y for x, y in centered)
    result = Fraction(1)
    for i, (x, y) in enumerate(centered):
        for u, v in centered[i+1:]:
            result *= Fraction(n*((x-u)**2+(y-v)**2), moment)
    return result


def decimal_interval(value, digits=12):
    scale = 10**digits
    lower = value.numerator*scale//value.denominator
    def fmt(integer):
        return str(integer//scale)+'.'+str(integer % scale).zfill(digits)
    return [fmt(lower), fmt(lower+1)]


def digest_integer(value):
    return hashlib.sha256(value.to_bytes((value.bit_length()+7)//8, 'big')).hexdigest()


def verify():
    baseline_raw = (ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes()
    need(hashlib.sha256(baseline_raw).hexdigest() == BASELINE_SHA256, 'pinned complete prior witness bytes')
    baseline = {r['n']: r['integer_points'] for r in json.loads(baseline_raw)['rows']}
    witnesses = json.loads((HERE/'integer_witnesses.json').read_text(encoding='utf-8'))
    need(witnesses['baseline_sha256'] == BASELINE_SHA256, 'candidate provenance agrees with pinned baseline')
    counts = [r['n'] for r in witnesses['rows']]
    need(counts == sorted(set(counts)) and counts and all(type(n) is int and 10 <= n <= 100 for n in counts), 'ordered unique point counts')
    examples = [[[0, 0], [1, 0]], [[2, 3], [-1, 7], [4, -2]], [[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1]]]
    for points in examples:
        need(integer_value(points) == independently_centered_value(points), 'independent normalization formula')
        transformed = [[7*x+11, 7*y-23] for x, y in points]
        need(integer_value(points) == integer_value(transformed), 'exact translation and dilation invariance')
    malformed = [[], [[0, 0]], [[0, 0], [0, 0]], [[0, 0], [True, 1]], [[0, 0], [1.0, 2]], [[0, 0], [1, 2, 3]]]
    rejected = 0
    for points in malformed:
        try:
            raw_value(points)
        except (ValueError, ArithmeticError):
            rejected += 1
    need(rejected == len(malformed), 'malformed witnesses rejected')
    rows = []
    for row in witnesses['rows']:
        n = row['n']
        points = row['integer_points']
        need(len(points) == n, 'actual witness point count')
        new_product, new_mass = raw_value(points)
        old_product, old_mass = raw_value(baseline[n])
        exponent = n*(n-1)//2
        numerator = new_product*old_mass**exponent
        denominator = old_product*new_mass**exponent
        need(numerator > denominator, 'strict new attainable lower bound')
        ratio = Fraction(numerator, denominator)
        need(ratio > 1+RELATIVE_THRESHOLD, 'material exact gain threshold')
        need(ratio == integer_value(points)/integer_value(baseline[n]), 'independent normalized-value ratio')
        if n <= 20:
            need(integer_value(points) == independently_centered_value(points), 'independent actual coordinate multiplication')
        rows.append({'n': n, 'ratio_new_over_baseline_interval': decimal_interval(ratio),
                     'relative_gain_percent_interval': decimal_interval(100*(ratio-1)),
                     'relative_gain_basis_points_floor': ((ratio-1)*10000).numerator//((ratio-1)*10000).denominator,
                     'strictly_improves_pinned_baseline': True, 'new_centered_mass_D': str(new_mass),
                     'baseline_centered_mass_D': str(old_mass), 'new_raw_distance_product_sha256': digest_integer(new_product),
                     'ratio_numerator_sha256': digest_integer(ratio.numerator),
                     'ratio_denominator_sha256': digest_integer(ratio.denominator)})
    return {'passed': True, 'date': '2026-11-05', 'optimization_safe': True,
            'baseline_path': 'symmetry-breaking/integer_witnesses.json', 'baseline_sha256': BASELINE_SHA256,
            'point_counts': counts, 'material_relative_gain_threshold': '1/10000',
            'independent_centering_examples': len(examples), 'malformed_inputs_rejected': rejected,
            'scope': 'Exact attainable lower bounds against the pinned saved witnesses only; no global optimality or exact stationary configuration is certified.',
            'rows': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true', help='write a new report after all exact checks pass')
    args = parser.parse_args()
    if not args.record_report:
        manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
        seen = set()
        for row in manifest['files']:
            name = row['path']
            need(type(name) is str and name not in seen and '/' not in name and '\\' not in name
                 and name not in ('.', '..'), 'safe unique local fixed input')
            seen.add(name)
            raw = (HERE/name).read_bytes()
            need(len(raw) == row['bytes'] and hashlib.sha256(raw).hexdigest() == row['sha256'], 'fixed input '+name)
    report = verify()
    path = HERE/'verification_report.json'
    if args.record_report:
        path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    else:
        need(report == json.loads(path.read_text(encoding='utf-8')), 'exact saved report reproduced')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
