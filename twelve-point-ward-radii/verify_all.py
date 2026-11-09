"""Replay the Ward row cuts and the exact radial population certificate."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import lcm, prod
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
U = F(121, 72)
DMIN = F(1, 3)
LEVEL = F(1, 4000)
BANDS = ((F(3, 25), F(6, 25)), (F(47, 50), F(71, 50)))
EXCLUDED = ((F(0), BANDS[0][0]), (BANDS[0][1], BANDS[1][0]),
            (BANDS[1][1], F(33, 20)))


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def initial(coefficients, shape, box, cover):
    denominator = lcm(*(a.denominator for a in coefficients.values()))
    size = prod(shape)
    tensor = [0]*size
    for (i, j, k), value in coefficients.items():
        tensor[(i*shape[1]+j)*shape[2]+k] = value.numerator*(denominator//value.denominator)
    for axis, dimension in enumerate(shape):
        matrix, divisor = cover.transform(dimension-1, *box[axis])
        stride = prod(shape[axis+1:])
        group = dimension*stride
        transformed = [0]*size
        for base in range(0, size, group):
            for offset in range(stride):
                values = tensor[base+offset:base+offset+group:stride]
                for i, row in enumerate(matrix):
                    transformed[base+i*stride+offset] = sum(x*y for x, y in zip(row, values))
        tensor = transformed
        denominator *= divisor
    packed = 0
    for value in reversed(tensor):
        lane = cover.S*value//denominator+cover.BIAS
        need(0 <= lane < 1 << 62, 'initial coefficient lane is safe')
        packed = (packed << cover.W) | lane
    return packed


def outside(box):
    (ul, uh), (vl, vh), (tl, th) = box
    min_abs = max(F(0), tl, -th)
    # The first radius is the marked point. Do not discard boxes with u > v.
    return min_abs**2 > uh*vh or uh+vh-2*tl < DMIN


def replay_rows(h, shape, cover, slack):
    raw = json.loads((HERE/'cuts.json').read_bytes())
    rows = [tuple(map(F, row)) for row in raw]
    trees = (HERE/'cover.txt').read_bytes().splitlines()
    need(len(rows) == len(trees) == 53, 'one tree per affine row cut')
    index = 0
    for start, end in EXCLUDED:
        cursor = start
        while cursor < end:
            lo, hi, _, _ = rows[index]
            need(lo == cursor and lo < hi <= end, 'complete contiguous exclusion interval')
            cursor = hi
            index += 1
    need(index == len(rows), 'every affine cut belongs to the exclusion set')
    cover.selftest()
    packed = cover.Packed(shape)
    totals = dict(nodes=0, splits=0, accepted=0, pruned=0, maximum_depth=0)
    for (lo, hi, a, b), tree in zip(rows, trees):
        box = ((lo, hi), (F(0), U), (-U, U))
        coefficients = dict(h)
        # Subtract a*(v-(11-u)/11) + b*(t+u/11) before checking the gap.
        coefficients[(0, 0, 0)] += a
        coefficients[(1, 0, 0)] -= (a+b)/11
        coefficients[(0, 1, 0)] -= a
        coefficients[(0, 0, 1)] -= b
        stack = [(box, initial(coefficients, shape, box, cover), 0)]
        splits = accepted = pruned = 0
        for action in tree:
            need(bool(stack), 'tree action has a pending box')
            box, control, depth = stack.pop()
            totals['maximum_depth'] = max(totals['maximum_depth'], depth)
            if action == ord('P'):
                need(outside(box), 'exact Gram or minimum distance exclusion')
                pruned += 1
            elif action == ord('A'):
                need(packed.dominates(control, box, slack-LEVEL), 'exact row cut leaf')
                accepted += 1
            else:
                need(action in b'012', 'valid split axis')
                axis = action-ord('0')
                left, right = packed.split(control, axis)
                lo_axis, hi_axis = box[axis]
                mid = (lo_axis+hi_axis)/2
                lbox, rbox = list(box), list(box)
                lbox[axis], rbox[axis] = (lo_axis, mid), (mid, hi_axis)
                stack.append((tuple(rbox), right, depth+1))
                stack.append((tuple(lbox), left, depth+1))
                splits += 1
        need(not stack and len(tree) == 2*splits+1 and accepted+pruned == splits+1,
             'complete row cover with zero pending boxes')
        totals['nodes'] += len(tree)
        totals['splits'] += splits
        totals['accepted'] += accepted
        totals['pruned'] += pruned
    return dict(totals, intervals=len(rows), pending=0)


def rational_value(obj):
    return F(int(obj['num']), int(obj['den']))


def population(data, budget, cover, ward):
    # Q_0 is positive definite by the bound's exact LDL checks.
    Q = [[ward.hex_rational(value) for value in row] for row in data['Q_blocks'][0]]
    means = list(map(ward.hex_rational, data['radial_moment_centers']))
    sub = [row[1:] for row in Q[1:]]
    native = json.loads((HERE/'jacobian_radial_solve.json').read_bytes())
    need(native['operation_id'] == 'matrix.rational_linear_system.solve', 'recorded native solve')
    result = native['output']
    need(result['outcome'] == 'UNIQUE', 'four-dimensional radial system')
    need([[rational_value(x) for x in row] for row in result['matrix']['entries']] == sub,
         'recorded native matrix matches the exact radial block')
    c = [F(-4), F(1), F(0), F(0)]
    y = list(map(rational_value, result['solution']))
    need(list(map(rational_value, result['rhs'])) == c and
         [sum(a*b for a, b in zip(row, y)) for row in sub] == c,
         'independent exact matrix-vector verification')
    dual_norm_squared = sum(a*b for a, b in zip(c, y))
    need(0 < dual_norm_squared < 20, 'dual radial quadratic form bounded by twenty')
    first = 11-12*means[0]
    first_error = abs(first*sum(a*b for a, b in zip(Q[0][1:], y)))
    need(first_error < F(1, 10000), 'first radial coordinate rounding accounted for')
    constant = 4*means[1]-means[2]
    # On each band, certify a lower supporting line for each signed cubic.
    lines = ((1, (F(347, 100), F(509, 100)), F(-26, 5)),
             (-1, (F(-199, 50), F(-263, 50)), F(531, 100)))
    population_bounds = {}
    for sign, intercepts, slope in lines:
        for (lo, hi), intercept in zip(BANDS, intercepts):
            polynomial = [sign*constant-intercept, -slope, -4*sign, F(sign)]
            for i in range(64):
                a = lo+(hi-lo)*i/64
                b = lo+(hi-lo)*(i+1)/64
                bernstein, _ = cover.transform(3, a, b)
                need(all(sum(x*y for x, y in zip(row, polynomial)) >= 0 for row in bernstein),
                     'exact supporting line for signed radial cubic')
        tested = range(3) if sign == 1 else range(4, 13)
        for k in tested:
            lower = k*intercepts[0]+(12-k)*intercepts[1]+11*slope
            need(lower > first_error and (lower-first_error)**2 > 2*budget*dual_norm_squared,
                 'radial moment penalty exceeds the whole Ward budget')
            population_bounds[str(k)] = str(lower)
    return {'possible_inner_population': [3], 'outer_population': 9,
            'signed_radial_moment_lower_bounds_excluded_populations': population_bounds,
            'radial_dual_norm_squared_upper': '20', 'first_coordinate_error_upper': '1/10000',
            'exact_native_linear_solve_verified_by_substitution': True,
            'cubic_bernstein_subintervals': 256}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for name, expected in manifest.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, 'fixed input '+name)
    ward = load('twelve-point-ward-upper/verify_all.py', 'ward_radii_prior')
    cover = ward.load('twelve-point-kernel-upper/verify_cover.py', 'ward_radii_cover')
    prior_bound = ward.load('twelve-point-kernel-upper/verify_all.py', 'ward_radii_witness')
    upper_report = json.loads((ROOT/'twelve-point-ward-upper/verification_report.json').read_bytes())
    data = json.loads((ROOT/'twelve-point-ward-upper/kernel.json').read_bytes())
    B = ward.hex_rational(data['B'])
    need(upper_report['passed'] and upper_report['full_continuous_cover_replayed'] and
         upper_report['cover']['pending'] == 0 and F(upper_report['cover']['c']) == ward.SLACK and
         F(upper_report['mass_eleven_logarithmic_global_upper']) == B+66*ward.SLACK and
         upper_report['exact_positive_definite_Gram_blocks'] == 12 and
         upper_report['combined_kernel_sum_nonnegative_at_complete_stationary_configurations'],
         'fixed previously accepted Ward theorem')
    V11, _ = prior_bound.witness()
    budget = B+66*ward.SLACK-prior_bound.log_lower(V11)
    need(0 < budget < 66*ward.SLACK < 11*LEVEL, 'one excluded row exceeds the global deficit')
    single_radius_cut = load('twelve-point-removal-domain/verify_all.py', 'ward_radii_single_cut')
    cutoff = single_radius_cut.verify_single_radius(V11)
    need(F(cutoff['strict_radius_upper']) <= EXCLUDED[-1][1], 'analytic radius cap closes the tail')
    h, shape = ward.combined_polynomial(data, ward.ward_matrices(data), cover)
    counts = replay_rows(h, shape, cover, ward.SLACK)
    radial = population(data, budget, cover, ward)
    report = {
        'passed': True, 'n': 12, 'mass': 11, 'scope': 'necessary conditions at all global maximizers',
        'previous_ward_cover_used_as_fixed_theorem': True,
        'previous_ward_cover_replayed_by_this_script': False,
        'slack_budget_upper': str(66*ward.SLACK), 'excluded_row_slack_lower': str(11*LEVEL),
        'row_cut_cover': counts, 'analytic_radius_upper': cutoff['strict_radius_upper'],
        'squared_radius_open_bands': [[str(x) for x in band] for band in BANDS],
        'population': radial, 'sharp_global_value_proved': False,
    }
    rendered = json.dumps(report, indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_bytes(rendered.encode('utf-8'))
    else:
        need(json.loads((HERE/'verification_report.json').read_bytes()) == report,
             'fresh replay agrees with the saved report')
    print(rendered, end='')


if __name__ == '__main__':
    main()
