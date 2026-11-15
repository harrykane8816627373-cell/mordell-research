"""Exact removal of the outermost ten-point radial band."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
UPPER = ROOT / 'ten-point-ward-upper'
BANDS = ((F(0), F(3, 10)), (F(41, 60), F(71, 50)), (F(149, 100), F(1539, 1000)))
ROW_FLOOR = {1: F(12, 1000), 2: F(8, 1000), 3: F(1, 1000)}
SLACK = F(3, 10000)
GAP = F(1, 10**6)
U = BANDS[2][1]
SIX_MASS = F(1986, 625)
GRAM_CERTIFICATES = {
    2: {
        'b': ('0', '.83460', '3.10791', '7.48675', '15.00212', '1.80211'),
        'slope': '-.07584',
        'floors': ('-.00039', '-.00008', '.00912'),
    },
    3: {
        'b': ('0', '1.25733', '4.66682', '11.17105', '22.25368', '2.78446'),
        'slope': '-.11380',
        'floors': ('-.00130', '-.00111', '.01661'),
    },
}


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_bytes())


def verify_inputs():
    manifest = read_json(HERE / 'source_manifest.json')
    for key, base in (('files', HERE), ('dependencies', ROOT)):
        for path, expected in manifest[key].items():
            need(digest(base / path) == expected, 'fixed input ' + path)
    collapse = read_json(ROOT / 'ten-point-collapse-domain/verification_report.json')
    ward = read_json(UPPER / 'verification_report.json')
    radii = read_json(ROOT / 'ten-point-ward-radii/verification_report.json')
    inner = read_json(ROOT / 'ten-point-ward-inner-count/verification_report.json')
    need(all(r['passed'] is True for r in (collapse, ward, radii, inner)), 'accepted prerequisites')
    need(F(collapse['radius_strict_upper']) == U and
         F(collapse['subset_mass_and_population']['6']['strict_uncentered_mass_lower']) == SIX_MASS,
         'fixed collapse theorem conclusions')
    need(9 - 4 * BANDS[2][0] < SIX_MASS, 'four outermost points violate six-point complement mass')
    need(F(ward['c']) == SLACK and ward['PSD_blocks'] == 10 and ward['unresolved'] == 0 and
         ward['kernel_sha256'] == digest(UPPER / 'kernel.json'),
         'accepted global Ward kernel and two-point deficit')
    excluded = [tuple(map(F, cover['excluded_closed_interval'])) for cover in radii['covers']]
    need(radii['applies_to_every_global_maximizer'] is True and
         excluded == [(BANDS[0][1], BANDS[1][0]), (BANDS[1][1], BANDS[2][0])],
         'fixed prior radial exclusions')
    need(inner['normalization_mass'] == 9 and inner['all_global_maximizers'] is True
         and inner['exact_inner_count'] == 2 and inner['outer_count'] == 8 and
         F(inner['inner_radius_strict_upper']) == BANDS[0][1],
         'two inner and eight outer points')


def row_polynomial(poly, multipliers, case, other, radii):
    a, b, force, inner_dual, tail_dual = multipliers
    floor = (ROW_FLOOR[case] / 9
             + inner_dual * (int(other == 0) - F(2, 9))
             + tail_dual * (int(other == 2) - F(case - 1, 9)))
    result = radii.rowpoly(poly, (a, b, force))
    # Replace the upstream row floor by a constant depending on the neighbor band.
    for index, d_coeff in (((1, 0, 0), F(1)), ((0, 1, 0), F(1)), ((0, 0, 1), F(-2))):
        result[index] = result.get(index, F(0)) + d_coeff * (radii.FLOOR - floor)
    return result


def prunable(box, ward, radii):
    (ul, uh), (vl, vh), (tl, th) = box
    if ul >= U or vl >= U or uh + vh <= F(1007, 5000):
        return True
    if any(a <= x <= y <= b for a, b in radii.EXCLUSIONS
           for x, y in ((ul, uh), (vl, vh))):
        return True
    if uh + vh - 2 * tl <= F(203, 500) and uh + vh + 2 * th <= F(1, 2500):
        return True
    return radii.prunable(box, ward.U, ward.DMIN, ward.L)


def replay_tree(ward, radii, packed, shape, poly, box, tree):
    stack = [(box, ward.initial(poly, shape, box), 0)]
    counts = {key: 0 for key in 'AP012'}
    maximum_depth = 0
    for action in tree:
        need(stack and action in counts, 'valid action in nonempty tail tree')
        current_box, control, depth = stack.pop()
        counts[action] += 1
        maximum_depth = max(maximum_depth, depth)
        if action == 'P':
            need(prunable(current_box, ward, radii), 'exact necessary-domain exclusion')
        elif action == 'A':
            need(packed.dominates(control, current_box, F(0)), 'exact tail row Bernstein bound')
        else:
            axis = int(action)
            low, high = radii.split(current_box, axis)
            left, right = packed.split(control, axis)
            stack.extend(((high, right, depth + 1), (low, left, depth + 1)))
    need(not stack and counts['A'] + counts['P'] == 1 + sum(counts[key] for key in '012'),
         'complete tail row partition')
    return counts, maximum_depth


def radial_bound(data, case, radial):
    record = GRAM_CERTIFICATES[case]
    b = tuple(map(F, record['b']))
    a = F(record['slope'])
    floors = tuple(map(F, record['floors']))
    poly, penalty = radial.polynomial(radial.radial_data(data), b, a)
    covers = [radial.strictly_above(poly, band, floor) for band, floor in zip(BANDS, floors)]
    lower = 2 * floors[0] + (8 - case) * floors[1] + case * floors[2] - penalty
    return lower, covers


def verify(progress=False):
    verify_inputs()
    ward = load('tail_upper', UPPER / 'verify_all.py')
    radii = load('tail_radii', ROOT / 'ten-point-ward-radii/verify_all.py')
    radial = load('tail_radial', ROOT / 'ten-point-ward-inner-count/verify_radial.py')
    data = read_json(UPPER / 'kernel.json')
    ward_blocks, forces = ward.corrections(data)
    poly, shape = ward.polynomial(data, ward_blocks, forces)
    need(ward.C == SLACK and ward.U_NEW == U and forces, 'same global force kernel')
    need(ward.log.log_lower(ward.witness()) > ward.hx(data['B']) - GAP, 'strict logarithmic value budget')
    ward.old.selftest()
    packed = ward.Packed(shape, True)
    cases = []
    for case in (1, 2, 3):
        entries = read_json(HERE / f'cover-{case}.json')
        neighbor_bands = BANDS[:2] if case == 1 else BANDS
        cursor = BANDS[2][0]
        counts = {key: 0 for key in 'AP012'}
        maximum_depth = 0
        for slab in entries:
            lo, hi = F(slab['lo']), F(slab['hi'])
            multipliers = tuple(map(F, slab['dual']))
            need(lo == cursor and lo < hi <= U and len(multipliers) == 5
                 and len(slab['trees']) == len(neighbor_bands), 'tail slabs and every neighbor band')
            cursor = hi
            for other, (tree, band) in enumerate(zip(slab['trees'], neighbor_bands)):
                adjusted = row_polynomial(poly, multipliers, case, other, radii)
                box = ((lo, hi), band, (-U, U))
                local, depth = replay_tree(ward, radii, packed, shape, adjusted, box, tree)
                for label, count in local.items():
                    counts[label] += count
                maximum_depth = max(maximum_depth, depth)
            if progress:
                print(f'case {case}, through {cursor}, {sum(counts.values())} nodes', file=sys.stderr)
        need(cursor == U, 'complete closed marked tail interval')
        result = {'tail_count': case, 'marked_row_lower': str(ROW_FLOOR[case]),
                  'slabs': len(entries), 'counts': counts, 'nodes': sum(counts.values()),
                  'maximum_depth': maximum_depth, 'unresolved': 0}
        if case == 1:
            contradiction = ROW_FLOOR[case] - 36 * SLACK
            need(contradiction > GAP, 'single marked row exceeds unordered-pair budget')
            result['budget_upper'] = str(GAP)
        else:
            gram_lower, covers = radial_bound(data, case, radial)
            contradiction = gram_lower + case * ROW_FLOOR[case] - (10 - case) * 9 * SLACK
            need(contradiction > 2 * GAP, 'marked rows and Gram exceed double budget')
            result.update({'budget_upper': str(2 * GAP), 'gram_cover': covers,
                           'gram_lower_decimal_display': float(gram_lower)})
        result['contradiction_lower_decimal_display'] = float(contradiction)
        cases.append(result)
    return {'passed': True, 'normalization_mass': 9, 'all_global_maximizers': True,
            'exact_inner_count': 2, 'outer_count': 8,
            'inner_squared_radius_strict_upper': '3/10',
            'outer_squared_radius_open_interval': ['41/60', '71/50'],
            'tail_excluded_closed_interval': list(map(str, BANDS[2])),
            'tail_counts_excluded_by_six_point_mass': '4 or more',
            'total_row_nodes': sum(case['nodes'] for case in cases),
            'ward_budget_minus_witness_log_strict_upper': str(GAP), 'covers': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    parser.add_argument('--progress', action='store_true')
    args = parser.parse_args()
    result = verify(args.progress)
    report = HERE / 'verification_report.json'
    raw = (json.dumps(result, indent=2) + '\n').encode()
    if args.record_report:
        report.write_bytes(raw)
    else:
        need(result == read_json(report), 'saved exact report matches replay')
    print(raw.decode(), end='')


if __name__ == '__main__':
    main()

