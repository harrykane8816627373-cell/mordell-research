"""Replay the continuation forest and account for every old unresolved box."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import argparse
import gzip
import hashlib
import json
import refinement as R

E = R.E
HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'joint-closure/certificates/partial_cover.json.gz'


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def check(path, replay_base=True):
    raw = BASE.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == R.BASE_SHA256, 'frozen base certificate')
    base = json.loads(gzip.decompress(raw))
    if replay_base:
        from check_cover import check as check_base
        result = check_base(BASE)
        need(result['passed'] and result['counts']['unresolved'] == 21, 'base tree replayed')
    data = json.loads(gzip.decompress(path.read_bytes()))
    need(data['base_sha256'] == R.BASE_SHA256, 'same frozen starting frontier')
    E.setup(base['precision_bits'])
    need(data['precision_bits'] == E.BITS, 'same exact grid')
    settings = E.model(base['settings']['scope'])
    need(data['settings'] == json.loads(json.dumps(settings)), 'same physical model')
    old = {i: [tuple(v) for v in box] for i, box in base['pending']}
    roots = data['roots']
    need(len(roots) == len(old) and {i for i, j in roots} == set(old)
        and len({j for i, j in roots}) == len(roots), 'all 21 base regions covered once')
    pending = {}
    for index, box in data['pending']:
        need(index not in pending, 'unique open leaf')
        pending[index] = [tuple(v) for v in box]
    need(data['complete_cover'] is (not pending)
        and data['full_eight_point_exclusion_proved'] is (not pending), 'honest initial status')
    stack = [(j, old[i]) for i, j in roots]
    nodes, seen, counts = data['nodes'], set(), Counter()
    taylor_bounds = []
    taylor_centers = Counter()
    while stack:
        index, original = stack.pop()
        need(type(index) is int and 0 <= index < len(nodes) and index not in seen,
            'forest without shared or unreachable branches')
        seen.add(index)
        node = nodes[index]
        if node['kind'] == 'unknown':
            need(index in pending and pending[index] == original, 'unresolved box remains explicit')
            counts['unresolved'] += 1
            continue
        need(index not in pending, 'processed node is not marked pending')
        box, reason, aux = R.contract(original, settings, passes=5)
        if node['kind'] == 'excluded':
            need(box is None and reason == node['reason'], 'whole leaf excluded with exact bound')
            counts[reason] += 1
            if reason == 'mass_log_Taylor':
                taylor_bounds.append(aux['log_Delta_over_target_upper'])
                taylor_centers[aux['center_source']] += 1
        elif node['kind'] == 'split':
            need(box is not None, 'partition follows sound contractions')
            dim, cut = node['dim'], node['cut']
            need(type(dim) is int and 0 <= dim < 13 and type(cut) is int
                and box[dim][0] < cut < box[dim][1], 'proper exact partition')
            left, right = list(box), list(box)
            left[dim], right[dim] = (box[dim][0], cut), (cut, box[dim][1])
            stack += [(node['left'], left), (node['right'], right)]
            counts['split'] += 1
        else:
            raise ArithmeticError('unknown continuation node')
    need(len(seen) == len(nodes) and counts['unresolved'] == len(pending), 'complete forest accounting')
    complete = not pending
    need(data['complete_cover'] is complete and data['full_eight_point_exclusion_proved'] is complete,
        'partial continuation is not an emptiness proof')
    excluded = sum(v for k, v in counts.items() if k not in ['split', 'unresolved'])
    return {'passed': True, 'base_certificate_sha256': R.BASE_SHA256,
        'continuation_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'precision_bits': E.BITS, 'scope': settings['scope'],
        'base_open_regions_accounted_for': len(old), 'continuation_nodes': len(nodes),
        'continuation_counts': dict(sorted(counts.items())),
        'combined_tree_nodes': 60021-len(old)+len(nodes),
        'combined_excluded_leaves': 29990+excluded,
        'unresolved_regions': len(pending),
        'new_Taylor_excluded_leaves': len(taylor_bounds),
        'Taylor_center_counts': dict(sorted(taylor_centers.items())),
        'largest_Taylor_upper_endpoint': max(taylor_bounds) if taylor_bounds else None,
        'largest_Taylor_log_upper': str(Fraction(max(taylor_bounds), E.S)) if taylor_bounds else None,
        'complete_cover': complete, 'full_eight_point_exclusion_proved': complete}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('certificate', type=Path)
    p.add_argument('--base-already-replayed', action='store_true',
        help='only for a caller that separately replayed the pinned base tree')
    args = p.parse_args()
    print(json.dumps(check(args.certificate, not args.base_already_replayed), indent=2))
