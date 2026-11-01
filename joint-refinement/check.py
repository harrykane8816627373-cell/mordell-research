"""Replay the continuation forest and account for every old unresolved box."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import argparse
import gzip
import hashlib
import json
import sys
import refinement as R

E = R.E
HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'joint-closure/certificates/partial_cover.json.gz'


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def check(path, replay_base=True, progress=False):
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
    algorithms, quadratic_centers, quadratic_schemes = Counter(), Counter(), Counter()
    quadratic_bounds = []
    quadratic_without_local_Taylor = 0
    vertex_bounds,vertex_centers=[],Counter()
    vertex_without_local_Taylor=0
    centered_bounds,centered_centers=[],Counter()
    centered_multiplier_choices=Counter()
    centered_mass_multipliers=[]
    while stack:
        index, original = stack.pop()
        need(type(index) is int and 0 <= index < len(nodes) and index not in seen,
            'forest without shared or unreachable branches')
        seen.add(index)
        if progress and len(seen)%5000==0:
            print('Replayed continuation nodes:',len(seen),'of',len(nodes),flush=True,file=sys.stderr)
        node = nodes[index]
        if node['kind'] == 'unknown':
            need(index in pending and pending[index] == original, 'unresolved box remains explicit')
            counts['unresolved'] += 1
            continue
        need(index not in pending, 'processed node is not marked pending')
        algorithm = node.get('algorithm', 1)
        need(type(algorithm) is int and algorithm in [1,2,3,4], 'supported node algorithm')
        algorithms[algorithm] += 1
        box, reason, aux = R.contract(original, settings, passes=5, algorithm=algorithm)
        if node['kind'] == 'excluded':
            need(box is None and reason == node['reason'], 'whole leaf excluded with exact bound')
            counts[reason] += 1
            if reason == 'mass_log_Taylor':
                taylor_bounds.append(aux['log_Delta_over_target_upper'])
                taylor_centers[aux['center_source']] += 1
            if reason == 'mass_quadratic_log':
                quadratic_bounds.append(aux['log_Delta_over_target_upper'])
                quadratic_centers[aux['center_source']] += 1
                quadratic_schemes[aux['weight_scheme']+'/'+aux['quadratic_box_bound']] += 1
                quadratic_without_local_Taylor += not aux['legacy_midpoint_Taylor_defined']
            if reason == 'mass_quadratic_vertices':
                vertex_bounds.append(aux['log_Delta_over_target_upper'])
                vertex_centers[aux['center_source']] += 1
                vertex_without_local_Taylor += not aux['legacy_midpoint_Taylor_defined']
            if reason == 'mass_centered_quadratic':
                centered_bounds.append(aux['log_Delta_over_target_upper'])
                centered_centers[aux['center_source']]+=1
                centered_multiplier_choices[str(aux['mass_multiplier_choices'])]+=1
                centered_mass_multipliers.append(aux['mass_multiplier'])
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
        'processed_node_algorithms': {str(k):v for k,v in sorted(algorithms.items())},
        'combined_tree_nodes': 60021-len(old)+len(nodes),
        'combined_excluded_leaves': 29990+excluded,
        'unresolved_regions': len(pending),
        'new_Taylor_excluded_leaves': len(taylor_bounds),
        'Taylor_center_counts': dict(sorted(taylor_centers.items())),
        'largest_Taylor_upper_endpoint': max(taylor_bounds) if taylor_bounds else None,
        'largest_Taylor_log_upper': str(Fraction(max(taylor_bounds), E.S)) if taylor_bounds else None,
        'quadratic_excluded_leaves': len(quadratic_bounds),
        'quadratic_center_counts': dict(sorted(quadratic_centers.items())),
        'quadratic_bound_schemes': dict(sorted(quadratic_schemes.items())),
        'quadratic_with_legacy_midpoint_Taylor_unavailable': quadratic_without_local_Taylor,
        'largest_quadratic_log_upper': str(Fraction(max(quadratic_bounds), E.S)) if quadratic_bounds else None,
        'vertex_excluded_leaves':len(vertex_bounds),
        'vertex_center_counts':dict(sorted(vertex_centers.items())),
        'vertex_with_legacy_midpoint_Taylor_unavailable':vertex_without_local_Taylor,
        'largest_vertex_log_upper':str(Fraction(max(vertex_bounds),E.S)) if vertex_bounds else None,
        'centered_excluded_leaves':len(centered_bounds),
        'centered_center_counts':dict(sorted(centered_centers.items())),
        'centered_multiplier_choice_counts':dict(sorted(centered_multiplier_choices.items())),
        'centered_mass_multiplier_range':[str(Fraction(min(centered_mass_multipliers),E.S)),
            str(Fraction(max(centered_mass_multipliers),E.S))] if centered_mass_multipliers else None,
        'largest_centered_log_upper':str(Fraction(max(centered_bounds),E.S)) if centered_bounds else None,
        'complete_cover': complete, 'full_eight_point_exclusion_proved': complete}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('certificate', type=Path)
    p.add_argument('--base-already-replayed', action='store_true',
        help='only for a caller that separately replayed the pinned base tree')
    p.add_argument('--progress',action='store_true',help='write replay progress to stderr')
    args = p.parse_args()
    print(json.dumps(check(args.certificate, not args.base_already_replayed,args.progress), indent=2))
