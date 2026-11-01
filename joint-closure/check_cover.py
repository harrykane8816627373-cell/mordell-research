"""Replay every partition, contraction and exclusion using integer intervals only."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import gzip
import json
import interval_engine as E


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def check(path, compact_output=None):
    raw = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
    data = json.loads(raw)
    E.setup(data['precision_bits'])
    settings = E.model(data['settings']['scope'])
    require(data['settings'] == json.loads(json.dumps(settings)), 'unchanged rigorous model')
    pending = {}
    for index, box in data['pending']:
        require(index not in pending, 'unique unresolved leaf')
        pending[index] = [tuple(v) for v in box]
    require(data['complete_cover'] is (not pending)
            and data['full_eight_point_exclusion_proved'] is (not pending),
            'unresolved regions forbid a full-exclusion claim')
    nodes = data['nodes']
    stack = [(0, settings['initial_box'])]
    seen = set()
    counts = Counter()
    newtons = 0
    removed_noop_preconditioners = 0
    while stack:
        index, input_box = stack.pop()
        require(type(index) is int and 0 <= index < len(nodes) and index not in seen,
                'partition is a rooted tree')
        seen.add(index)
        record = nodes[index]
        if record['kind'] == 'unknown':
            require(index in pending and pending[index] == input_box,
                    'every unresolved region remains explicitly present')
            counts['unresolved'] += 1
            continue
        require(index not in pending, 'processed node is not unresolved')
        passes = record.get('passes', 5)
        require(type(passes) is int and 1 <= passes <= 10, 'bounded exact propagation passes')
        box, reason, aux = E.contract(input_box, settings, passes=passes)
        if 'newton' in record:
            R = record['newton']
            require(box is not None and len(R) == 13 and all(len(row) == 13 for row in R),
                    'valid interval-Newton preconditioner')
            require(all(type(v) is int and abs(v) < 1 << 100 for row in R for v in row),
                    'exact rational preconditioner entries')
            try:
                before = box
                after, _ = E.newton_contract(box, R, settings['radius'],
                    record.get('newton_type', 'polynomial'),
                    aux[0][7] if record.get('use_force_slab', False) else None)
                box, reason, aux = E.contract(after, settings, passes=passes)
                if compact_output is not None and after == before:
                    for key in ['newton', 'newton_type', 'use_force_slab']:
                        record.pop(key, None)
                    if box != before:
                        record['repeat_basic'] = True
                    removed_noop_preconditioners += 1
            except ValueError:
                box, reason, aux = None, 'interval_newton', None
            newtons += 1
        elif record.get('repeat_basic', False):
            require(box is not None, 'repeat propagation starts on a live box')
            box, reason, aux = E.contract(box, settings, passes=passes)
        if 'poly_rayleigh' in record:
            c = record['poly_rayleigh']
            require(box is not None and len(c) == 8 and all(len(v) == 2 for v in c),
                    'valid polynomial full-Hessian coefficient vector')
            require(all(type(v) is int for z in c for v in z) and any(v for z in c for v in z),
                    'nonzero exact Gaussian-integer polynomial coefficients')
            require(E.polynomial_rayleigh(aux[0], c)[1] < 0,
                    'strict polynomial full-Hessian negative curvature')
            box, reason = None, 'full_polynomial_Hessian'
        if 'rayleigh' in record:
            h = record['rayleigh']
            require(box is not None and len(h) == 8 and all(len(v) == 2 for v in h),
                    'valid full-Hessian direction')
            require(all(type(v) is int for z in h for v in z) and any(v for z in h for v in z),
                    'nonzero exact Gaussian-integer direction')
            require(E.rayleigh(aux[0], aux[1], h)[1] < 0, 'strict full-Hessian negative curvature')
            box, reason = None, 'full_Hessian_Rayleigh'
        if record['kind'] == 'excluded':
            require(box is None and reason == record['reason'], 'entire leaf rigorously excluded')
            counts[reason] += 1
        elif record['kind'] == 'split':
            require(box is not None, 'split follows sound contractions')
            dim, cut = record['dim'], record['cut']
            require(type(dim) is int and 0 <= dim < 13 and type(cut) is int
                    and box[dim][0] < cut < box[dim][1], 'proper exact partition')
            left, right = list(box), list(box)
            left[dim], right[dim] = (box[dim][0], cut), (cut, box[dim][1])
            stack += [(record['right'], right), (record['left'], left)]
            counts['split'] += 1
        else:
            raise ArithmeticError('known certificate node kind')
    require(len(seen) == len(nodes), 'no unreachable certificate nodes')
    require(counts['unresolved'] == len(pending), 'all unresolved leaves accounted for')
    complete = counts['unresolved'] == 0
    require(data['complete_cover'] is complete, 'honest coverage status')
    # A partial cover is never accepted as an empty global feasible set.
    require(data['full_eight_point_exclusion_proved'] is complete, 'no incomplete full-proof claim')
    if compact_output is not None:
        for key in ['seconds', 'processed_this_run', 'counts_this_run']:
            data.pop(key, None)
        data['compaction'] = {'removed_exactly_noop_preconditioners': removed_noop_preconditioners}
        encoded = (json.dumps(data, separators=(',', ':'))+'\n').encode('utf-8')
        compact_output.write_bytes(gzip.compress(encoded, compresslevel=9, mtime=0))
    return {'passed': True, 'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'precision_bits': E.BITS, 'scope': settings['scope'], 'tree_nodes': len(nodes),
            'counts': dict(counts), 'interval_Newton_enclosures_replayed': newtons,
            'removed_exactly_noop_preconditioners': removed_noop_preconditioners,
            'all_partitions_and_contractions_replayed': True,
            'complete_cover': complete, 'full_eight_point_exclusion_proved': complete}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--compact-output', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.certificate, args.compact_output), indent=2))
