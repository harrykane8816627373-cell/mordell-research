"""Continue the frozen 21-box frontier; acceptance contains no floating point."""
from pathlib import Path
from collections import Counter
import argparse
import gzip
import hashlib
import json
import time
import refinement as R

E = R.E
HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'joint-closure/certificates/partial_cover.json.gz'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes', type=int, default=10000)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--resume', type=Path)
    args = p.parse_args()
    source = BASE.read_bytes()
    if hashlib.sha256(source).hexdigest() != R.BASE_SHA256:
        raise ArithmeticError('frozen starting frontier')
    base = json.loads(gzip.decompress(source))
    E.setup(base['precision_bits'])
    settings = E.model(base['settings']['scope'])
    if args.resume:
        data = json.loads(gzip.decompress(args.resume.read_bytes()))
        nodes = data['nodes']
        roots = data['roots']
        pending = [(i, [tuple(v) for v in box]) for i, box in data['pending']]
    else:
        nodes = [{'kind': 'unknown'} for _ in base['pending']]
        roots = [[i, j] for j, (i, _) in enumerate(base['pending'])]
        pending = [(j, [tuple(v) for v in box]) for j, (_, box) in enumerate(base['pending'])]
    counts = Counter()
    start = time.monotonic()
    for processed in range(args.nodes):
        if not pending:
            break
        index, original = pending.pop()
        box, reason, aux = R.contract(original, settings, passes=5)
        if box is None:
            record = {'kind': 'excluded', 'reason': reason}
            counts[reason] += 1
        else:
            dim = max(range(13), key=lambda k: box[k][1]-box[k][0])
            lo, hi = box[dim]
            cut = (lo+hi)//2
            if cut <= lo or cut >= hi:
                pending.insert(0, (index, original))
                break
            left, right = list(box), list(box)
            left[dim], right[dim] = (lo, cut), (cut, hi)
            li, ri = len(nodes), len(nodes)+1
            nodes += [{'kind': 'unknown'}, {'kind': 'unknown'}]
            record = {'kind': 'split', 'dim': dim, 'cut': cut, 'left': li, 'right': ri}
            pending += [(ri, right), (li, left)]
            counts['split'] += 1
        nodes[index] = record
        if (processed+1) % 2000 == 0:
            print('processed', processed+1, 'pending', len(pending),
                'counts', dict(counts), 'seconds', round(time.monotonic()-start, 2), flush=True)
    data = {'base_sha256': R.BASE_SHA256, 'precision_bits': E.BITS,
        'settings': settings, 'roots': roots, 'nodes': nodes, 'pending': pending,
        'complete_cover': not pending, 'full_eight_point_exclusion_proved': not pending}
    args.output.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),
        compresslevel=9,mtime=0))
    print(json.dumps({'counts': dict(counts), 'pending': len(pending),
        'full_eight_point_exclusion_proved': not pending,
        'seconds': time.monotonic()-start},indent=2))


if __name__ == '__main__':
    main()
