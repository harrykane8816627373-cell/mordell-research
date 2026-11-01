"""Propose an integer-checkable partial cover. Unknown leaves are never discarded."""
from pathlib import Path
from collections import Counter
import argparse
import json
import gzip
import time
import numpy as np
import interval_engine as E


def midpoint_objects(box, midpoint=None):
    midpoint = midpoint or [(a+b)//2 for a, b in box]
    qp = E.positions([E.point(v) for v in midpoint])
    q = np.array([complex(z[0][0]/E.S, z[1][0]/E.S) for z in qp])
    delta = q[:, None]-q[None, :]
    np.fill_diagonal(delta, 1)
    if min(abs(delta[i, j]) for i, j in E.PAIRS) < 1e-5:
        return None
    inv = 1/delta
    np.fill_diagonal(inv, 0)
    b = -(inv*inv)
    np.fill_diagonal(b, -b.sum(axis=1))
    h = np.block([[4*np.eye(8)+b.real, -b.imag], [-b.imag, 4*np.eye(8)-b.real]])
    L = np.array([[complex(*z) for z in row] for row in E.L])
    realL = np.r_[L.real, L.imag]
    return h, realL.T@h@realL


def numeric_columns(midpoint):
    qp = E.positions([E.point(v) for v in midpoint])
    q = np.array([complex(z[0][0]/E.S, z[1][0]/E.S) for z in qp])
    V = np.vander(q, 9, increasing=True)
    sums = V.sum(axis=0)
    sums[1] = 0
    F = np.zeros((8, 8), dtype=complex)
    F[:, 1] = 8*q.conjugate()
    t = abs(q)**2
    for k in range(2, 8):
        F[:, k] = 8*k*t*V[:, k-2]+k*(k-1)*V[:, k-2]
        for r in range(k-1):
            F[:, k] -= 2*(k-1-r)*sums[r]*V[:, k-2-r]
    return V[:, :8], F, sums[2:8]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--nodes', type=int, default=5000)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--scope', default='global_maximum')
    p.add_argument('--resume', type=Path)
    args = p.parse_args()
    E.setup(40)
    settings = E.model(args.scope)
    if args.resume:
        raw = gzip.decompress(args.resume.read_bytes()) if args.resume.suffix == '.gz' else args.resume.read_bytes()
        old = json.loads(raw)
        nodes, pending = old['nodes'], [(i, [tuple(v) for v in box]) for i, box in old['pending']]
        if old['settings'] != json.loads(json.dumps(settings)):
            raise ValueError('resume model mismatch')
    else:
        nodes = [{'kind': 'unknown'}]
        pending = [(0, settings['initial_box'])]
    counts = Counter()
    start = time.monotonic()
    processed = 0
    while pending and processed < args.nodes:
        index, input_box = pending.pop()
        record = {'passes': 3}
        box, reason, aux = E.contract(input_box, settings, passes=3)
        if box is not None and max(b-a for a, b in box) < E.S//2:
            try:
                midpoint, qbox = E.slab_center(box, settings['radius'], aux[0][7])
                try:
                    E.invert_edges(qbox)
                    kind = 'force'
                    J = midpoint_objects(box, midpoint)[1]
                except ZeroDivisionError:
                    kind = 'polynomial'
                    Jexact = E.polynomial_jacobian(E.positions([E.point(v) for v in midpoint]))
                    J = np.array([[(lo+hi)/(2*E.S) for lo, hi in row] for row in Jexact])
                if np.linalg.cond(J) < 1e9:
                    R = np.linalg.inv(J)
                    R = [[int(round(v*E.S)) for v in row] for row in R]
                    # The proposer floats choose R only; the enclosure is integer arithmetic.
                    record['newton'] = R
                    record['newton_type'] = kind
                    record['use_force_slab'] = True
                    newbox, _ = E.newton_contract(box, R, settings['radius'], kind, aux[0][7])
                    if any(nhi-nlo < (hi-lo)*.9 for (lo, hi), (nlo, nhi) in zip(box, newbox)):
                        box, reason, aux = E.contract(newbox, settings, passes=3)
                    else:
                        record.pop('newton', None)
                        record.pop('newton_type', None)
                        record.pop('use_force_slab', None)
            except (ZeroDivisionError, np.linalg.LinAlgError):
                record.pop('newton', None)
            except ValueError:
                if 'newton' not in record:
                    raise
                box, reason, aux = None, 'interval_newton', None
        if box is not None and max(b-a for a, b in box) < E.S:
            midpoint, _ = E.slab_center(box, settings['radius'], aux[0][7])
            V, F, sv = numeric_columns(midpoint)
            M = 64*(V[:, 2:].conj().T@V[:, 2:])-8*np.outer(sv.conj(), sv)-F[:, 2:].conj().T@F[:, 2:]
            values, vectors = np.linalg.eigh((M+M.conj().T)/2)
            if values[0] < -.1:
                coeff = np.zeros(8, dtype=complex)
                coeff[2:] = vectors[:, 0]
                coeff[0] = -sv@coeff[2:]/8
                coeff = [[int(round(z.real*1024)), int(round(z.imag*1024))] for z in coeff]
                if E.polynomial_rayleigh(aux[0], coeff)[1] < 0:
                    record['poly_rayleigh'] = coeff
                    box, reason = None, 'full_polynomial_Hessian'
        if box is not None and max(b-a for a, b in box) < E.S//2:
            center = midpoint_objects(box)
            if center is not None:
                values, vectors = np.linalg.eigh(center[0])
                if values[0] < -.1:
                    vector = vectors[:, 0]
                    direction = [[int(round(vector[i]*1024)), int(round(vector[8+i]*1024))]
                                 for i in range(8)]
                    if E.rayleigh(aux[0], aux[1], direction)[1] < 0:
                        record['rayleigh'] = direction
                        box, reason = None, 'full_Hessian_Rayleigh'
        if box is None:
            record.update(kind='excluded', reason=reason)
            counts[reason] += 1
        else:
            dim = max(range(13), key=lambda k: box[k][1]-box[k][0])
            lo, hi = box[dim]
            cut = (lo+hi)//2
            if cut <= lo or cut >= hi:
                pending.insert(0, (index, input_box))
                nodes[index] = {'kind': 'unknown', 'reason': 'mesh_limit'}
                break
            left, right = list(box), list(box)
            left[dim], right[dim] = (lo, cut), (cut, hi)
            li, ri = len(nodes), len(nodes)+1
            nodes += [{'kind': 'unknown'}, {'kind': 'unknown'}]
            record.update(kind='split', dim=dim, cut=cut, left=li, right=ri)
            pending += [(ri, right), (li, left)]
            counts['split'] += 1
        nodes[index] = record
        processed += 1
        if processed % 10000 == 0:
            print('processed', processed, 'pending', len(pending), 'seconds', round(time.monotonic()-start, 1),
                  'counts', dict(counts), flush=True)
    report = {'precision_bits': E.BITS, 'settings': settings,
              'nodes': nodes, 'pending': pending, 'processed_this_run': processed,
              'counts_this_run': dict(counts), 'complete_cover': not pending,
              'full_eight_point_exclusion_proved': not pending,
              'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps({k: report[k] for k in ['processed_this_run', 'counts_this_run',
                     'complete_cover', 'full_eight_point_exclusion_proved', 'seconds']}, indent=2))


if __name__ == '__main__':
    main()
