"""Exploration only: floating stationary scouts, never a completeness certificate."""
from pathlib import Path
import json
import math
import numpy as np

OUT = Path(__file__).resolve().parent/'scout.json'
GEN = np.random.default_rng(20261101)
KEEP = np.array([*range(8), *range(9, 16)])


def objects(q):
    delta = q[:, None]-q[None, :]
    np.fill_diagonal(delta, 1)
    inv = 1/delta
    np.fill_diagonal(inv, 0)
    force = inv.sum(axis=1)
    g = 4*q-force.conjugate()
    b = -(inv*inv)
    np.fill_diagonal(b, -b.sum(axis=1))
    h = np.block([[4*np.eye(8)+b.real, -b.imag], [-b.imag, 4*np.eye(8)-b.real]])
    return g, h


def fix_gauge(q):
    pivot = np.argmax(np.abs(q))
    q = q.copy()
    q[[0, pivot]] = q[[pivot, 0]]
    q *= np.exp(-1j*np.angle(q[0]))
    return q


def root(q):
    q = fix_gauge(q)
    for _ in range(100):
        q = fix_gauge(q)
        g, h = objects(q)
        fullg = np.r_[g.real, g.imag]
        norm = np.linalg.norm(fullg)
        if norm < 2e-11:
            return q
        if not np.isfinite(norm) or np.max(np.abs(q)) > 10:
            return None
        matrix = h[np.ix_(KEEP, KEEP)]
        try:
            step = np.linalg.solve(matrix, fullg[KEEP])
        except np.linalg.LinAlgError:
            return None
        sn = np.linalg.norm(step)
        if sn > 1:
            step /= sn
        vector = np.zeros(16)
        vector[KEEP] = step
        direction = vector[:8]+1j*vector[8:]
        accepted = False
        for power in range(15):
            new = q-direction/(2**power)
            try:
                gn, _ = objects(new)
            except FloatingPointError:
                continue
            if np.linalg.norm(gn) < norm*(1-1e-5/(2**power)):
                q = new
                accepted = True
                break
        if not accepted:
            return None
    return None


roots = {}
successes = 0
for trial in range(3000):
    q = GEN.normal(size=8)+1j*GEN.normal(size=8)
    q -= q.mean()
    q *= np.sqrt(7/np.vdot(q, q).real)
    q = root(q)
    if q is not None:
        successes += 1
        pairs = np.array([abs(q[i]-q[j])**2 for i in range(8) for j in range(i)])
        key = tuple(np.round(np.sort(pairs), 7))
        if key not in roots:
            g, h = objects(q)
            moments = [sum(q**k) for k in range(2, 9)]
            eps = sum(abs(q)**2*(abs(q)**2-1)**2)
            vals, vectors = np.linalg.eigh(h)
            roots[key] = {'q': [[float(z.real), float(z.imag)] for z in q],
                'force_residual': float(np.linalg.norm(g)),
                'log_Delta_over_target': float(np.log(pairs).sum()-7*np.log(7)),
                'E': float(eps), 'ellipse': float(abs(moments[0])**2+54/175*abs(moments[1])**2),
                's2_abs': float(abs(moments[0])), 's3_abs': float(abs(moments[1])),
                'minimum_H_eigenvalue': float(vals[0]), 'H_eigenvalues': vals.tolist(),
                'negative_direction': vectors[:, 0].tolist()}
    if trial % 250 == 249:
        print('trials', trial+1, 'converged', successes, 'unlabeled distance fingerprints', len(roots), flush=True)

rows = sorted(roots.values(), key=lambda r: -r['log_Delta_over_target'])
OUT.write_text(json.dumps({'exploration_only': True, 'complete_stationary_catalogue': False,
    'global_exclusion_proved': False, 'trials': 3000, 'convergences': successes, 'roots': rows}, indent=2)+'\n')
print('TOP ROOTS', flush=True)
for i, row in enumerate(rows[:30]):
    print(i, {k: round(row[k], 8) for k in ['log_Delta_over_target', 'E', 'ellipse',
                                          's2_abs', 's3_abs', 'minimum_H_eigenvalue']}, flush=True)
