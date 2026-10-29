"""Optional NumPy exploration; only verify_witnesses.py certifies results.

Deterministic seeds and a bounded L-BFGS search, without symmetry constraints.
Floating termination means neither exact stationarity nor global optimality.
Saved integer coordinates are independently accepted with Python integers.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def objective(x):
    z = x.reshape(-1, 2)
    n = len(z)
    diff = z[:, None, :]-z[None, :, :]
    d2 = np.sum(diff*diff, axis=2)
    np.fill_diagonal(d2, 1.0)
    if np.min(d2) <= 0:
        return math.inf, np.zeros_like(x)
    f = (n-1)/2*np.sum(z*z)-np.sum(np.log(d2))/2
    g = (n-1)*z-2*np.sum(diff/d2[:, :, None], axis=1)
    return float(f), g.ravel()


def optimize(seed, max_steps=1200):
    x = np.asarray(seed, dtype=float).ravel()
    f, g = objective(x)
    history = []
    status = 'iteration_budget'
    for step in range(max_steps):
        if np.max(np.abs(g)) < 2e-8:
            status = 'floating_gradient_tolerance'
            break
        q = g.copy()
        alphas = []
        for s, y, rho in reversed(history):
            a = rho*np.dot(s, q)
            alphas.append(a)
            q -= a*y
        scale = (np.dot(history[-1][0], history[-1][1]) /
                 np.dot(history[-1][1], history[-1][1])) if history else 1/(len(x)//2)
        p = scale*q
        for (s, y, rho), a in zip(history, reversed(alphas)):
            p += s*(a-rho*np.dot(y, p))
        p = -p
        slope = float(np.dot(g, p))
        if slope >= 0:
            history.clear()
            p = -g/(len(x)//2)
            slope = float(np.dot(g, p))
        alpha = 1.0
        for backtrack in range(40):
            trial = x+alpha*p
            nf, ng = objective(trial)
            if nf <= f+1e-4*alpha*slope:
                break
            alpha *= 0.5
        else:
            status = 'floating_line_search_stop'
            break
        s, y = trial-x, ng-g
        sy = float(np.dot(s, y))
        if sy > 1e-14*np.linalg.norm(s)*np.linalg.norm(y):
            history.append((s, y, 1/sy))
            history = history[-15:]
        x, f, g = trial, nf, ng
    z = x.reshape(-1, 2)
    z -= np.mean(z, axis=0)
    z *= math.sqrt(len(z)/np.sum(z*z))
    return z, {'objective': f, 'gradient_max': float(np.max(np.abs(g))),
               'steps': step+1, 'status': status}


def regular_seed(row):
    n = row['n']
    points = [[0.0, 0.0]] if row['center'] else []
    for (m, b), sign in zip(row['rings'], row['signs']):
        r = math.sqrt(b/(n-1))
        angles = (2*np.pi*np.arange(m)+(np.pi if sign < 0 else 0))/m
        points.extend(np.column_stack((r*np.cos(angles), r*np.sin(angles))).tolist())
    return np.array(points)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--min-n', type=int, default=10)
    parser.add_argument('--max-n', type=int, default=100)
    parser.add_argument('--steps', type=int, default=1200)
    parser.add_argument('--output', type=Path, default=HERE/'integer_witnesses.json')
    args = parser.parse_args()
    if not 10 <= args.min_n <= args.max_n <= 100 or args.steps < 1:
        raise ValueError('require 10<=min<=max<=100 and positive steps')
    previous = json.loads((ROOT/'radial-angular/population_search_verification.json').read_text(encoding='utf-8'))
    rows = {r['n']: r for r in previous['rows']}
    rng = np.random.default_rng(20261028)
    output = []
    last = None
    for n in range(args.min_n, args.max_n+1):
        base = regular_seed(rows[n])
        seeds = [('regular_rings', base), ('perturbed_rings', base+0.04*rng.normal(size=base.shape))]
        if last is not None:
            angle = n*2.399963229728653
            extended = np.vstack((last, [1.45*math.cos(angle), 1.45*math.sin(angle)]))
            seeds.append(('previous_n_plus_point', extended))
        best = None
        runs = []
        for label, seed in seeds:
            points, info = optimize(seed, args.steps)
            info['seed'] = label
            runs.append(info)
            if best is None or info['objective'] < best[0]:
                best = (info['objective'], points, label)
        last = best[1]
        integer_points = np.rint(last*10**8).astype(np.int64).tolist()
        output.append({'n': n, 'integer_points': integer_points,
                       'selected_seed': best[2], 'exploratory_runs': runs})
        args.output.write_text(json.dumps({'date': '2026-10-28',
            'scope': 'Floating discovery only. Integer coordinates require independent exact acceptance.',
            'random_seed': 20261028, 'grid_scale': 10**8, 'rows': output},
            ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print(n, best[2], 'log_value', n*(n-1)/2-best[0], flush=True)


if __name__ == '__main__':
    main()
