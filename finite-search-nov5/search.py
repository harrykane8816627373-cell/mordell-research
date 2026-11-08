"""Deterministic optional NumPy discovery. Exact acceptance is separate."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SPEC = importlib.util.spec_from_file_location('prior_search', ROOT/'symmetry-breaking/explore.py')
PRIOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRIOR)


def normalize(points):
    z = np.array(points, dtype=float).reshape(-1, 2)
    z -= z.mean(axis=0)
    return z*math.sqrt(len(z)/np.sum(z*z))


def log_product(points):
    z = normalize(points)
    d = z[:, None, :]-z[None, :, :]
    values = np.sum(d*d, axis=2)
    return float(np.log(values[np.triu_indices(len(z), 1)]).sum())


def seeds_for(n, baseline, last, rng):
    base = baseline[n]
    seeds = [('baseline', base)]
    for amplitude in (.04, .15, .35):
        for k in range(2):
            seeds.append((f'baseline_noise_{amplitude}_{k}', base+amplitude*rng.normal(size=base.shape)))
    index = np.arange(n)
    radius = np.sqrt(2*(index+.5)/n)
    angle = index*2.399963229728653
    spiral = np.column_stack((radius*np.cos(angle), radius*np.sin(angle)))
    seeds.append(('golden_spiral', spiral))
    for k in range(3):
        radius = np.sqrt(2*rng.uniform(size=n))
        angle = rng.uniform(0, 2*np.pi, n)
        seeds.append((f'uniform_disk_{k}', np.column_stack((radius*np.cos(angle), radius*np.sin(angle)))))
    for k in range(2):
        seeds.append((f'gaussian_{k}', rng.normal(size=(n, 2))))
    if n+1 in baseline:
        larger = baseline[n+1]
        radii = np.sum(larger*larger, axis=1)
        for label, index in [('inner', int(np.argmin(radii))), ('outer', int(np.argmax(radii)))]:
            seeds.append((f'baseline_next_remove_{label}', np.delete(larger, index, axis=0)))
    if last is not None:
        for k in range(4):
            a = n*2.399963229728653+k*np.pi/2
            seeds.append((f'previous_best_add_outer_{k}', np.vstack((last, [1.38*np.cos(a), 1.38*np.sin(a)]))))
    return seeds


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--min-n', type=int, default=10)
    p.add_argument('--max-n', type=int, default=50)
    p.add_argument('--steps', type=int, default=600)
    p.add_argument('--random-seed', type=int, default=20261105)
    p.add_argument('--initial-discovery', type=Path, help='warm-start all same/neighbor-n seeds from a saved discovery wave')
    p.add_argument('--output', type=Path, default=HERE/'discovery.json')
    args = p.parse_args()
    if not 10 <= args.min_n <= args.max_n <= 100 or args.steps < 1:
        raise ValueError('10<=min<=max<=100, positive steps required')
    raw = (ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes()
    old = json.loads(raw)
    baseline = {r['n']: normalize(r['integer_points']) for r in old['rows']}
    seed_baseline = dict(baseline)
    initial_hash = None
    if args.initial_discovery is not None:
        initial_raw = args.initial_discovery.read_bytes()
        initial_hash = hashlib.sha256(initial_raw).hexdigest()
        for row in json.loads(initial_raw)['rows']:
            seed_baseline[row['n']] = normalize(row['integer_points'])
    rng = np.random.default_rng(args.random_seed)
    rows = []
    last = None
    started = time.monotonic()
    for n in range(args.min_n, args.max_n+1):
        old_log = log_product(baseline[n])
        best = seed_baseline[n]
        best_log = log_product(best)
        chosen = 'initial_saved_coordinates' if initial_hash else 'unchanged_baseline'
        runs = []
        for label, seed in seeds_for(n, seed_baseline, last, rng):
            if initial_hash:
                label = label.replace('baseline', 'warm_saved')
            points, info = PRIOR.optimize(normalize(seed), args.steps)
            actual_log = log_product(points)
            info.update(seed=label, normalized_log_product=actual_log)
            runs.append(info)
            if actual_log > best_log:
                best, best_log, chosen = points, actual_log, label
        last = best
        row = {'n': n, 'baseline_log_float': old_log, 'selected_seed': chosen,
               'new_log_float': best_log, 'log_gain_float': best_log-old_log,
               'integer_points': np.rint(best*10**10).astype(np.int64).tolist(), 'runs': runs}
        rows.append(row)
        report = {'date': '2026-11-05', 'scope': 'floating discovery only; use verify_all.py for exact acceptance',
                  'random_seed': args.random_seed, 'steps_per_run': args.steps, 'grid_scale': 10**10,
                  'numpy_version': np.__version__, 'baseline_path': 'symmetry-breaking/integer_witnesses.json',
                  'initial_discovery_sha256': initial_hash,
                  'baseline_sha256': hashlib.sha256(raw).hexdigest(),
                  'optimizer_path': 'symmetry-breaking/explore.py',
                  'optimizer_sha256': hashlib.sha256((ROOT/'symmetry-breaking/explore.py').read_bytes()).hexdigest(),
                  'elapsed_seconds': round(time.monotonic()-started, 3), 'rows': rows}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print(n, 'log_gain', best_log-old_log, chosen, 'runs', len(runs), flush=True)


if __name__ == '__main__':
    main()
