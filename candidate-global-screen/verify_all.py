"""Replay saved numerical screening statistics and exact candidate symmetries.

The numerical replay is an integrity/summary check, not a new search or a
global upper-bound proof. Only verify_symmetry.py uses exact root intervals.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def summarize(rows, pinned_log, tolerance=1e-7):
    # Kept identical to search.py without importing NumPy or the optimizer.
    valid = [r for r in rows if r['gradient_max'] < 1e-5]
    values = sorted(valid, key=lambda r: -r['normalized_log_value'])
    clusters = []
    for r in values:
        if not clusters or abs(r['normalized_log_value'] - clusters[-1]['highest_log_value']) > tolerance:
            clusters.append({'highest_log_value': r['normalized_log_value'],
                             'lowest_log_value': r['normalized_log_value'], 'count': 0,
                             'strata': {}, 'seed_indices': []})
        c = clusters[-1]
        c['count'] += 1
        c['lowest_log_value'] = min(c['lowest_log_value'], r['normalized_log_value'])
        c['strata'][r['stratum']] = c['strata'].get(r['stratum'], 0) + 1
        c['seed_indices'].append(r['seed_index'])
    strata = {}
    for name in sorted(set(r['stratum'] for r in rows)):
        sample = [r for r in rows if r['stratum'] == name]
        strata[name] = {
            'runs': len(sample),
            'gradient_below_1e_5': sum(r['gradient_max'] < 1e-5 for r in sample),
            'at_pinned_value_within_1e_7': sum(
                abs(r['normalized_log_value'] - pinned_log) <= tolerance and r['gradient_max'] < 1e-5
                for r in sample)}
    return {'cluster_log_tolerance': tolerance, 'gradient_screen_tolerance': 1e-5,
            'pinned_log_value': pinned_log, 'runs': len(rows),
            'gradient_screened_runs': len(valid),
            'best_value_repeat_count': clusters[0]['count'] if clusters else 0,
            'clusters': clusters, 'stratified_frequencies': strata,
            'runnerup_log_gap': clusters[0]['highest_log_value'] - clusters[1]['highest_log_value']
                if len(clusters) > 1 else None,
            'runnerup_gap_scope': 'distinct observed gradient-screened logvalueclusters only; absence is not an exclusion theorem'}


def normalized_log_value(points):
    """Independent stdlib floating cross-check; deliberately not a proof bound."""
    n = len(points)
    cx = sum(x for x, y in points) / n
    cy = sum(y for x, y in points) / n
    z = [(x - cx, y - cy) for x, y in points]
    scale = math.sqrt(n / sum(x*x + y*y for x, y in z))
    z = [(x*scale, y*scale) for x, y in z]
    return sum(math.log((x-a)**2 + (y-b)**2)
               for i, (x, y) in enumerate(z) for a, b in z[i+1:])


def replay_screen():
    old = json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())['rows']
    newer = json.loads((ROOT/'finite-search-nov5/integer_witnesses.json').read_bytes())['rows']
    best = {r['n']: r['integer_points'] for r in old}
    best.update({r['n']: r['integer_points'] for r in newer})
    cert = json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    certified = {r['n']: r for r in cert['rows']}
    geometry = {r['n']: r for r in json.loads((HERE/'geometry_report.json').read_bytes())['rows']}
    table = []
    expected_strata = {'anchors': 5, 'perturbations': 30, 'uniform_disk': 30,
                       'gaussian': 20, 'continuation': 10, 'spirals': 5}
    for n in range(10, 20):
        data = json.loads((HERE/f'data/n{n}.json').read_bytes())
        need(data['n'] == n and data['random_seed'] == 2026110800+n, 'fixed search identity')
        need(data['max_steps_per_run'] == 500, 'fixed optimization budget')
        for name, sha in data['source_sha256'].items():
            need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha,
                 'archived search source ' + name)
        rows = data['runs']
        need(len(rows) == 100 and [r['seed_index'] for r in rows] == list(range(100)),
             'all 100 archived runs in original order')
        need({s: sum(r['stratum'] == s for r in rows) for s in expected_strata} == expected_strata,
             'all fixed independent and seeded strata present')
        for r in rows:
            need(r['stratum'] in expected_strata and 0 <= r['steps'] <= 500,
                 'run stratum and bounded step count')
            need(math.isfinite(r['normalized_log_value']) and
                 math.isfinite(r['gradient_max']) and r['gradient_max'] >= 0,
                 'finite saved floating observations')
            points = r['terminal_integer_points_grid_1e10']
            need(len(points) == n and all(len(p) == 2 and all(type(v) is int for v in p) for p in points),
                 'complete terminal coordinates')
        pinned = data['summary']['pinned_log_value']
        summary = summarize(rows, pinned)
        need(summary == data['summary'], 'saved floating statistics reproduced without optimizer')
        selected = max([pinned] + [r['normalized_log_value'] for r in rows])
        need(selected == data['best_log_value'] and selected-pinned == data['best_log_gain_over_pinned'],
             'saved best log value and gain reproduced')
        if n in certified:
            v = [q/2**cert['center_exponent'] for q in certified[n]['center_numerators']]
            flat = v[:1] + [0.0] + v[1:]
            points = list(zip(flat[::2], flat[1::2]))
        else:
            points = best[n]
        need(abs(normalized_log_value(points)-pinned) < 1e-10,
             'pinned source value independent floating cross-check')
        row = {'n': n, 'runs': 100, 'gradient_screened_runs': summary['gradient_screened_runs'],
               'observed_best_cluster_count': summary['best_value_repeat_count'],
               'observed_value_clusters': len(summary['clusters']),
               'runnerup_log_gap': summary['runnerup_log_gap'],
               'independent_disk_and_gaussian_runs': 50,
               'independent_best_cluster_hits': sum(summary['clusters'][0]['strata'].get(s, 0)
                                                    for s in ('uniform_disk', 'gaussian')),
               'best_log_gain_over_pinned': data['best_log_gain_over_pinned'],
               'pinned_root_classification': certified[n]['classification'] if n in certified else 'not_certified',
               'pinned_min_curvature_float': geometry[n]['sphere_tangent_min_curvature'] if n in geometry else None}
        table.append(row)
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for path, sha in manifest['files'].items():
        need(hashlib.sha256((HERE/path).read_bytes()).hexdigest() == sha, 'fixed own input ' + path)
    for path, sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, 'fixed prior input ' + path)
    table = replay_screen()
    # This is the exact interval part, separate from all floating diagnostics.
    process = subprocess.run([sys.executable, '-B', '-S', '-O', str(HERE/'verify_symmetry.py')],
                             capture_output=True, text=True, check=True)
    symmetry = json.loads(process.stdout)
    need(symmetry['passed'] is True, 'exact candidate symmetry replay accepted')
    report = {'passed': True, 'date': '2026-11-08',
              'scope': 'archived numerical summary integrity plus exact candidate symmetries',
              'search_rerun': False, 'saved_search_runs_checked': 1000,
              'screening_counts_prove_global_optimality': False,
              'global_upper_bound_proved': False,
              'numerical_replay_arithmetic': 'Python standard-library binary float; not interval proof',
              'pinned_log_value_crosscheck_absolute_tolerance': 1e-10,
              'proof_priority': [12, 10, 19], 'screening_table': table,
              'exact_symmetry': symmetry, 'input_sha256': manifest['files'],
              'dependency_sha256': manifest['dependencies']}
    if args.record_report:
        (HERE/'verification_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    else:
        need(report == json.loads((HERE/'verification_report.json').read_bytes()), 'complete saved report reproduced')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
