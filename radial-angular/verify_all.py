"""Replay this research supplement with exact fixed inputs in a temporary copy.

From the repository or extracted package root:
    python -B -S -O radial-angular/verify_all.py
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def main():
    began = time.monotonic()
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    paths = [r['path'] for r in manifest['files']]
    require(len(paths) == len(set(paths)), 'unique manifest paths')
    for row in manifest['files']:
        rel = Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts, 'safe relative input path')
        data = (ROOT/rel).read_bytes()
        require(len(data) == row['bytes'], 'byte count '+row['path'])
        require(hashlib.sha256(data).hexdigest() == row['sha256'], 'sha256 '+row['path'])
    entries = [('verify_radial_layers.py', [], 'radial_layers_verification.json'),
               ('verify_mixed_energy.py', [], 'mixed_energy_verification.json'),
               ('verify_population_gap.py', ['--write-report'], 'population_gap_verification.json'),
               ('verify_population_search.py', [], 'population_search_verification.json'),
               ('verify_joint_budget.py', [], 'joint_budget_verification.json')]
    runs = []
    with tempfile.TemporaryDirectory(prefix='mordell-radial-angular-') as td:
        stage = Path(td)
        for rel in paths:
            target = stage/rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/rel, target)
        for name, args, report_name in entries:
            start = time.monotonic()
            proc = subprocess.run([sys.executable, '-X', 'int_max_str_digits=0', '-B', '-S', '-O', name, *args],
                                  cwd=stage/'radial-angular', capture_output=True, text=True,
                                  encoding='utf-8', timeout=120)
            require(proc.returncode == 0, name+' failed: '+proc.stdout[-1200:]+proc.stderr[-1800:])
            got = json.loads((stage/'radial-angular'/report_name).read_text(encoding='utf-8'))
            expected = json.loads((HERE/report_name).read_text(encoding='utf-8'))
            require(got == expected, 'deterministic report replay '+report_name)
            require(got.get('passed') is True or got.get('status') == 'PASS', 'PASS report '+name)
            runs.append({'entry': name, 'passed': True, 'seconds': round(time.monotonic()-start, 3)})
            print('PASS '+name, flush=True)
    populations = json.loads((HERE/'population_search_verification.json').read_text(encoding='utf-8'))
    improved = populations['strict_explicit_complete_previous_envelope_improvements']
    require(improved == sorted(set(range(13, 101))-{17}), '87 complete-envelope improvements')
    require(max(r['continuous_phase_relative_gap_ppm_ceiling'] for r in populations['rows']) <= 1,
            'all 99 selected fixed families have <=1ppm continuous-phase gap')
    row95 = next(r for r in populations['rows'] if r['n'] == 95)
    require(row95['explicit_gain_basis_points_floor'] >= 1204022, 'n95 exact gain')
    layers = json.loads((HERE/'radial_layers_verification.json').read_text(encoding='utf-8'))
    require(layers['finite_layer_constant'] == '3/22' and layers['asymptotic_layer_constant'] == '101/250',
            'general layer constants')
    joint = json.loads((HERE/'joint_budget_verification.json').read_text(encoding='utf-8'))
    hundred = next(r for r in joint['rows'] if r['n'] == 100)
    require(hundred['new_necessary_layers'] >= 4 and F(hundred['simple_h1_three_layer_delta_lower']) == F(9, 4),
            'all global hundred-point optimizers have at least four radii')
    report = {'passed': True, 'date': '2026-10-28', 'manifest_inputs': len(paths), 'runs': runs,
              'finite_layer_lower': 'K > (3/22) sqrt(n), all n>=2',
              'asymptotic_layer_lower': 'liminf K_n/sqrt(n) >= 101/250',
              'hundred_point_global_optimum_minimum_radii': 4,
              'hundred_point_gap_upper': hundred['new_global_optimizer_gap_upper_decimal'],
              'strict_explicit_lower_envelope_improvements': len(improved), 'improved_n': improved,
              'selected_continuous_phase_families_with_gap_at_most_one_ppm': 99,
              'independent_even_compositions_checked': populations['independent_ordered_even_compositions_checked'],
              'scope': 'Necessary configuration restrictions, constructive lower bounds and conditional regional bounds. No general sharp value, full eight-point cover, exact gamma or globally optimal ring model.',
              'elapsed_seconds': round(time.monotonic()-began, 3)}
    (HERE/'verification_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
