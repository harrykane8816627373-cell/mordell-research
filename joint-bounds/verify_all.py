"""Verify fixed bytes and replay the new exact certificates in isolation.

Run from a repository or offline-package root:
    python -B -S -O joint-bounds/verify_all.py
No network, installation, or nonstandard dependency is required.
"""
from pathlib import Path
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
    started = time.monotonic()
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    paths = [row['path'] for row in manifest['files']]
    require(len(paths) == len(set(paths)), 'unique manifest paths')
    for row in manifest['files']:
        rel = Path(row['path'])
        require(not rel.is_absolute() and '..' not in rel.parts, 'safe manifest relative path')
        data = (ROOT/rel).read_bytes()
        require(len(data) == row['bytes'], 'byte size '+row['path'])
        require(hashlib.sha256(data).hexdigest() == row['sha256'], 'sha256 '+row['path'])
    entries = [('verify_angular_bounds.py', [], 'angular_verification.json'),
               ('verify_phase_search.py', [], 'phase_search_verification.json'),
               ('verify_subset_moment.py', ['--write-report'], 'subset_moment_comparison.json')]
    runs = []
    with tempfile.TemporaryDirectory(prefix='mordell-joint-') as td:
        stage = Path(td)
        for rel in paths:
            target = stage/rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/rel, target)
        for name, args, report_name in entries:
            then = time.monotonic()
            proc = subprocess.run([sys.executable, '-X', 'int_max_str_digits=0', '-B', '-S', '-O', name, *args],
                                  cwd=stage/'joint-bounds', capture_output=True,
                                  text=True, encoding='utf-8', timeout=120)
            require(proc.returncode == 0, name+' failed: '+proc.stdout[-1500:]+proc.stderr[-1500:])
            got = json.loads((stage/'joint-bounds'/report_name).read_text(encoding='utf-8'))
            expected = json.loads((HERE/report_name).read_text(encoding='utf-8'))
            require(got == expected, 'deterministic report replay '+report_name)
            require(got.get('passed') is True or got.get('status') == 'PASS', 'PASS '+report_name)
            runs.append({'entry': name, 'passed': True, 'seconds': round(time.monotonic()-then, 3)})
            print('PASS '+name, flush=True)
    phase = json.loads((HERE/'phase_search_verification.json').read_text(encoding='utf-8'))
    improved = phase['strict_complete_previous_envelope_improvements']
    require(len(improved) == 45, '45 complete-envelope improvements')
    row91 = next(row for row in phase['rows'] if row['n'] == 91)
    require(row91['gain_basis_points_floor_over_complete_previous'] >= 2483, 'n91 certified 24.83 percent')
    report = {'passed': True, 'date': '2026-10-28', 'manifest_inputs': len(paths), 'runs': runs,
              'strict_global_disk_correction': 'All n>=2 by analytic proof; finite consistency checked n=2..100.',
              'uniform_continuous_phase_family_gain': 'All n>=13: log gain > 1/10; product gain > 11/10.',
              'strict_lower_envelope_improvement_count': len(improved), 'improved_n': improved,
              'n91_gain_basis_points_floor': row91['gain_basis_points_floor_over_complete_previous'],
              'subset_formula_comparisons': 147,
              'known_limits': 'No full eight-point cover, unknown sharp finite-n values, or linear gamma improvement. Subset recurrence remains above disk in all 147 cases.',
              'elapsed_seconds': round(time.monotonic()-started, 3)}
    (HERE/'verification_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
