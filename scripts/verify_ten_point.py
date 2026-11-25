"""Replay the ten-point exact proof chain in dependency order."""
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
SUITES = ('maximizer-domain', 'ten-point-collapse-domain', 'ten-point-ward-upper',
          'ten-point-ward-radii', 'ten-point-ward-population', 'ten-point-ward-inner-count',
          'ten-point-counted-radii', 'ten-point-tail-exclusion', 'ten-point-banded-kernel',
          'ten-point-inner-pair', 'ten-point-banded-kernel-refined', 'ten-point-curved-kernel',
          'ten-point-refined-radii', 'ten-point-conditioned-inner', 'ten-point-radial-dispersion',
          'ten-point-outer-angles', 'ten-point-outer-moments', 'ten-point-contact-geometry',
          'ten-point-local-basin', 'ten-point-pair-midpoints', 'ten-point-tight-pairing',
          'ten-point-angular-upper', 'ten-point-midpoint-budget', 'ten-point-weighted-pair-basin', 'ten-point-d2-cover',
          'ten-point-saddle-exclusion', 'ten-point-angular-radii', 'ten-point-outer-gap',
          'ten-point-inner-midpoint', 'ten-point-inner-edge', 'ten-point-quartic-midpoints',
          'ten-point-half-turn-chart', 'ten-point-half-turn-dispersion', 'ten-point-phase-sweep',
          'ten-point-inner-field', 'ten-point-outer-recovery', 'ten-point-pair-perturbation')


def replay(suite,entry='verify_all.py',record='verification_report.json'):
    process = subprocess.run([sys.executable, '-B', '-S', '-O', str(ROOT/suite/entry)],
                             cwd=ROOT, capture_output=True, encoding='utf-8', timeout=1200)
    if process.returncode:
        raise ArithmeticError(suite + ': ' + process.stderr[-3000:] + process.stdout[-2000:])
    report = json.loads(process.stdout)
    if report['passed'] is not True:
        raise ArithmeticError(suite + ': acceptance flag')
    saved = json.loads((ROOT/suite/record).read_bytes())
    if report != saved:
        raise ArithmeticError(suite + ': replay differs from the recorded report')
    print('PASS '+suite+'/'+entry, flush=True)


def main():
    for suite in SUITES:
        replay(suite)
    replay('ten-point-sharp-contacts','verify_sharp.py','global_verification_report.json')
    print('PASS n10 sharp global value and unique maximizing orbit.', flush=True)


if __name__ == '__main__':
    main()
