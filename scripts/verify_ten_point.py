"""Replay the ten-point exact proof chain in dependency order."""
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
SUITES = ('maximizer-domain', 'ten-point-collapse-domain', 'ten-point-ward-upper',
          'ten-point-ward-radii', 'ten-point-ward-population', 'ten-point-ward-inner-count',
          'ten-point-counted-radii', 'ten-point-tail-exclusion', 'ten-point-banded-kernel',
          'ten-point-inner-pair', 'ten-point-contact-geometry',
          'ten-point-local-basin', 'ten-point-pair-midpoints')


def main():
    for suite in SUITES:
        process = subprocess.run([sys.executable, '-B', '-S', '-O', str(ROOT/suite/'verify_all.py')],
                                 cwd=ROOT, capture_output=True, encoding='utf-8', timeout=1200)
        if process.returncode:
            raise ArithmeticError(suite + ': ' + process.stderr[-3000:] + process.stdout[-2000:])
        report = json.loads(process.stdout)
        if report['passed'] is not True:
            raise ArithmeticError(suite + ': acceptance flag')
        saved = json.loads((ROOT/suite/'verification_report.json').read_bytes())
        if report != saved:
            raise ArithmeticError(suite + ': replay differs from the recorded report')
        print('PASS '+suite+'/verify_all.py', flush=True)
    print('PASS n10 strict global ratio, radial restrictions, exact inner count, and local sharpness.', flush=True)


if __name__ == '__main__':
    main()
