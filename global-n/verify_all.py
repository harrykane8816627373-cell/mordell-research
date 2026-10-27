"""Replay exact finite checks. Analytic proofs remain separate obligations."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

from grid_bounds import (certified_record, rational_value, grid_radius,
                         relative_grid_plan, rounding_exponent,
                         passes_rounded_subset_cuts, subset_floor)
from verify_recursion_variance import run as variance_run

ROOT = Path(__file__).resolve().parent


def require(condition, name):
    if not condition:
        raise ArithmeticError(name)


def independent_value(points):
    """Center coordinate Fractions directly; no n^2/D formula."""
    n = len(points)
    means = [sum(F(p[j]) for p in points) / n for j in range(2)]
    z = [(F(p[0]) - means[0], F(p[1]) - means[1]) for p in points]
    Q = sum(x*x+y*y for x, y in z)
    value = prod((x-u)**2+(y-v)**2 for (x, y), (u, v) in combinations(z, 2))
    return value * (F(n) / Q) ** (n*(n-1)//2)


def main():
    manifest = json.loads((ROOT/'source_manifest.json').read_text(encoding='utf-8-sig'))
    for row in manifest['files']:
        p = ROOT/row['file']
        require(p.parent == ROOT and p.is_file(), 'invalid manifest member')
        data = p.read_bytes()
        require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], 'source hash: '+row['file'])
    stored = json.loads((ROOT/'grid_values.json').read_text(encoding='utf-8-sig'))
    total = 0
    for row in stored['rows']:
        replay = certified_record(row['n'], row['q'])
        require(row == replay, 'complete grid replay mismatch')
        require(F(row['value']) == independent_value(row['witness']), 'independent witness normalization')
        require(row['configurations'] == comb(row['grid_points'], row['n']), 'coverage count')
        total += row['configurations']
    independent_cases = 0
    points = [(a, b) for a in range(-1, 2) for b in range(-1, 2)]
    for n in range(2, 6):
        for block in combinations(points, n):
            require(rational_value(block) == independent_value(block), 'independent rational objective')
            independent_cases += 1
    plans = json.loads((ROOT/'grid_plans.json').read_text(encoding='utf-8-sig'))
    for row in plans['rows']:
        require(relative_grid_plan(row['n'], row['epsilon']) == row, 'stored plan mismatch')
        require(row['search_executed'] is False, 'plan incorrectly marked as executed')
    boundaries = 0
    for n in range(3, 101):
        # Exercise the exact boundary that caught the dyadic overestimate.
        from grid_bounds import ceil_sqrt_ratio
        q = (n+1)*ceil_sqrt_ratio(n-1, 2)
        simple = F(7*n*(n-1), 12*q*q)
        epsilon = simple/(1-simple)
        plan = relative_grid_plan(n, epsilon)
        require(plan['q'] == q, 'exact boundary plan')
        eta, tau = rounding_exponent(n, q)
        require(2*q*q*tau*tau >= n-1 and 0 <= tau < 1, 'rational square root enclosure')
        require(eta <= simple, 'boundary exponent enclosure')
        boundaries += 1
    require(relative_grid_plan(9, F(21, 179))['q'] == 20, 'n=9 boundary regression')
    try:
        certified_record(8, 58, max_configurations=1000)
    except ValueError as e:
        require('above budget' in str(e), 'unexpected refusal')
    else:
        raise ArithmeticError('huge search was not refused before execution')
    for q in range(1, 21):
        square = [(q, 0), (0, q), (-q, 0), (0, -q)]
        require(passes_rounded_subset_cuts(square, q), 'true square round image cut away')
        require(grid_radius(4, q) >= q, 'integer box omitted square')
    # Young's inequality used by subset pruning, checked independently in Q^2.
    young_cases = 0
    for q in range(1, 10):
        for ax, ay, bx, by in ((1,2,3,4), (3,-5,-2,1), (-2,4,7,-3), (0,0,1,-1)):
            lhs = F((ax+bx)**2+(ay+by)**2)
            rhs = (1-F(1,q))*(ax*ax+ay*ay)-(q-1)*(bx*bx+by*by)
            require(lhs >= rhs, 'rounded subset Young bound')
            young_cases += 1
    variance = variance_run()
    with tempfile.TemporaryDirectory(prefix='mordell-grid-check-') as tmp:
        copy = Path(tmp)/'verify_second_order_grid.py'
        shutil.copy2(ROOT/'verify_second_order_grid.py', copy)
        result = subprocess.run([sys.executable, '-B', '-S', '-O', str(copy)], capture_output=True, text=True)
        require(result.returncode == 0, 'second-order replay failed: '+result.stderr)
        second = json.loads((Path(tmp)/'second_order_grid_verification.json').read_text(encoding='utf-8'))
        require(second['passed'], 'second-order report failed')
    report = {'passed': True, 'source_files_checked':len(manifest['files']),
              'grid_configurations_replayed':total,
              'grid_cases':len(stored['rows']),
              'independent_normalization_cases':independent_cases,
              'exact_stopping_boundaries':boundaries,
              'subset_cut_square_roundings':20, 'young_inequality_cases':young_cases,
              'variance':variance, 'second_order':second,
              'scope':'Finite exact checks; all-n theorems require GRID_PROOF.md and RECURRENCE_PROOF.md. Large plans were not executed.'}
    (ROOT/'verification_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
