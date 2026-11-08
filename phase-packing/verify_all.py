"""Independent vertex enumeration, direct coordinates and tamper rejection."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import random
import sys

from packing import (build_certificate, check_certificate, model,
                     packing_check, solve_packing)
from build_certificates import input_paths


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def square_solve(matrix, rhs):
    """Independent Gauss-Jordan implementation for vertex enumeration."""
    size = len(rhs)
    a = [[F(x) for x in row]+[F(y)] for row, y in zip(matrix, rhs)]
    for col in range(size):
        pivot = next((row for row in range(col, size) if a[row][col]), None)
        if pivot is None:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        q = a[col][col]
        a[col] = [x/q for x in a[col]]
        for row in range(size):
            if row != col:
                q = a[row][col]
                a[row] = [x-q*y for x, y in zip(a[row], a[col])]
    return [a[i][-1] for i in range(size)]


def enumerate_vertices(A, c):
    n = len(c)
    if n == 0:
        return F(0), 1, 1
    constraints = [list(row) for row in A]
    rhs = [F(1)]*len(A)
    constraints += [[F(i == j) for j in range(n)] for i in range(n)]
    rhs += [F(0)]*n
    best = None
    tried = feasible = 0
    for subset in combinations(range(len(constraints)), n):
        tried += 1
        x = square_solve([constraints[i] for i in subset], [rhs[i] for i in subset])
        if x is None or any(t < 0 for t in x):
            continue
        if any(sum((a*b for a, b in zip(row, x)), F(0)) > 1 for row in A):
            continue
        feasible += 1
        value = sum((a*b for a, b in zip(c, x)), F(0))
        best = value if best is None else max(best, value)
    need(best is not None, 'bounded packing polytope has no feasible vertex')
    return best, tried, feasible


def direct_squares_value(radii, rotations):
    points = []
    for r, (a, b) in zip(radii, rotations):
        need(a*a+b*b == 1, 'exact unit rotation required')
        x, y = r*a, r*b
        for _ in range(4):
            points.append((x, y))
            x, y = -y, x
    need(sum(x for x, y in points) == sum(y for x, y in points) == 0, 'direct centered points')
    Q = sum((x*x+y*y for x, y in points), F(0))
    D = F(1)
    for (x, y), (a, b) in combinations(points, 2):
        D *= (x-a)**2+(y-b)**2
    n = len(points)
    return D*(F(n)/Q)**(n*(n-1)//2)


def load_inputs(base):
    manifest = json.loads((base/'source_manifest.json').read_text(encoding='utf-8'))
    need(set(manifest) == {'algorithm', 'files'} and manifest['algorithm'] == 'sha256', 'manifest schema')
    actual = {p.relative_to(base).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths(base)}
    need(manifest['files'] == actual, 'input file missing, added or changed')
    names = ('four_squares_triangles.json', 'four_squares_radius_box.json',
             'mixed_twentyfive_triangles.json', 'mixed_twentyfive_cycles.json', 'compatible_fourteen.json')
    certs = {}
    for name in names:
        data = (base/'certificates'/name).read_bytes()
        need(len(data) <= 100000, 'certificate input budget')
        certs[name] = json.loads(data)
    return certs, actual


def verify(base):
    certs, hashes = load_inputs(base)
    total_vertex_subsets = 0
    total_feasible_vertices = 0
    for name, cert in certs.items():
        result = check_certificate(cert)
        ms, rs, edges, candidates, A, c = model(cert['populations'], cert['radius_intervals'], cert['maximum_cycle_length'])
        brute, count, feasible = enumerate_vertices(A, c)
        need(result == brute, 'simplex vs independent complete vertex enumeration: '+name)
        total_vertex_subsets += count
        total_feasible_vertices += feasible
        # Re-running discovery must yield an accepted certificate, but its
        # basis is not prescribed: primal/dual equality is what is trusted.
        regenerated = build_certificate(ms, rs, cert['maximum_cycle_length'], cert['center'])
        need(check_certificate(regenerated) == result, 'independent fixed certificate vs rerun optimizer')

    sq = certs['four_squares_triangles.json']
    box = certs['four_squares_radius_box.json']
    tri = certs['mixed_twentyfive_triangles.json']
    cyc = certs['mixed_twentyfive_cycles.json']
    P, G = F(sq['optimal_log_gap']), F(sq['greedy_log_gap'])
    need(sq['primal'] == ['1/2']*4, 'four-square exact fractional witness')
    need(P > F(189, 100)*G, 'four-square log gap improves by more than 89 percent')
    old_greedy_upper = F(sq['pair_upper'])/(1+G)
    need(F(sq['all_phase_upper']) < F(669, 1000)*old_greedy_upper, 'four-square rational upper improves over greedy')
    need(F(box['all_phase_upper']) < F(7, 10)*old_greedy_upper, 'entire radius box improves even fixed-point old bound')
    need(F(cyc['optimal_log_gap']) > F(144, 100)*F(tri['optimal_log_gap']), 'four-cycle strictly improves optimal triangle packing')
    selected = [(vertices, F(weight)) for vertices, weight in zip(cyc['cycles'], cyc['primal']) if F(weight)]
    need(selected == [([0, 1, 2, 3], F(1))], 'mixed family uses specified odd four-cycle')
    need(F(certs['compatible_fourteen.json']['optimal_log_gap']) == 0, 'compatible family must retain zero gap')

    # No root finding or trig evaluation: 16 radius corners, four exact
    # rotations patterns, all 120 distances multiplied for each 16-point set.
    dirs = ((F(1), F(0)), (F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(-8, 17), F(15, 17)))
    direct_count = 0
    rs = [[F(a), F(b)] for a, b in box['radius_intervals']]
    for radii in product(*rs):
        for shift in range(4):
            rotations = [dirs[(i+shift) % 4] for i in range(4)]
            actual_value = direct_squares_value(radii, rotations)
            need(actual_value <= F(box['all_phase_upper']), 'direct coordinate value violates all-box upper')
            direct_count += 1

    # Additional bounded 0/1 packing LPs exercise zero reduced costs,
    # degeneracy and varying column incidences independently of geometry.
    rng = random.Random(11052026)
    random_count = 0
    for _ in range(16):
        A = [[F(rng.randrange(2)) for _ in range(3)] for _ in range(4)]
        for j in range(3):
            A[j][j] = F(1)
        c = [F(rng.randrange(1, 9), rng.randrange(1, 8)) for _ in range(3)]
        x, y, steps = solve_packing(A, c)
        got = packing_check(A, c, x, y)
        wanted, count, feasible = enumerate_vertices(A, c)
        need(got == wanted, 'simplex disagrees with full random-LP vertex enumeration')
        total_vertex_subsets += count
        total_feasible_vertices += feasible
        random_count += 1

    mutations = []
    def mutate(label, fn):
        c = deepcopy(sq)
        fn(c)
        mutations.append((label, c))
    mutate('over-spent primal edge', lambda c: c['primal'].__setitem__(0, '1'))
    mutate('missing dual payment', lambda c: c.__setitem__('dual', ['0']*6))
    mutate('forged geometric coefficient', lambda c: c['gamma'].__setitem__(0, str(F(c['gamma'][0])+1)))
    mutate('unbound radius change', lambda c: c['radius_intervals'][0].__setitem__(0, '8999/10000'))
    mutate('missing cycle', lambda c: c['cycles'].pop())
    mutate('duplicated cycle', lambda c: c['cycles'].__setitem__(1, c['cycles'][0]))
    mutate('forged objective', lambda c: c.__setitem__('optimal_log_gap', str(P+1)))
    mutate('forged pair product', lambda c: c.__setitem__('pair_upper', '1'))
    mutate('forged final bound', lambda c: c.__setitem__('all_phase_upper', '1'))
    mutate('negative dual weight', lambda c: c['dual'].__setitem__(0, '-1'))
    mutate('float endpoint', lambda c: c['radius_intervals'][0].__setitem__(0, 0.9))
    mutate('overlapping radius intervals', lambda c: c['radius_intervals'][0].__setitem__(1, '1'))
    mutate('unbound edge ordering', lambda c: c['dual_edges'].reverse())
    mutate('scope escalation', lambda c: c.__setitem__('scope', 'unrestricted M_n'))
    mutate('unknown field', lambda c: c.__setitem__('trust_saved_result', True))
    mutate('undefined rational', lambda c: c['primal'].__setitem__(0, '1/0'))
    rejected = []
    for label, cert in mutations:
        try:
            check_certificate(cert)
        except ValueError:
            rejected.append(label)
        else:
            raise ArithmeticError('tampered certificate accepted: '+label)
    scope = 'fixed regular-ring populations and stated radius intervals, all phases; no unrestricted optimum claim'
    return {'passed': True, 'python': sys.version.split()[0], 'optimization_safe': True,
            'fixed_certificates_checked': len(certs),
            'all_vertex_constraint_subsets_examined': total_vertex_subsets,
            'feasible_vertex_records_examined': total_feasible_vertices,
            'additional_independent_packing_LPs': random_count,
            'exact_direct_sixteen_point_coordinate_values': direct_count,
            'distances_per_direct_value': 120,
            'tamper_cases_rejected': rejected,
            'four_square_log_gap_optimum': str(P), 'four_square_greedy_gap': str(G),
            'four_square_new_to_old_rational_upper': str(F(sq['all_phase_upper'])/old_greedy_upper),
            'radius_box_new_to_fixed_old_rational_upper': str(F(box['all_phase_upper'])/old_greedy_upper),
            'mixed_triangle_LP_optimum': tri['optimal_log_gap'],
            'mixed_cycle_LP_optimum': cyc['optimal_log_gap'],
            'optimality_scope': 'optimal fractional packing over every enumerated eligible cycle; not exact physical minimum phase deficit',
            'mathematical_scope': scope,
            'continuum_coverage': 'analytic cycle and endpoint inequalities; finite point tests only check implementation',
            'source_sha256': hashes}


if __name__ == '__main__':
    base = Path(__file__).resolve().parent
    report = verify(base)
    (base/'verification_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({
        'passed': report['passed'],
        'optimization_safe': report['optimization_safe'],
        'fixed_certificates_checked': report['fixed_certificates_checked'],
        'vertex_constraint_subsets_examined': report['all_vertex_constraint_subsets_examined'],
        'direct_coordinate_values': report['exact_direct_sixteen_point_coordinate_values'],
        'tamper_cases_rejected': len(report['tamper_cases_rejected']),
        'exact_results': 'verification_report.json',
        'scope': report['mathematical_scope'],
    }, indent=2))
