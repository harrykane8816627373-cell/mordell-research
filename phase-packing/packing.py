"""Exact fractional cycle packing for fixed regular-ring phase families.

The certificate checker reconstructs every claimed LP coefficient and
checks primal/dual feasibility and equality using Fraction only.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, lcm, prod


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    require(type(value) in (str, int, F), 'exact integer, rational string or Fraction required')
    try:
        return F(value)
    except (ZeroDivisionError, ValueError) as error:
        raise ValueError('malformed rational coefficient') from error


def ring_box(populations, intervals):
    ms = tuple(populations)
    require(1 <= len(ms) <= 8, 'one through eight rings supported')
    require(all(type(m) is int and 2 <= m <= 100 for m in ms), 'integer populations 2..100 required')
    require(sum(ms) <= 100 and lcm(*ms) <= 240, 'supported arithmetic budget: total population <=100, lcm <=240')
    require(len(intervals) == len(ms), 'one radius interval per ring')
    rs = tuple(tuple(map(rational, pair)) for pair in intervals)
    require(all(len(pair) == 2 for pair in rs), 'two interval endpoints required')
    require(all(0 < lo <= hi for lo, hi in rs), 'positive ordered endpoints required')
    require(all(max(q.numerator.bit_length(), q.denominator.bit_length()) <= 128
                for pair in rs for q in pair), 'radius endpoint bit budget exceeded')
    require(all(a[1] < b[0] for a, b in zip(rs, rs[1:])), 'strictly separated increasing radius intervals required')
    return ms, rs


def all_cycles(K, maximum_length):
    require(type(maximum_length) is int and 3 <= maximum_length <= 8,
            'maximum cycle length must be integer 3..8')
    out = []
    for size in range(3, min(K, maximum_length)+1):
        for vertices in combinations(range(K), size):
            first = vertices[0]
            for tail in permutations(vertices[1:]):
                if tail[0] < tail[-1]:
                    out.append((first,) + tail)
                    require(len(out) <= 5000, 'cycle enumeration exceeds supported budget')
    return out


def cycle_data(ms, rs, vertices):
    L = lcm(*(ms[i] for i in vertices))
    edges = [tuple(sorted((i, j))) for i, j in zip(vertices, vertices[1:]+vertices[:1])]
    multipliers = [L//lcm(ms[i], ms[j]) for i, j in edges]
    if sum(multipliers) % 2 == 0:
        return None
    inverse_gap = F(0)
    for (i, j), a in zip(edges, multipliers):
        # r_i/r_j < 1; eta is increasing in that ratio.
        lij = lcm(ms[i], ms[j])
        q = (rs[i][0]/rs[j][1])**lij
        eta_lower = 4*q/(1+q)**2
        inverse_gap += F(a*a, gcd(ms[i], ms[j]))/eta_lower
    return {'vertices': vertices, 'edges': tuple(edges),
            'multipliers': tuple(multipliers), 'gamma': 1/inverse_gap}


def model(populations, intervals, maximum_length=3):
    ms, rs = ring_box(populations, intervals)
    edges = tuple(combinations(range(len(ms)), 2))
    candidates = []
    for vertices in all_cycles(len(ms), maximum_length):
        record = cycle_data(ms, rs, vertices)
        if record is not None:
            candidates.append(record)
    require(len(candidates) <= 300, 'LP exceeds 300-cycle supported budget')
    A = [[F(edge in item['edges']) for item in candidates] for edge in edges]
    c = [item['gamma'] for item in candidates]
    return ms, rs, edges, candidates, A, c


def solve_packing(A, c):
    """Bland-rule rational simplex, starting at the all-slack solution.

    The returned certificate is checked independently; optimality never
    depends on trusting the pivot implementation or its stopping flag.
    """
    rows, cols = len(A), len(c)
    require(rows <= 100 and cols <= 300, 'supported LP budget exceeded')
    require(all(len(row) == cols and all(x in (0, 1) for x in row) for row in A), '0/1 incidence matrix required')
    require(all(x > 0 for x in c), 'strictly positive objective coefficients required')
    require(all(any(A[i][j] for i in range(rows)) for j in range(cols)), 'each cycle needs an edge')
    tab = [[F(x) for x in A[i]] + [F(i == j) for j in range(rows)] for i in range(rows)]
    b = [F(1)]*rows
    cost = list(c) + [F(0)]*rows
    basis = [cols+i for i in range(rows)]
    pivots = 0
    while True:
        cb = [cost[k] for k in basis]
        reduced = [cost[j]-sum((cb[i]*tab[i][j] for i in range(rows)), F(0))
                   for j in range(cols+rows)]
        enter = next((j for j, value in enumerate(reduced) if j not in basis and value > 0), None)
        if enter is None:
            break
        eligible = [(b[i]/tab[i][enter], basis[i], i) for i in range(rows) if tab[i][enter] > 0]
        require(eligible, 'unexpected unbounded LP')
        leave = min(eligible)[2]
        p = tab[leave][enter]
        tab[leave] = [x/p for x in tab[leave]]
        b[leave] /= p
        for i in range(rows):
            if i != leave:
                q = tab[i][enter]
                if q:
                    tab[i] = [x-q*y for x, y in zip(tab[i], tab[leave])]
                    b[i] -= q*b[leave]
        basis[leave] = enter
        pivots += 1
        require(pivots <= 10000, 'simplex pivot budget exceeded')
        require(all(x >= 0 for x in b), 'simplex lost primal feasibility')
    primal = [F(0)]*cols
    for i, basic in enumerate(basis):
        if basic < cols:
            primal[basic] = b[i]
    dual = [sum((cost[basis[i]]*tab[i][cols+e] for i in range(rows)), F(0)) for e in range(rows)]
    return primal, dual, pivots


def packing_check(A, c, primal, dual):
    require(len(primal) == len(c) and len(dual) == len(A), 'certificate dimension mismatch')
    require(all(x >= 0 for x in primal+dual), 'negative certificate weight')
    for row in A:
        require(sum((a*x for a, x in zip(row, primal)), F(0)) <= 1, 'primal edge capacity exceeded')
    for j, value in enumerate(c):
        require(sum((A[i][j]*dual[i] for i in range(len(A))), F(0)) >= value, 'dual cycle constraint violated')
    objective = sum((x*y for x, y in zip(c, primal)), F(0))
    require(objective == sum(dual, F(0)), 'primal/dual objective mismatch')
    return objective


def greedy_gap(candidates):
    used = set()
    value = F(0)
    chosen = []
    for item in sorted(candidates, key=lambda row: (-row['gamma'], row['vertices'])):
        if not any(edge in used for edge in item['edges']):
            used.update(item['edges'])
            value += item['gamma']
            chosen.append(item['vertices'])
    return value, chosen


def normalized_pair_upper(ms, rs, center=0):
    require(type(center) is int and center in (0, 1), 'center must be integer 0 or 1')
    n = sum(ms)+center
    Q_lower = sum((m*lo*lo for m, (lo, hi) in zip(ms, rs)), F(0))
    value = (F(n)/Q_lower)**(n*(n-1)//2)
    for m, (lo, hi) in zip(ms, rs):
        value *= m**m*hi**(m*(m-1)+2*center*m)
    for i, j in combinations(range(len(ms)), 2):
        L = lcm(ms[i], ms[j])
        value *= (rs[i][1]**L + rs[j][1]**L)**(2*gcd(ms[i], ms[j]))
    return value


def build_certificate(populations, intervals, maximum_length=3, center=0):
    ms, rs, edges, candidates, A, c = model(populations, intervals, maximum_length)
    primal, dual, pivots = solve_packing(A, c)
    optimum = packing_check(A, c, primal, dual)
    greedy, chosen = greedy_gap(candidates)
    U = normalized_pair_upper(ms, rs, center)
    return {'format': 'mordell-phase-packing-v1', 'populations': list(ms),
            'radius_intervals': [[str(a), str(b)] for a, b in rs],
            'center': center, 'maximum_cycle_length': maximum_length,
            'cycles': [list(item['vertices']) for item in candidates],
            'gamma': [str(x) for x in c], 'primal': [str(x) for x in primal],
            'dual_edges': [list(edge) for edge in edges], 'dual': [str(x) for x in dual],
            'optimal_log_gap': str(optimum), 'greedy_log_gap': str(greedy),
            'greedy_cycles': [list(x) for x in chosen],
            'pair_upper': str(U), 'all_phase_upper': str(U/(1+optimum)),
            'pivots': pivots,
            'scope': 'specified regular-ring populations and radius intervals, all continuous phases; not unrestricted M_n'}


def check_certificate(cert):
    expected = {'format', 'populations', 'radius_intervals', 'center',
                'maximum_cycle_length', 'cycles', 'gamma', 'primal', 'dual_edges',
                'dual', 'optimal_log_gap', 'greedy_log_gap', 'greedy_cycles',
                'pair_upper', 'all_phase_upper', 'pivots', 'scope'}
    require(type(cert) is dict and set(cert) == expected, 'unknown or missing certificate field')
    require(cert['format'] == 'mordell-phase-packing-v1', 'wrong certificate format')
    ms, rs, edges, candidates, A, c = model(cert['populations'], cert['radius_intervals'], cert['maximum_cycle_length'])
    require(cert['cycles'] == [list(row['vertices']) for row in candidates], 'missing, reordered or tampered cycle list')
    require([rational(x) for x in cert['gamma']] == c, 'tampered geometric gamma coefficient')
    require(cert['dual_edges'] == [list(e) for e in edges], 'tampered edge indexing')
    primal, dual = [rational(x) for x in cert['primal']], [rational(x) for x in cert['dual']]
    optimum = packing_check(A, c, primal, dual)
    require(rational(cert['optimal_log_gap']) == optimum, 'tampered optimal objective')
    greedy, chosen = greedy_gap(candidates)
    require(rational(cert['greedy_log_gap']) == greedy and cert['greedy_cycles'] == [list(x) for x in chosen], 'tampered greedy comparison')
    U = normalized_pair_upper(ms, rs, cert['center'])
    require(rational(cert['pair_upper']) == U, 'tampered pair upper')
    require(rational(cert['all_phase_upper']) == U/(1+optimum), 'tampered all-phase upper')
    require(type(cert['pivots']) is int and 0 <= cert['pivots'] <= 10000, 'invalid pivot metadata')
    require(cert['scope'] == 'specified regular-ring populations and radius intervals, all continuous phases; not unrestricted M_n', 'wrong scope')
    return optimum
