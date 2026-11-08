"""Nonlinear angular deficit certificates; exact acceptance, no trig calls."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, gcd, lcm

PI_LO = F(3141592653589793, 10**15)
PI_HI = F(3141592653589794, 10**15)
LOG_ORDER = 3


def need(ok, message):
    if not ok:
        raise ValueError(message)


def rational(x):
    need(type(x) in (str, int, F), 'exact rational input required')
    try:
        q = F(x)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('malformed rational') from error
    need(max(q.numerator.bit_length(), q.denominator.bit_length()) <= 16000, 'rational arithmetic budget')
    return q


def atan_bounds(q, terms=20):
    lo = sum(((-1)**j*q**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    next_value = (-1)**terms*q**(2*terms+1)/F(2*terms+1)
    return min(lo, lo+next_value), max(lo, lo+next_value)


def check_pi():
    a, b = atan_bounds(F(1, 5))
    c, d = atan_bounds(F(1, 239))
    t2 = 2*F(1, 5)/(1-F(1, 5)**2)
    t4 = 2*t2/(1-t2*t2)
    need((t4-F(1, 239))/(1+t4/F(239)) == 1, 'Machin tangent addition identity')
    need(0 < 4*a-d < 4*b-c < 1, 'Machin angle range excludes tangent periods')
    lo, hi = 16*a-4*d, 16*b-4*c
    need(PI_LO < lo < hi < PI_HI, 'Machin alternating-series pi bracket')
    return lo, hi


def poly_add(a, b):
    out = [F(0)]*max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    return out


def poly_mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def sine_minorant():
    p = [F(0)]*9
    p[2] = PI_LO**2/4
    p[4] = -PI_HI**4/48
    p[6] = PI_LO**6/1440
    p[8] = -PI_HI**8/80640
    # Ignoring the positive x^6 term gives this elementary positivity bound.
    need(p[2]+p[4]+p[8] > 0, 'sine polynomial nonnegative on unit interval')
    return p


def edge_minorant(g, eta):
    need(type(g) is int and g > 0 and 0 < eta < 1, 'positive edge weight and eta in (0,1)')
    y = [eta*v for v in sine_minorant()]
    power = [F(1)]
    out = [F(0)]
    for k in range(1, LOG_ORDER+1):
        power = poly_mul(power, y)
        out = poly_add(out, [g*v/F(k) for v in power])
    return out


def bernstein_integer_coefficients(poly):
    degree = len(poly)-1
    values = [sum((poly[j]*F(comb(k, j), comb(degree, j)) for j in range(k+1)), F(0))
              for k in range(degree+1)]
    D = lcm(*(v.denominator for v in values))
    return [v.numerator*(D//v.denominator) for v in values], D


def subdivide(coeff):
    """de Casteljau at 1/2, returning integer coefficients up to positive scale."""
    degree = len(coeff)-1
    levels = [coeff]
    for _ in range(degree):
        row = levels[-1]
        levels.append([a+b for a, b in zip(row, row[1:])])
    left = [levels[k][0] << (degree-k) for k in range(degree+1)]
    right = [levels[degree-k][k] << k for k in range(degree+1)]
    return left, right


def positive_cover(poly, max_depth=22):
    root, D = bernstein_integer_coefficients(poly)
    leaves = []
    def visit(row, path):
        if min(row) >= 0:
            leaves.append(path)
            return
        need(max(row) >= 0, 'negative polynomial subinterval')
        need(len(path) < max_depth, 'Bernstein depth budget exhausted')
        need(len(leaves) < 20000, 'Bernstein leaf budget exhausted')
        left, right = subdivide(row)
        visit(left, path+'0')
        visit(right, path+'1')
    visit(root, '')
    return leaves


def check_cover(poly, paths):
    need(type(paths) is list and 0 < len(paths) <= 20000, 'nonempty bounded cover list')
    need(all(type(p) is str and len(p) <= 22 and all(q in '01' for q in p) for p in paths), 'binary dyadic paths required')
    need(len(paths) == len(set(paths)), 'duplicate Bernstein leaf')
    leaves = set(paths)
    prefixes = {p[:k] for p in leaves for k in range(len(p)+1)}
    root, D = bernstein_integer_coefficients(poly)
    visited = set()
    def visit(row, path):
        need(path in prefixes, 'incomplete continuum cover')
        if path in leaves:
            need(min(row) >= 0, 'uncertified polynomial leaf')
            visited.add(path)
            return
        left, right = subdivide(row)
        visit(left, path+'0')
        visit(right, path+'1')
    visit(root, '')
    need(visited == leaves, 'hidden descendant or unused Bernstein leaf')
    return len(visited)


def validate_geometry(populations, intervals, center):
    ms = tuple(populations)
    need(3 <= len(ms) <= 6 and all(type(m) is int and 2 <= m <= 30 for m in ms), 'three through six integer ring populations 2..30')
    need(sum(ms) <= 100 and lcm(*ms) <= 240, 'geometry resource budget')
    need(type(center) is int and center in (0, 1), 'center must be integer 0 or 1')
    need(len(intervals) == len(ms) and all(len(v) == 2 for v in intervals), 'one two-endpoint radius interval per ring')
    rs = tuple(tuple(map(rational, v)) for v in intervals)
    need(all(0 < lo <= hi for lo, hi in rs), 'positive ordered endpoints')
    need(all(a[1] < b[0] for a, b in zip(rs, rs[1:])), 'strictly disjoint increasing radius intervals')
    return ms, rs


def cycle_geometry(ms, rs, vertices):
    need(type(vertices) in (list, tuple) and 3 <= len(vertices) <= len(ms), 'simple cycle required')
    vertices = tuple(vertices)
    need(all(type(i) is int and 0 <= i < len(ms) for i in vertices) and len(set(vertices)) == len(vertices), 'distinct valid cycle vertices')
    L = lcm(*(ms[i] for i in vertices))
    out = []
    for i, j in zip(vertices, vertices[1:]+vertices[:1]):
        i, j = sorted((i, j))
        lij = lcm(ms[i], ms[j])
        q = (rs[i][0]/rs[j][1])**lij
        out.append({'edge': (i, j), 'a': L//lij, 'g': gcd(ms[i], ms[j]),
                    'eta': 4*q/(1+q)**2})
    need(sum(v['a'] for v in out) % 2 == 1, 'cycle integer multiplier sum must be odd')
    return out


def old_gamma(edges):
    return 1/sum((F(e['a']**2, e['g'])/e['eta'] for e in edges), F(0))


def jensen_reinforcement(edges, order=6):
    gamma = old_gamma(edges)
    G = sum(e['g'] for e in edges)
    need(0 <= gamma < G, 'log Jensen argument in range')
    return G*sum(((gamma/G)**k/F(k) for k in range(1, order+1)), F(0))


def pair_box_upper(ms, rs, center):
    n = sum(ms)+center
    Q = sum((m*lo*lo for m, (lo, hi) in zip(ms, rs)), F(0))
    value = (F(n)/Q)**(n*(n-1)//2)
    for m, (lo, hi) in zip(ms, rs):
        value *= m**m*hi**(m*(m-1)+2*center*m)
    for i, j in combinations(range(len(ms)), 2):
        L = lcm(ms[i], ms[j])
        value *= (rs[i][1]**L+rs[j][1]**L)**(2*gcd(ms[i], ms[j]))
    return value


def support_polynomial(g, eta, a, slope, intercept):
    poly = edge_minorant(g, eta)
    poly[0] += intercept
    poly[1] -= slope*a
    return poly


def check_certificate(cert):
    fields = {'format', 'date', 'populations', 'radius_intervals', 'center', 'cycles',
              'old_packing_gap', 'jensen_gap', 'nonlinear_gap', 'pair_box_upper',
              'all_phase_upper', 'scope'}
    need(type(cert) is dict and set(cert) == fields, 'certificate schema')
    need(cert['format'] == 'mordell-nonlinear-phase-v1' and cert['date'] == '2026-11-08', 'fixed certificate format/date')
    ms, rs = validate_geometry(cert['populations'], cert['radius_intervals'], cert['center'])
    check_pi()
    need(type(cert['cycles']) is list and 0 < len(cert['cycles']) <= 20, 'bounded cycle list')
    capacities = {edge: F(0) for edge in combinations(range(len(ms)), 2)}
    old = jensen = nonlinear = F(0)
    leaf_count = 0
    for item in cert['cycles']:
        need(type(item) is dict and set(item) == {'vertices', 'weight', 'slope', 'supports'}, 'cycle certificate schema')
        edges = cycle_geometry(ms, rs, item['vertices'])
        weight, slope = rational(item['weight']), rational(item['slope'])
        need(weight > 0 and slope >= 0, 'positive packing weight and nonnegative dual slope')
        need(type(item['supports']) is list and len(item['supports']) == len(edges), 'one support certificate per cycle edge')
        total_intercepts = F(0)
        for edge, support in zip(edges, item['supports']):
            need(type(support) is dict and set(support) == {'edge', 'intercept', 'cover'}, 'support certificate schema')
            need(support['edge'] == list(edge['edge']), 'support edge binding')
            b = rational(support['intercept'])
            need(b >= 0, 'nonnegative support intercept')
            poly = support_polynomial(edge['g'], edge['eta'], edge['a'], slope, b)
            leaf_count += check_cover(poly, support['cover'])
            total_intercepts += b
            capacities[edge['edge']] += weight
        deficit = slope-total_intercepts
        need(deficit > old_gamma(edges), 'claimed nonlinear cycle must improve the old quadratic gap')
        old += weight*old_gamma(edges)
        jensen += weight*jensen_reinforcement(edges)
        nonlinear += weight*deficit
    need(all(value <= 1 for value in capacities.values()), 'packing edge deficit spent twice')
    need(rational(cert['old_packing_gap']) == old and rational(cert['jensen_gap']) == jensen, 'baseline gaps recomputed')
    need(rational(cert['nonlinear_gap']) == nonlinear, 'nonlinear objective recomputed')
    U = pair_box_upper(ms, rs, cert['center'])
    need(rational(cert['pair_box_upper']) == U and rational(cert['all_phase_upper']) == U/(1+nonlinear), 'final box upper binding')
    need(cert['scope'] == 'specified regular-ring populations and radius intervals, all continuous phases; not unrestricted M_n', 'scope escalation')
    return {'old_gap': old, 'jensen_gap': jensen, 'nonlinear_gap': nonlinear,
            'bernstein_leaves': leaf_count, 'upper_ratio_to_old': (1+old)/(1+nonlinear)}
