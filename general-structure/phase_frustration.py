"""Exact certificates for continuous-phase frustration of regular rings.

See CONTINUOUS_PHASES.md. Bounds apply to the specified populations and
radii, for ALL continuous phases; they are not upper bounds for unrestricted
M_n. Only Python's standard library is used. Checks survive python -O.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm, prod


def v2(m):
    if type(m) is not int or m < 1:
        raise ValueError('positive integer required')
    return (m & -m).bit_length() - 1


def populations_valid(populations):
    ms = tuple(populations)
    if not ms or any(type(m) is not int or m < 2 for m in ms):
        raise ValueError('nonempty populations, each an integer >= 2')
    return ms


def simultaneous_pair_maxima(populations):
    """Complete criterion, for continuous phases, not only sign phases."""
    ms = populations_valid(populations)
    es = [v2(m) for m in ms]
    smallest = min(es)
    return all(es.count(e) <= (2 if e == smallest else 1) for e in set(es))


def maximizing_signs(populations):
    """A simultaneous equality witness when one exists; None otherwise."""
    ms = populations_valid(populations)
    if not simultaneous_pair_maxima(ms):
        return None
    signs = [-1] * len(ms)
    es = [v2(m) for m in ms]
    if es.count(min(es)) == 2:
        signs[es.index(min(es))] = 1
    return tuple(signs)


def bad_triangles(populations):
    """Cycle certificates: all three integer L/L_ij multipliers are odd."""
    ms = populations_valid(populations)
    out = []
    for tri in combinations(range(len(ms)), 3):
        L = lcm(*(ms[i] for i in tri))
        edges = tuple(combinations(tri, 2))
        multipliers = tuple(L // lcm(ms[i], ms[j]) for i, j in edges)
        if all(a % 2 for a in multipliers):
            out.append((tri, edges, multipliers))
    return out


def validate_rings(rings, center):
    if type(center) is not int or center not in (0, 1):
        raise ValueError('center must be integer 0 or 1')
    rings = tuple((m, F(r)) for m, r in rings)
    populations_valid(m for m, r in rings)
    if any(r <= 0 for m, r in rings):
        raise ValueError('positive rational radii required')
    if any(a >= b for (_, a), (_, b) in zip(rings, rings[1:])):
        raise ValueError('radii must be strictly increasing')
    return rings


def pairwise_relaxation(rings, center=0):
    """Normalized product U obtained by independently maximizing each pair."""
    rings = validate_rings(rings, center)
    n = center + sum(m for m, r in rings)
    Q = sum((m*r*r for m, r in rings), F(0))
    N = n*(n-1)//2
    value = F(n, 1)**N / Q**N
    for m, r in rings:
        value *= m**m * r**(m*(m-1) + 2*center*m)
    for (m, r), (k, s) in combinations(rings, 2):
        g = gcd(m, k)
        L = lcm(m, k)
        value *= (r**L + s**L)**(2*g)
    return value


def sign_product(rings, signs, center=0):
    """Exact normalized value for X**m = sign*r**m, for arbitrary m>=2."""
    rings = validate_rings(rings, center)
    signs = tuple(signs)
    if len(signs) != len(rings) or any(type(s) is not int or s not in (-1, 1) for s in signs):
        raise ValueError('one integer +/-1 sign per ring')
    n = center + sum(m for m, r in rings)
    Q = sum((m*r*r for m, r in rings), F(0))
    value = (F(n, 1)/Q)**(n*(n-1)//2)
    for m, r in rings:
        value *= m**m * r**(m*(m-1)+2*center*m)
    for i, j in combinations(range(len(rings)), 2):
        m, r = rings[i]
        k, s = rings[j]
        g = gcd(m, k)
        value *= (signs[i]**(k//g)*r**lcm(m, k)
                  - signs[j]**(m//g)*s**lcm(m, k))**(2*g)
    return value


def triangle_gap(rings, triangle):
    """Positive rational log gap for one incompatible triangle."""
    rings = validate_rings(rings, 0)
    tri = tuple(triangle)
    if len(tri) != 3 or len(set(tri)) != 3 or any(type(i) is not int or i < 0 or i >= len(rings) for i in tri):
        raise ValueError('three distinct valid indices required')
    L = lcm(*(rings[i][0] for i in tri))
    inverse = F(0)
    for i, j in combinations(tri, 2):
        m, r = rings[i]
        k, s = rings[j]
        lij = lcm(m, k)
        a = L//lij
        if a % 2 == 0:
            raise ValueError('triangle has no odd-cycle obstruction')
        A, B = r**lij, s**lij
        eta = 4*A*B/(A+B)**2
        inverse += F(a*a, gcd(m, k))/eta
    return 1/inverse


def family_bound(rings, center=0):
    """All-phase upper bound with a deterministic edge-disjoint packing.

    Each selected triangle spends three distinct edge deficits; reuse of
    an edge is forbidden, avoiding double counting. Greedy selection is
    not claimed to solve a maximum packing problem.
    """
    rings = validate_rings(rings, center)
    ms = tuple(m for m, r in rings)
    candidates = [(triangle_gap(rings, tri), tri, edges, multipliers)
                  for tri, edges, multipliers in bad_triangles(ms)]
    candidates.sort(key=lambda row: (-row[0], row[1]))
    spent = set()
    selected = []
    gap = F(0)
    for gamma, tri, edges, multipliers in candidates:
        if any(edge in spent for edge in edges):
            continue
        spent.update(edges)
        gap += gamma
        selected.append({'triangle': tri, 'edges': edges,
                         'odd_multipliers': multipliers, 'log_gap': gamma})
    upper = pairwise_relaxation(rings, center)
    return {'populations': ms, 'radii': tuple(r for m, r in rings), 'center': center,
            'all_pairs_simultaneously_maximal': simultaneous_pair_maxima(ms),
            'bad_triangle_count': len(candidates), 'selected_triangles': selected,
            'log_gap_lower': gap, 'pairwise_relaxation': upper,
            'all_continuous_phase_upper': upper/(1+gap),
            'scope': 'fixed populations and radii, all continuous phases; not unrestricted M_n'}


def three_equal_ring_value(m, radii):
    """Value at polynomial phases 0,2*pi/3,4*pi/3, NOT a claimed optimum."""
    if type(m) is not int or m < 2 or len(radii) != 3:
        raise ValueError('integer m>=2 and exactly three radii required')
    rings = validate_rings([(m, r) for r in radii], 0)
    n = 3*m
    Q = sum((m*r*r for mm, r in rings), F(0))
    value = (F(n, 1)/Q)**(n*(n-1)//2)
    for mm, r in rings:
        value *= m**m * r**(m*(m-1))
    for (_, r), (_, s) in combinations(rings, 2):
        A, B = r**m, s**m
        value *= (A*A+A*B+B*B)**m
    return value


def examples():
    rings = ((4, F(4, 5)), (4, F(1)), (4, F(6, 5)))
    bound = family_bound(rings)
    sign_values = [(sign_product(rings, s), s) for s in product((-1, 1), repeat=3)]
    best, signs = max(sign_values)
    continuous = three_equal_ring_value(4, [r for m, r in rings])
    return {'square_triple': bound,
            'continuous_vs_sign': {'rings': rings, 'best_signs': signs,
                                   'best_sign_value': best,
                                   'continuous_polynomial_phases': ['0', '2*pi/3', '4*pi/3'],
                                   'continuous_value': continuous,
                                   'exact_ratio': continuous/best,
                                   'ratio_strictly_above_two': continuous > 2*best,
                                   'claim': 'improves this fixed sign family, not the known unrestricted n=12 lower bound'},
            'mixed_ten': family_bound(((2, F(1, 2)), (4, F(1)), (4, F(3, 2)))),
            'compatible_eight': family_bound(((3, F(4, 5)), (3, F(1)), (2, F(6, 5))))}


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


if __name__ == '__main__':
    import json
    print(json.dumps(encode(examples()), ensure_ascii=False, indent=2))
