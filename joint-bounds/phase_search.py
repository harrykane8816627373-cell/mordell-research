"""A uniform radial improvement of the original all-n regular-ring family.

All acceptance decisions use integers/Fraction. For n >= 13, contracting
the squared radius of the innermost six-ring by 7/8 increases every old
continuous-phase configuration by more than 11/10 after recentering and
moment normalization. The implementation additionally optimizes all +/-
phases at the new prescribed radii. It does not optimize all radii/phases
or identify any unrestricted planar maximum.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod
from pathlib import Path
import sys

LEGACY = Path(__file__).resolve().parents[1] / 'configuration-bounds'
if str(LEGACY) not in sys.path:
    sys.path.insert(0, str(LEGACY))
from phase_graph import polar_data, optimize_signs, optimized_polar


def normalized_ring_value(n, rings, signs):
    """Exact value for squared radii b_i, with at most one singleton.

    Center the configuration, then scale its total squared norm to n.
    Ring populations >=2 have zero centroid. The possible singleton
    contributes b/n to the removed centroid moment, independently of phase.
    Cross least common multiples must be even for this rational formula.
    """
    rings = tuple((m, F(b)) for m, b in rings)
    signs = tuple(signs)
    if type(n) is not int or n < 2 or sum(m for m, b in rings) != n:
        raise ValueError('integer n >= 2 must equal the total population')
    if any(type(m) is not int or m < 1 or b <= 0 for m, b in rings):
        raise ValueError('positive integer populations and rational radii required')
    if len(signs) != len(rings) or any(type(s) is not int or s not in (-1, 1) for s in signs):
        raise ValueError('one integer +/- sign per ring required')
    if sum(m == 1 for m, b in rings) > 1:
        raise ValueError('at most one singleton ring is supported')
    if any(a >= b for (_, a), (_, b) in zip(rings, rings[1:])):
        raise ValueError('strictly increasing positive squared radii required')
    N = n * (n - 1) // 2
    Q = sum((m * b for m, b in rings), F(0))
    Q -= sum((b for m, b in rings if m == 1), F(0)) / n
    D = prod(m**m * b**(m * (m - 1) // 2) for m, b in rings)
    for i, j in combinations(range(len(rings)), 2):
        mi, a = rings[i]
        mj, b = rings[j]
        g = gcd(mi, mj)
        L = mi * mj // g
        if L % 2:
            raise ValueError('even cross least common multiples required')
        C = signs[j]**(mi // g) * b**(L // 2) - signs[i]**(mj // g) * a**(L // 2)
        if C == 0:
            raise ValueError('collision')
        D *= C**(2 * g)
    return (F(n) / Q)**N * D


def worst_phase_gain(n):
    """Rational lower bound for ALL continuous ring phases, n >= 13.

    This takes the worst interaction for each pair separately; these
    worst phases need not be simultaneously attainable with many rings.
    Thus this value is a certified lower bound, not a phase minimum.
    """
    if type(n) is not int or n < 13:
        raise ValueError('integer n >= 13 required')
    rings = polar_data(n)[2]
    if rings[0] != (6, 6):
        raise ArithmeticError('unexpected original inner-ring data')
    N = n * (n - 1) // 2
    t = F(7, 8)
    Q = F(n * n) - sum((F(b, n) for m, b in rings if m == 1), F(0))
    gain = t**15 * (Q / (Q - F(9, 2)))**N
    for m, b in rings[1:]:
        g = gcd(6, m)
        p = 3 * m // g
        gain *= ((F(b)**p + (6 * t)**p) / (F(b)**p + 6**p))**(2 * g)
    return gain


def analytic_log_gain_floor(n):
    """All-n proof's rational logarithmic lower bound, n >= 24.

    No logarithm is evaluated. RADIAL_CONTRACTION.md proves this is a
    lower bound for the actual logarithmic gain for every continuous phase.
    """
    if type(n) is not int or n < 24:
        raise ValueError('integer n >= 24 required')
    return F(9, 4) * F(n - 1, n) - F(225, 112) - F(1, 336) - F(12, 343)


def contracted_polar(n, max_states=1000000):
    """Exact new +/- maximum at the prescribed contracted squared radii."""
    if type(n) is not int or n < 13:
        raise ValueError('integer n >= 13 required')
    original = optimized_polar(n, max_states)
    rings = list((m, F(b)) for m, b in original['rings'])
    rings[0] = (6, F(21, 4))
    # Multiply all squared radii by four for the existing integer graph
    # solver. Global scaling preserves every phase comparison.
    graph_rings = tuple((m, int(4 * b)) for m, b in rings)
    graph = optimize_signs(graph_rings, max_states)
    candidate = normalized_ring_value(n, rings, graph['signs'])
    old_sign_value = normalized_ring_value(n, rings, original['signs'])
    if candidate < old_sign_value or old_sign_value / original['value'] < worst_phase_gain(n):
        raise ArithmeticError('radial gain or exact phase comparison failed')
    if candidate <= F(11, 10) * original['value']:
        raise ArithmeticError('uniform rational gain failed')
    return {'n': n, 'rings': tuple(rings), 'signs': graph['signs'],
            'value': candidate, 'old_phase_value': original['value'],
            'old_sign_value_at_new_radii': old_sign_value,
            'phase_states': graph['states'],
            'scope': 'Exact +/- phase maximum at prescribed contracted radii only.'}


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('n', type=int)
    args = p.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    print(json.dumps(contracted_polar(args.n), default=str, ensure_ascii=False, indent=2))
