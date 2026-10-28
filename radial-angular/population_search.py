"""Exact population search guided by the optimized angular log average.

The dynamic program maximizes an explicit angular-average surrogate over
ALL ordered even ring populations, plus a central point when n is odd.
For its selected population vector, the averaged-log optimum radii have
the closed form b_i=2*S_i+m_i-1, |z|^2=b_i/(n-1). The final candidate
optimizes the entire +/- phase family at these radii by the existing graph
solver. Neither optimization identifies the unrestricted planar maximum.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod
from pathlib import Path
import sys

LEGACY = Path(__file__).resolve().parents[1] / 'configuration-bounds'
if str(LEGACY) not in sys.path:
    sys.path.insert(0, str(LEGACY))
from phase_graph import optimize_signs


def optimized_radii(populations, center=0):
    """Return the exact optimum radii for the ANGULAR LOG AVERAGE only."""
    populations = tuple(populations)
    if type(center) is not int or center not in (0, 1):
        raise ValueError('center must be integer zero or one')
    if not populations or any(type(m) is not int or m < 2 for m in populations):
        raise ValueError('nonempty integer ring populations >=2 required')
    inside = center
    rings = []
    for m in populations:
        rings.append((m, 2 * inside + m - 1))
        inside += m
    n = inside
    if sum(m*b for m, b in rings) != n*(n-1):
        raise ArithmeticError('angular-average moment identity failed')
    return tuple(rings)


def proxy_value(n, rings):
    """Maximum exponential angular log average for given populations."""
    if type(n) is not int or n < 2:
        raise ValueError('integer n>=2 required')
    if any(type(m) is not int or m < 2 or type(b) is not int or b < 1 or (m*b) % 2
           for m, b in rings):
        raise ValueError('integer populations, positive integer b and even m*b required')
    if sum(m*b for m, b in rings) != n*(n-1):
        raise ValueError('prescribed moment identity required')
    return F(prod(m**m * b**(m*b//2) for m, b in rings), (n-1)**(n*(n-1)//2))


def population_dp(n, max_transitions=1000000):
    """Exact DP optimum over all ordered even compositions of n-h.

    h=n mod 2 central points. A state is total population s already placed.
    Appending a ring m gives A=2*s+m-1 and integer multiplicative reward
    m**m*A**(m*A//2). All final normalizing denominators are identical.
    The budget is checked BEFORE any partial enumeration is returned.
    """
    if type(n) is not int or n < 2:
        raise ValueError('integer n>=2 required')
    h = n % 2
    q = (n-h)//2
    planned = q*(q+1)//2
    if max_transitions is not None and (type(max_transitions) is not int or max_transitions < 1):
        raise ValueError('positive integer transition budget or None required')
    if max_transitions is not None and planned > max_transitions:
        raise ValueError('complete population DP needs %s transitions, above budget %s' % (planned, max_transitions))
    best = {h: 1}
    paths = {h: ()}
    transitions = 0
    for s in range(h, n+1, 2):
        for m in range(2, n-s+1, 2):
            A = 2*s+m-1
            candidate = best[s]*m**m*A**(m*A//2)
            transitions += 1
            if s+m not in best or candidate > best[s+m]:
                best[s+m] = candidate
                paths[s+m] = paths[s]+(m,)
    if transitions != planned:
        raise ArithmeticError('complete DP transition count mismatch')
    rings = optimized_radii(paths[n], h)
    N = n*(n-1)//2
    surrogate = F(best[n], (n-1)**N)
    if proxy_value(n, rings) != surrogate:
        raise ArithmeticError('DP witness does not reproduce its surrogate')
    return {'n': n, 'center': h, 'populations': paths[n], 'rings': rings,
            'surrogate': surrogate, 'integer_dp_value': best[n],
            'transitions': transitions,
            'scope': 'Exact optimum of angular-log-average surrogate over all ordered even compositions, not of actual distance product.'}


def ring_value(n, rings, signs, center=0):
    """Exact normalized product for even regular rings and optional origin."""
    rings = tuple((m, F(b)) for m, b in rings)
    signs = tuple(signs)
    if type(n) is not int or n < 2 or type(center) is not int or center not in (0, 1):
        raise ValueError('integer n>=2 and center zero or one required')
    if any(type(m) is not int or m < 2 or m % 2 or b <= 0 for m, b in rings):
        raise ValueError('positive squared radii and even ring populations >=2 required')
    if sum(m for m, b in rings)+center != n:
        raise ValueError('total ring plus center population must equal n')
    if len(signs) != len(rings) or any(type(s) is not int or s not in (-1, 1) for s in signs):
        raise ValueError('one integer +/-1 sign per ring required')
    if any(a >= b for (_, a), (_, b) in zip(rings, rings[1:])):
        raise ValueError('strictly increasing squared ring radii required')
    N = n*(n-1)//2
    Q = sum((m*b for m, b in rings), F(0))
    D = prod(m**m*b**(m*(m-1)//2+center*m) for m, b in rings)
    for i, j in combinations(range(len(rings)), 2):
        mi, a = rings[i]
        mj, b = rings[j]
        g = gcd(mi, mj)
        L = mi*mj//g
        cross = signs[j]**(mi//g)*b**(L//2)-signs[i]**(mj//g)*a**(L//2)
        if cross == 0:
            raise ValueError('collision')
        D *= cross**(2*g)
    return (F(n)/Q)**N*D


def continuous_phase_upper(n, rings, center=0):
    """Rational upper bound for ALL continuous rotations at these radii.

    Each cross resultant is bounded by the sum of the two moduli. The
    pairwise maximizing phases need not be simultaneously attainable.
    This is not an upper bound for other radii or other populations.
    """
    rings = tuple((m, F(b)) for m, b in rings)
    # Reuse all geometric/normalization validation. This value is not used
    # in the upper-bound formula and does not assert an optimal phase.
    ring_value(n, rings, (1,)*len(rings), center)
    N = n*(n-1)//2
    Q = sum((m*b for m, b in rings), F(0))
    D = prod(m**m*b**(m*(m-1)//2+center*m) for m, b in rings)
    for i, j in combinations(range(len(rings)), 2):
        mi, a = rings[i]
        mj, b = rings[j]
        g = gcd(mi, mj)
        L = mi*mj//g
        D *= (b**(L//2)+a**(L//2))**(2*g)
    return (F(n)/Q)**N*D


def population_search(n, max_transitions=1000000, max_phase_states=1000000):
    """Produce an exact candidate plus a separately proved existence bound.

    The surrogate is an attainable lower bound by continuous angular
    averaging. The actual +/- candidate is an explicit configuration.
    Their relative order is NOT assumed by the implementation.
    """
    result = population_dp(n, max_transitions)
    phase = optimize_signs(result['rings'], max_phase_states)
    actual = ring_value(n, result['rings'], phase['signs'], result['center'])
    phase_upper = continuous_phase_upper(n, result['rings'], result['center'])
    if phase_upper < actual or phase_upper < result['surrogate']:
        raise ArithmeticError('continuous-phase upper bound consistency failed')
    return {**result, 'signs': phase['signs'], 'phase_states': phase['states'],
            'value': actual, 'lower_bound': max(actual, result['surrogate']),
            'continuous_phase_upper': phase_upper,
            'scope': 'DP selects the best angular-average surrogate; only signs are then globally optimized at those selected radii. No complete actual population/radius/continuous-phase optimum is claimed.'}


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('n', type=int)
    args = p.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    print(json.dumps(population_search(args.n), default=str, ensure_ascii=False, indent=2))
