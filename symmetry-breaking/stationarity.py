"""Exact regular-ring stationarity obstructions; see STATIONARITY.md.

The general population theorem allows arbitrary points outside one regular
polygon. The implementation evaluates specified even regular-ring families
with rational squared radii and real polynomial phases +/-1.
"""
from fractions import Fraction as F
from math import prod


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def binomial_remainder(poly, m, c):
    """Coefficients modulo X**m-c, in increasing degree order."""
    if type(m) is not int or m < 1 or F(c) == 0:
        raise ValueError('positive integer degree and nonzero constant required')
    out = [F(0)]*m
    for k, coefficient in enumerate(poly):
        out[k % m] += F(coefficient)*F(c)**(k//m)
    return out


def regular_component_possible(m, other_nonzero, center=0):
    """Necessary count condition ONLY, never a sufficient existence test."""
    if (type(m) is not int or m < 2 or type(other_nonzero) is not int or
            other_nonzero < 0 or type(center) is not int or center not in (0, 1)):
        raise ValueError('m>=2, other_nonzero>=0, center in {0,1}, all integers')
    return other_nonzero == 0 or m <= other_nonzero


def ring_stationarity(rings, signs, center=0):
    rings = tuple((m, F(b)) for m, b in rings)
    signs = tuple(signs)
    if type(center) is not int or center not in (0, 1) or not rings:
        raise ValueError('nonempty rings and center in {0,1} required')
    if any(type(m) is not int or m < 2 or m % 2 or b <= 0 for m, b in rings):
        raise ValueError('even populations >=2 and positive rational squared radii required')
    if len(signs) != len(rings) or any(type(s) is not int or s not in (-1, 1) for s in signs):
        raise ValueError('one integer +/-1 sign per ring required')
    if any(a >= b for (_, a), (_, b) in zip(rings, rings[1:])):
        raise ValueError('strictly increasing squared radii required')
    n = center+sum(m for m, b in rings)
    Q = sum((m*b for m, b in rings), F(0))
    kappa = F(n*(n-1), 2)/Q
    rows = []
    residual_lower = F(0)
    for i, (m, b) in enumerate(rings):
        c = signs[i]*b**(m//2)
        external = [F(1)]
        for j, (mj, bj) in enumerate(rings):
            if i == j:
                continue
            factor = [-signs[j]*bj**(mj//2)]+[F(0)]*(mj-1)+[F(1)]
            external = multiply(external, factor)
        C = F(m-1, 2)+center-kappa*b
        numerator = [(k+C)*q for k, q in enumerate(external)]
        remainder = binomial_remainder(numerator, m, c)
        coefficient_norm = sum((a*a*b**k for k, a in enumerate(remainder)), F(0))
        denominator_upper = b*prod((b**(mj//2)+bj**(mj//2))**2
                                  for j, (mj, bj) in enumerate(rings) if j != i)
        lower = m*coefficient_norm/denominator_upper
        residual_lower += lower
        rows.append({'population': m, 'squared_radius': b, 'constant_C': C,
                     'remainder': remainder, 'coefficient_norm': coefficient_norm,
                     'residual_norm_squared_lower': lower,
                     'stationary_on_this_ring': coefficient_norm == 0})
    # sin(pi/m)>=2/m, and (sqrt(b)-sqrt(a))^2 >= (b-a)^2/(2(a+b)).
    separation_candidates = [16*b/F(m*m) for m, b in rings]
    separation_candidates += [(b-a)**2/(2*(a+b)) for i, (_, a) in enumerate(rings)
                              for _, b in rings[i+1:]]
    if center:
        separation_candidates.append(rings[0][1])
    separation_lower = min(separation_candidates)
    gain_lower = residual_lower/(2*n*n/separation_lower+kappa)
    return {'n': n, 'center': center, 'moment': Q, 'multiplier': kappa,
            'stationary': all(r['stationary_on_this_ring'] for r in rows),
            'residual_norm_squared_lower': residual_lower,
            'scale_invariant_residual_lower': Q*residual_lower,
            'pair_separation_squared_lower': separation_lower,
            'log_global_gap_lower': gain_lower, 'rings': rows}
