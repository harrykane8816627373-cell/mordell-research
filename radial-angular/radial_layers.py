"""Exact radial-gap certificates with variable disk smearing.

See RADIAL_LAYERS.md for the all-n proof. Squared-radius bands use the
Coulomb normalization sum |x_i|^2=n-1; gap endpoints are radii, not squares.
"""
from fractions import Fraction as F
from math import isqrt


def require_n(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n >= 2 required')


def sqrt_upper(q, scale=10**12):
    q = F(q)
    if q < 0 or type(scale) is not int or scale < 1:
        raise ValueError('nonnegative rational and positive integer scale required')
    a, b = q.numerator*scale*scale, q.denominator
    v = isqrt(a//b)
    if v*v*b < a:
        v += 1
    return F(v, scale)


def log_bounds(q, terms=24):
    """Outward rational natural logarithm, with exact power-of-two reduction."""
    q = F(q)
    if q <= 0 or type(terms) is not int or terms < 1:
        raise ValueError('positive rational and positive series length required')
    if q < 1:
        lo, hi = log_bounds(1/q, terms)
        return -hi, -lo
    exponent = 0
    while q >= 2:
        q /= 2
        exponent += 1
    def core(v):
        y = (v-1)/(v+1)
        lo = 2*sum((y**(2*j+1)/(2*j+1) for j in range(terms)), F(0))
        return lo, lo + 2*y**(2*terms+1)/((2*terms+1)*(1-y*y))
    lo, hi = core(q)
    a, b = core(F(2))
    return lo + exponent*a, hi + exponent*b


def smearing_penalty_upper(h, terms=24):
    h = F(h)
    if h <= 0:
        raise ValueError('positive rational smearing multiplier required')
    lo, hi = log_bounds(h, terms)
    return max(F(0), (h*h-1)/2-lo)


def radial_cover_delta_lower(n, layers, radial_width_upper=F(0), h=F(1), scale=10**12, terms=24):
    """Lower bound when center radii fit K bands of total radial width W.

    Holds for every normalized collision-free configuration. No global
    optimality, radius cap, or regular angular distribution is assumed.
    W is a width in RADIUS, not squared radius; W=0 gives exact layers.
    """
    require_n(n)
    if type(layers) is not int or not 1 <= layers <= n:
        raise ValueError('integer layer count in [1,n] required')
    h = F(h)
    width = F(radial_width_upper)
    if h <= 0 or width < 0:
        raise ValueError('positive smearing multiplier and nonnegative radial width required')
    length_fraction = max(F(0), 1-width*sqrt_upper(F(1,2),scale)
                          -2*h*layers*sqrt_upper(F(1, n), scale))
    radial = F(9*n*n, 64*(layers+1)**2)*length_fraction**4
    penalty = n*smearing_penalty_upper(h, terms)
    return {'n': n, 'layers_upper': layers, 'smearing_multiplier': h,
            'radial_band_width_upper': width,
            'remaining_length_fraction_lower': length_fraction,
            'radial_delta_lower': radial, 'smearing_penalty_upper': penalty,
            'delta_lower': max(F(0), radial-penalty),
            'scope': 'All normalized configurations whose center radii fit the stated band count and total radial width.'}


def layer_delta_lower(n, layers, h=F(1), scale=10**12, terms=24):
    """Exact-radius case of radial_cover_delta_lower; no radius cap needed."""
    return radial_cover_delta_lower(n,layers,F(0),h,scale,terms)


def minimum_layers(n, proven_gap_upper, h_values=(F(1), F(99,100)), scale=10**12):
    """Necessary layer count under a separately PROVED upper bound on delta.

    For a global optimizer the gap upper bound must come from an attainable
    lower bound on M_n. Monotonicity permits an exact binary search in K.
    """
    require_n(n)
    gap = F(proven_gap_upper)
    hs = tuple(F(h) for h in h_values)
    if gap < 0 or not hs or any(h <= 0 for h in hs):
        raise ValueError('nonnegative proven gap and positive smearing choices required')
    def certificate(k):
        return max((layer_delta_lower(n, k, h, scale) for h in hs),
                   key=lambda row: row['delta_lower'])
    lo, hi = 1, n
    while lo < hi:
        mid = (lo+hi)//2
        if certificate(mid)['delta_lower'] > gap:
            lo = mid+1
        else:
            hi = mid
    return {'n': n, 'minimum_layers': lo, 'proven_gap_upper': gap,
            'last_excluded_certificate': certificate(lo-1) if lo > 1 else None,
            'scope': 'Necessary condition only; no existence or optimality certification.'}


def free_gap_polynomial_lower(a, b):
    a, b = F(a), F(b)
    if not 0 <= a < b or b*b > 2:
        raise ValueError('radius endpoints 0 <= a < b <= sqrt(2) required')
    return (a+b)*(b-a)**3/24


def free_gap_energy_bounds(a, b, terms=24):
    """Exact-log interval for min_c integral_a^b (c-t^2/2)^2 dt/t."""
    a, b = F(a), F(b)
    fallback = free_gap_polynomial_lower(a, b)
    if a == 0:
        return b**4/16, b**4/16
    lo, hi = log_bounds(b/a, terms)
    constant, numerator = (b**4-a**4)/16, (b*b-a*a)**2/16
    return max(fallback, constant-numerator/lo), constant-numerator/hi


def population_gap_energy_bounds(n, population, a, b, terms=24):
    """Energy of an empty gap with exactly population/n mass to its left."""
    require_n(n)
    a, b = F(a), F(b)
    if type(population) is not int or not 0 <= population <= n:
        raise ValueError('integer cumulative population in [0,n] required')
    fallback = free_gap_polynomial_lower(a, b)
    c = F(population, n)
    if a == 0:
        if population:
            raise ValueError('nonzero constant CDF at radius zero has infinite energy')
        return b**4/16, b**4/16
    lo, hi = log_bounds(b/a, terms)
    constant = -c*(b*b-a*a)/2+(b**4-a**4)/16
    return max(fallback, c*c*lo+constant), c*c*hi+constant


def upper_support_le(q, r2, bound2):
    c = bound2-q-r2
    return c >= 0 and c*c >= 4*q*r2


def lower_support_ge(q, r2, bound2):
    if bound2 == 0:
        return True
    c = q-r2-bound2
    return c >= 0 and c*c >= 4*r2*bound2


def radial_gap_certificate(n, bands, gaps, h=F(1), terms=24):
    """Assigned squared-radius bands plus disjoint PROVED empty radial gaps.

    Each band is (qlo,qhi,population). Every band's smeared support must be
    provably wholly before or after each gap; unknown overlaps are rejected.
    The known cumulative populations give exact weighted gap integrals.
    """
    require_n(n)
    bands = [(F(a), F(b), m) for a, b, m in bands]
    if not bands or any(not 0 <= a <= b or type(m) is not int or m < 1
                        for a, b, m in bands):
        raise ValueError('nonnegative assigned bands and positive populations required')
    if sum(m for a, b, m in bands) != n:
        raise ValueError('population sum must equal n')
    if not sum(m*a for a, b, m in bands) <= n-1 <= sum(m*b for a, b, m in bands):
        raise ValueError('bands incompatible with Coulomb quadratic moment n-1')
    h = F(h)
    if h <= 0:
        raise ValueError('positive smearing multiplier required')
    gaps = [(F(a), F(b)) for a, b in gaps]
    if not gaps:
        raise ValueError('at least one nonempty gap required')
    for a, b in gaps:
        free_gap_polynomial_lower(a, b)
    if any(b > c for (a, b), (c, d) in zip(gaps, gaps[1:])):
        raise ValueError('ordered nonoverlapping gaps required')
    r2, total = 2*h*h/n, F(0)
    rows = []
    for a, b in gaps:
        cumulative = 0
        for qlo, qhi, m in bands:
            if upper_support_le(qhi, r2, a*a):
                cumulative += m
            elif not lower_support_ge(qlo, r2, b*b):
                raise ValueError('a band may intersect the claimed empty gap')
        lo, hi = population_gap_energy_bounds(n, cumulative, a, b, terms)
        total += lo
        rows.append({'a': a, 'b': b, 'cumulative_population': cumulative,
                     'energy_lower': lo, 'energy_upper': hi})
    penalty = n*smearing_penalty_upper(h, terms)
    return {'n': n, 'smearing_multiplier': h, 'radial_energy_lower': total,
            'radial_delta_lower': n*n*total, 'smearing_penalty_upper': penalty,
            'delta_lower': max(F(0), n*n*total-penalty), 'gaps': rows,
            'scope': 'Conditional on the assigned bands; exact cumulative populations are used.'}


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('n', type=int)
    p.add_argument('layers', type=int)
    p.add_argument('--h', default='1')
    args = p.parse_args()
    print(json.dumps(layer_delta_lower(args.n, args.layers, F(args.h)),
                     default=str, indent=2))
