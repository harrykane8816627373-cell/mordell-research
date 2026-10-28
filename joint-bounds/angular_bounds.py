"""Exact angular moments added to the existing radial Coulomb certificate.

The analytic theorem is ANGULAR_PROOF.md. All returned numbers are rational.
Small moments are necessary, never sufficient, for global optimality.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'configuration-bounds'))
from radial_distribution import band_certificate


def require_n(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n >= 2 required')


def sqrt_upper(q, scale=10**8):
    q = F(q)
    if q < 0 or type(scale) is not int or scale < 1:
        raise ValueError('nonnegative rational and positive integer scale required')
    a, b = q.numerator * scale**2, q.denominator
    k = isqrt(a // b)
    if k*k*b < a:
        k += 1
    return F(k, scale)


def log2_bounds(terms=24):
    if type(terms) is not int or terms < 1:
        raise ValueError('positive integer series length required')
    y = F(1, 3)
    lo = 2 * sum((y**(2*j+1) / (2*j+1) for j in range(terms)), F(0))
    return lo, lo + 2*y**(2*terms+1) / ((2*terms+1)*(1-y*y))


def log_bounds(q, terms=32):
    """Rational logarithm enclosure with binary range reduction."""
    q = F(q)
    if q <= 0 or type(terms) is not int or terms < 1:
        raise ValueError('positive rational and positive integer terms required')
    if q < 1:
        lo, hi = log_bounds(1/q, terms)
        return -hi, -lo
    exponent = q.numerator.bit_length() - q.denominator.bit_length()
    scaled = q / (1 << exponent)
    if scaled < 1:
        exponent -= 1
        scaled *= 2
    y = (scaled-1)/(scaled+1)
    lo = 2*sum((y**(2*j+1)/(2*j+1) for j in range(terms)), F(0))
    hi = lo + 2*y**(2*terms+1)/((2*terms+1)*(1-y*y))
    l2, u2 = log2_bounds(terms)
    return exponent*l2 + lo, exponent*u2 + hi


def gap_upper_from_lower(n, attained_lower, terms=32):
    """Requires a separately proved positive lower value L <= M_n."""
    require_n(n)
    llo, lhi = log_bounds(attained_lower, terms)
    logn_hi = log_bounds(n, terms)[1]
    cross_hi = log_bounds(F(2*n, n-1), terms)[1]
    pairs = n*(n-1)//2
    return F(n, 2)*logn_hi + pairs*cross_hi - F(pairs, 2) - llo


def strict_disk_correction(n, scale=10**8, terms=24):
    """Positive rational c with M_n <= U_n exp(-c), for every n>=2.

    Uses centered GLOBAL maximizers. The exact analytic correction is
    (log 2)^2 / sum_{k=2}^n (1+sqrt(2/(n-1)))^(2*k)/k.
    No assertion of a best known envelope or linear asymptotic gain.
    """
    require_n(n)
    b = (1 + sqrt_upper(F(2, n-1), scale))**2
    denominator = sum((b**k / k for k in range(2, n+1)), F(0))
    loglo, loghi = log2_bounds(terms)
    return {'n': n, 'ratio_square_upper': b, 'log2_lower': loglo,
            'log2_upper': loghi, 'denominator_upper': denominator,
            'delta_lower': loglo**2 / denominator,
            'scope': 'Unconditional correction of disk U_n only; not a sharp finite-n value.'}


def angular_cost(moment_abs2_lower, support_r2_upper):
    """Requires independent bounds |p_k(x)|^2 >= supplied rational values.

    R^2 must bound (max_i |x_i| + sqrt(2/n))^2, the smeared point
    support. The rotational background need not be contained in this disk.
    """
    r2 = F(support_r2_upper)
    if r2 <= 0:
        raise ValueError('positive support squared radius required')
    total = F(0)
    for k, value in moment_abs2_lower.items():
        value = F(value)
        if type(k) is not int or k < 1 or value < 0:
            raise ValueError('positive integer moment indices and nonnegative squared lower bounds required')
        total += value / (k * r2**k)
    return total


def joint_certificate(n, bands, partition, moment_abs2_lower, scale=10**8):
    """Entire-class certificate with BOTH assigned radial bands and moments.

    bands are Coulomb-scale squared-radius bands (lo,hi,population).
    No global coverage or existence of a class is asserted. Radial and angular
    terms are additive because they occupy orthogonal Dirichlet modes.
    """
    require_n(n)
    bands = [(F(a), F(b), m) for a, b, m in bands]
    radial = band_certificate(n, bands, partition)
    qmax = max(b for a, b, m in bands)
    r2 = qmax + F(2, n) + 2*sqrt_upper(qmax*F(2, n), scale)
    angular = angular_cost(moment_abs2_lower, r2)
    return {'n': n, 'support_r2_upper': r2,
            'radial_delta_lower': radial['delta_lower'],
            'angular_delta_lower': angular,
            'delta_lower': radial['delta_lower'] + angular,
            'radial_cells': radial['cells'],
            'scope': 'Conditional on every supplied radial band and every moment lower bound; arbitrary angles allowed.'}


def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def rational_shape_data(points, degree):
    """Center rational points and compute normalized Coulomb squared moments.

    x=sqrt((n-1)/Q)*v. Although x can be irrational, its squared radii and
    |p_k(x)|^2 are rational, so certificates need no complex floating point.
    """
    points = [(F(a), F(b)) for a, b in points]
    n = len(points)
    require_n(n)
    if type(degree) is not int or degree < 1 or len(set(points)) != n:
        raise ValueError('positive degree and distinct points required')
    mean = (sum((a for a, b in points), F(0))/n,
            sum((b for a, b in points), F(0))/n)
    v = [(a-mean[0], b-mean[1]) for a, b in points]
    q = sum((a*a+b*b for a, b in v), F(0))
    factor = F(n-1)/q
    qs = [factor*(a*a+b*b) for a, b in v]
    powers = [(F(1), F(0)) for _ in v]
    moments = {}
    for k in range(1, degree+1):
        powers = [cmul(a, b) for a, b in zip(powers, v)]
        re = sum((a for a, b in powers), F(0))
        im = sum((b for a, b in powers), F(0))
        moments[k] = factor**k*(re*re+im*im)
    return {'n': n, 'quadratic_moment_before_scaling': q,
            'coulomb_radii_squared': qs, 'moment_abs2': moments}


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('n', type=int)
    args = p.parse_args()
    print(json.dumps(strict_disk_correction(args.n), default=str, indent=2))
