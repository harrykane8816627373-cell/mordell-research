"""Small-denominator logarithm enclosures and integer-point certificates."""
from fractions import Fraction as F
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'joint-bounds'))
from subset_moment import log_interval, outward_decimal


def integer_value(points):
    """Center and normalize implicitly, with one final rational reduction."""
    points = tuple(tuple(p) for p in points)
    if len(points) < 2 or any(len(p) != 2 or any(type(x) is not int for x in p) for p in points):
        raise ValueError('at least two integer coordinate pairs required')
    if len(set(points)) != len(points):
        raise ValueError('distinct points required')
    n = len(points)
    N = n*(n-1)//2
    sx, sy = (sum(p[k] for p in points) for k in (0, 1))
    D = n*sum(x*x+y*y for x, y in points)-sx*sx-sy*sy
    distances = [(x-u)**2+(y-v)**2 for i, (x, y) in enumerate(points)
                 for u, v in points[i+1:]]
    return F(n**(2*N)*prod(distances), D**N)


def bounded_log(value, bits=64, terms=28):
    """Outward log enclosure without high powers of giant input integers.

    Reduce value=2**e*r, 1<=r<2, bracket r by adjacent dyadics of
    denominator 2**bits, then use monotonicity and rigorous log series.
    The input rational is exact, not a rounded floating-point product.
    """
    value = F(value)
    if value <= 0 or type(bits) is not int or bits < 2:
        raise ValueError('positive rational value and integer bits>=2 required')
    if type(terms) is not int or terms < 1:
        raise ValueError('positive integer terms required')
    if value < 1:
        lo, hi = bounded_log(1/value, bits, terms)
        return -hi, -lo
    e = max(0, value.numerator.bit_length()-value.denominator.bit_length())
    numerator, denominator = value.numerator, value.denominator << e
    if numerator < denominator:
        e -= 1
        denominator >>= 1
    scale = 1 << bits
    cell, remainder = divmod(numerator*scale, denominator)
    lo = log_interval(F(cell, scale), terms)[0]
    hi = log_interval(F(cell+bool(remainder), scale), terms)[1]
    log2 = log_interval(F(2), terms)
    return lo+e*log2[0], hi+e*log2[1]
