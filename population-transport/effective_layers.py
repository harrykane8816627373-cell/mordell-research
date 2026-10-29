"""Effective finite-n corollary of spectral capacity and radial transport."""
from fractions import Fraction as F
from math import isqrt


def sqrt_interval(q,scale=10**12):
    q=F(q)
    if q<0 or type(scale)is not int or scale<1:
        raise ValueError('nonnegative rational and positive integer scale required')
    num,den=q.numerator*scale*scale,q.denominator
    s=isqrt(num//den)
    return F(s,scale),F(s+int(s*s*den<num),scale)


def effective_layers(n,scale=10**12):
    """For GLOBAL optimizers, n>=10000: K>sqrt(n)/2-11*n**(1/4)/20.

    The integer output uses an outward lower endpoint of this strict bound.
    No claim about the exact layer count, or smaller n, is made here.
    """
    if type(n)is not int or n<10000:
        raise ValueError('integer n>=10000 required for this corollary')
    lo,hi=sqrt_interval(n,scale)
    fourth_upper=sqrt_interval(hi,scale)[1]
    strict_lower=lo/2-F(11,20)*fourth_upper
    integer_lower=max(1,strict_lower.numerator//strict_lower.denominator+1)
    return {'n':n,'strict_layer_lower_endpoint':strict_lower,
            'minimum_layers':integer_lower,
            'sqrt_interval':[lo,hi],'fourth_root_upper':fourth_upper,
            'scope':'Necessary count for every global optimizer; the closed corollary requires n>=10000.'}
