"""Rational lower quadrature for the exact unordered smearing pair loss.

The proof of convexity is in PAIR_CONVEXITY.md.  Quadrature uses monotone
factor bounds, never an unproved monotonicity of their product.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'overlap-budget'))
from overlap_cost import certified_pi_upper


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def sqrt_bounds(value,mesh=10**10):
    q=F(value)
    if q<0 or type(mesh)is not int or mesh<1:raise ValueError('nonnegative root and positive integer mesh')
    k=isqrt(q.numerator*mesh*mesh//q.denominator)
    lo=F(k,mesh)
    return lo,lo if lo*lo==q else F(k+1,mesh)


def acos_lower(value,terms=12,mesh=10**10):
    x=F(value)
    if not 0<=x<=1 or type(terms)is not int or terms<2 or terms%2:
        raise ValueError('0<=x<=1 and positive even number of atan terms')
    t=sqrt_bounds((1-x)/(1+x),mesh)[0]
    return 2*sum(((-1)**j*t**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))


def radial_h_lower(radius,terms=24):
    s=F(radius)
    if not 0<s<=1 or type(terms)is not int or terms<2:raise ValueError('0<s<=1 and at least two log terms')
    q=1-s*s
    return sum((q**j/F(2*j) for j in range(2,terms+1)),F(0))


@lru_cache(maxsize=16)
def pair_integral_lower(cutoff,bins=4096,log_terms=24,atan_terms=12,mesh=10**10):
    """P(c)>=returned fraction for sqrt(2)<=c<2; includes both ordered terms."""
    c=F(cutoff)
    if c<=0 or c*c<2 or c>=2 or type(bins)is not int or not 1<=bins<=16384:
        raise ValueError('rational sqrt(2)<=cutoff<2 and 1..16384 bins')
    # Validate even if every evaluation were at a zero endpoint.
    acos_lower(F(1,2),atan_terms,mesh);radial_h_lower(F(1,2),log_terms)
    step=(2-c)/bins;total=F(0)
    for j in range(bins):
        left=c-1+j*step;right=left+step
        cosine=(c*c+left*left-1)/(2*c*left)
        require(0<=cosine<=1,'admissible lens angle')
        total+=step*left*acos_lower(cosine,atan_terms,mesh)*radial_h_lower(right,log_terms)
    return 4*total/certified_pi_upper()


FIXED_COSTS={F(3151,2000):F(1003,100000),F(8,5):F(13,1600),F(44,25):F(1,750)}


@lru_cache(maxsize=1)
def fixed_certificates():
    rows=[]
    for cutoff,cost in FIXED_COSTS.items():
        raw=pair_integral_lower(cutoff)
        require(raw>cost,'fixed exact integral lower '+str(cutoff))
        rows.append({'cutoff_force':cutoff,'accepted_cost_lower':cost,
                     'quadrature_lower':raw,'strict_margin':raw-cost})
    return rows
