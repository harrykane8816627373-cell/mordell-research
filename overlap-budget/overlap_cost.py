"""Certified parabolic-overlap lower cost for short pairs at force-one scale."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'collective-stability'))
from ring_population import pi_interval


@lru_cache(maxsize=1)
def certified_pi_upper():
    lower,upper=pi_interval()
    if not 3<lower<=upper<F(355,113):raise ArithmeticError('Machin pi enclosure')
    return F(355,113)


def pair_cost_squared_lower(cutoff):
    """P(c)^2 lower, where each UNORDERED pair contributes two ordered terms.

    Requires sqrt(2)<=actual distance<=cutoff<2. The decreasing cost can
    be evaluated at the cutoff. Constants never use floating point.
    """
    c=F(cutoff)
    if c<=0 or c*c<2 or c>=2:
        raise ValueError('rational sqrt(2)<=cutoff<2 required')
    eps=1-c/2
    return F(512**2,105**2)*2*(c-1)*eps**7/certified_pi_upper()**2


def pair_cost_lower(cutoff,scale=10**12):
    if type(scale)is not int or scale<1:raise ValueError('positive root mesh required')
    q=pair_cost_squared_lower(cutoff)
    return F(isqrt(q.numerator*scale*scale//q.denominator),scale)


def short_edge_capacity(n,proven_gap_upper,cutoff=F(8,5),dirichlet_lower=0,scale=10**12):
    """Global necessary edge count under independently proved energy inputs.

    dirichlet_lower is ONE bound on n^2D(nu), or an orthogonal radial plus
    angular projection. It is not the previous total linear-gap bound.
    """
    if type(n)is not int or n<2:raise ValueError('integer n>=2 required')
    gap,used=F(proven_gap_upper),F(dirichlet_lower)
    if gap<0 or used<0:raise ValueError('nonnegative proved bounds required')
    cost=pair_cost_lower(cutoff,scale)
    if cost<=0:raise ValueError('root mesh too coarse for positive pair cost')
    remaining=gap-used
    maximum=None if remaining<0 else min(n*(n-1)//2,(remaining/cost).numerator//(remaining/cost).denominator)
    return {'n':n,'cutoff_force':F(cutoff),'per_unordered_pair_cost_lower':cost,
            'remaining_gap':remaining,'maximum_short_edges':maximum,
            'entire_class_excluded':remaining<0,
            'scope':'Counts unordered pairs at force-one distance <=cutoff, assuming global geometry and the stated proved total-gap/Dirichlet inputs.'}


def frobenius_edge_cap(n,cutoff=F(8,5)):
    if type(n)is not int or n<2 or F(cutoff)<=0:raise ValueError('n>=2 and c>0 required')
    bound=F(n-1,2)*F(cutoff)**4
    return min(n*(n-1)//2,bound.numerator//bound.denominator)
