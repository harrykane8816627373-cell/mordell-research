"""Rational radial-population and conditional energy certificates.

The analytic proof is radial_audit.md. A small discrepancy is necessary,
never sufficient, for global optimality. Regional certificates do not
constitute a full covering of the unrestricted optimization domain.
"""
from fractions import Fraction as F
from math import isqrt


def require_n(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n>=2 required')


def int_root_floor(x, degree):
    if x < 0 or degree not in (2,3,4):
        raise ValueError('invalid integer root')
    if degree == 2:
        return isqrt(x)
    if degree == 4:
        return isqrt(isqrt(x))
    lo, hi = 0, 1 << ((x.bit_length()+degree-1)//degree)
    while lo < hi:
        mid=(lo+hi+1)//2
        if mid**degree<=x:lo=mid
        else:hi=mid-1
    return lo


def root_upper(x, degree=2, scale=10**12):
    x=F(x)
    if x<0 or type(scale) is not int or scale<1:
        raise ValueError('invalid rational root')
    a=x.numerator*scale**degree;b=x.denominator
    r=int_root_floor(a//b,degree)
    if r**degree*b<a:r+=1
    return F(r,scale)


def polar_gap_upper(n, scale=10**12):
    """Existing all-n polar lower bound + disk upper, valid n>=6.

    Returns an upper bound on log(U_n/M_n), not on arbitrary candidates.
    """
    require_n(n)
    if n<6:raise ValueError('effective polar remainder uses n>=6')
    return F(n,100)+F(77,6)*root_upper(n,2,scale)+F(35,4)


def cdf_error_upper(n, gap_upper, scale=10**12):
    """Requires a separately proved nonnegative bound on log(U_n/Delta).

    For global optima one can use polar_gap_upper(n), or a tighter bound
    obtained from a proved attainable lower value.
    """
    require_n(n);gap_upper=F(gap_upper)
    if gap_upper<0:raise ValueError('nonnegative proven gap upper bound required')
    return min(F(1),root_upper(6*gap_upper/(n*n),3,scale)
               +2*root_upper(F(1,n),2,scale)+F(1,n))


def ceil_f(x):
    return -((-x.numerator)//x.denominator)


def disk_population(n, s, gap_upper, scale=10**12):
    require_n(n);s=F(s)
    if s<0:raise ValueError('nonnegative original-scale squared radius required')
    e=cdf_error_upper(n,gap_upper,scale)
    target=min(s/2,F(1))
    lo=max(0,ceil_f(n*(target-e)))
    hi=min(n,(n*(target+e)).numerator//(n*(target+e)).denominator)
    return {'lower':lo,'upper':hi,'cdf_error_upper':e,
            'scope':'Necessary population interval under the supplied proven energy-gap upper bound.'}


def annulus_population(n, a, b, gap_upper, scale=10**12):
    a,b=F(a),F(b)
    if not 0<=a<=b:raise ValueError('need 0<=a<=b')
    left=disk_population(n,a,gap_upper,scale)
    right=disk_population(n,b,gap_upper,scale)
    return {'lower':max(0,right['lower']-left['upper']),
            'upper':min(n,right['upper']-left['lower']),
            'scope':'Counts a < |z_i|^2 <= b, in the original moment-n scale.'}


def inner_radius_upper(n,gap_upper,scale=10**12):
    """Necessary upper bound on min |z_i|, in original moment-n scale."""
    require_n(n);gap_upper=F(gap_upper)
    if gap_upper<0:raise ValueError('nonnegative proven gap upper bound required')
    alpha_inverse=root_upper(F(n,n-1),2,scale)
    r=root_upper(F(2,n),2,scale)
    fourth=root_upper(gap_upper/(n*n),4,scale)
    return min(F(1),alpha_inverse*(r+2*fourth))


def upper_support_le(q, r2, a):
    """(sqrt(q)+sqrt(r2))^2 <= a, by safe rational comparisons."""
    q,r2,a=F(q),F(r2),F(a)
    if q<0 or r2<0:raise ValueError('negative squared radius')
    c=a-q-r2
    return c>=0 and c*c>=4*q*r2


def lower_support_le(q, r2, b):
    """max(sqrt(q)-sqrt(r2),0)^2 <= b."""
    q,r2,b=F(q),F(r2),F(b)
    if q<0 or r2<0 or b<0:raise ValueError('negative squared radius')
    if q<=r2:return True
    c=q+r2-b
    return c<=0 or c*c<=4*q*r2


def cell_cost(a,b,lo,hi):
    """Rational lower bound from F(s) in [lo,hi] on [a,b] subset [0,2]."""
    a,b,lo,hi=map(F,(a,b,lo,hi))
    if not 0<=a<b<=2 or not 0<=lo<=hi<=1:
        raise ValueError('invalid energy cell')
    integral=F(0)
    end=min(b,2*lo)
    if end>a:
        integral+=((end-2*lo)**3-(a-2*lo)**3)/12
    start=max(a,2*hi)
    if b>start:
        integral+=((b-2*hi)**3-(start-2*hi)**3)/12
    return integral/(2*b)


def band_certificate(n,bands,partition):
    """Conditional lower bound on log(U_n/Delta), for an entire radial class.

    `bands` contains (qlo,qhi,population) in the COULOMB scale x,
    with sum population=n. Each point is assigned to one stated band.
    Degenerate bands give exact squared-radius data. This input is a
    condition on configurations, not a claim that the class is nonempty.
    """
    require_n(n)
    bands=[(F(a),F(b),m) for a,b,m in bands]
    if any(not 0<=a<=b or type(m) is not int or m<1 for a,b,m in bands):
        raise ValueError('valid rational bands and positive integer populations required')
    if sum(m for a,b,m in bands)!=n:raise ValueError('population sum must be n')
    if not sum(m*a for a,b,m in bands)<=n-1<=sum(m*b for a,b,m in bands):
        raise ValueError('bands incompatible with prescribed quadratic moment')
    partition=list(map(F,partition))
    if len(partition)<2 or not 0<=partition[0]<partition[-1]<=2 or any(a>=b for a,b in zip(partition,partition[1:])):
        raise ValueError('strictly increasing partition in [0,2] required')
    r2=F(2,n);rows=[];total=F(0)
    for a,b in zip(partition,partition[1:]):
        lo=F(sum(m for qlo,qhi,m in bands if upper_support_le(qhi,r2,a)),n)
        hi=F(sum(m for qlo,qhi,m in bands if lower_support_le(qlo,r2,b)),n)
        value=cell_cost(a,b,lo,hi)
        total+=value
        rows.append({'a':str(a),'b':str(b),'cdf_lower':str(lo),'cdf_upper':str(hi),'energy_lower':str(value)})
    return {'n':n,'delta_lower':n*n*total,'radial_energy_lower':total,'cells':rows,
            'scope':'For every moment-n configuration whose Coulomb-scale squared radii satisfy these assigned bands; no global domain coverage is asserted.'}


def exact_radii_certificate(qs,partition):
    qs=list(map(F,qs));n=len(qs)
    if sum(qs)!=n-1:raise ValueError('Coulomb-scale radii must sum to n-1')
    return band_certificate(n,[(q,q,1) for q in qs],partition)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('n',type=int)
    p.add_argument('--s',default='1',help='original-scale squared radius for a disk population')
    args=p.parse_args()
    gap=polar_gap_upper(args.n)
    result=disk_population(args.n,F(args.s),gap)
    result.update(n=args.n,s=args.s,proven_global_gap_upper=str(gap))
    print(json.dumps(result,default=str,indent=2))
