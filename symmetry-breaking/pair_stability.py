"""Radial/angle pair restrictions implied by the saturated stability operator.

Only necessary conditions for local maxima; see PAIR_STABILITY.md.
"""
from fractions import Fraction as F


def pair_slack(n, radius_squared_i, radius_squared_j, distance_squared):
    if type(n) is not int or n < 2:
        raise ValueError('integer n>=2 required')
    a,b,s=map(F,(radius_squared_i,radius_squared_j,distance_squared))
    if a<0 or b<0 or s<=0:
        raise ValueError('nonnegative squared radii and positive squared distance required')
    s0=F(2*n,(n-1)**2)
    return (n-1)*(s-s0)*(1+F(n-1,2*n)*(n-1-a-b)*s)-(a-b)**2


def pair_box_upper(n, radii_i, radii_j, distance):
    """Upper bound on pair_slack over a rational box with s>=s0.

    A negative result excludes the entire box from every local maximum.
    Independent interval products intentionally overestimate correlations.
    """
    if type(n) is not int or n < 2:
        raise ValueError('integer n>=2 required')
    a0,a1=map(F,radii_i)
    b0,b1=map(F,radii_j)
    s0,s1=map(F,distance)
    if not 0<=a0<=a1 or not 0<=b0<=b1 or not F(2*n,(n-1)**2)<=s0<=s1:
        raise ValueError('ordered nonnegative radius intervals and s>=basic separation required')
    dif_lower=max(F(0),a0-b1,b0-a1)
    v0,v1=n-1-a1-b1,n-1-a0-b0
    w=[v0*s0,v0*s1,v1*s0,v1*s1]
    h0,h1=1+F(n-1,2*n)*min(w),1+F(n-1,2*n)*max(w)
    t0,t1=s0-F(2*n,(n-1)**2),s1-F(2*n,(n-1)**2)
    return (n-1)*max(t0*h0,t0*h1,t1*h0,t1*h1)-dif_lower**2
