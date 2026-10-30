"""Certified Gaussian disk means and rotation-filtered Dirichlet probes.

No wave local-smoothing theorem is used. Implementation limits are explicit;
mathematical identities are not restricted to these arithmetic budgets.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,isqrt

MESH=10**40


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def outward(lo,hi):
    lo,hi=F(lo),F(hi)
    if lo>hi:raise ValueError('ordered interval')
    a,b=lo*MESH,hi*MESH
    return F(a.numerator//a.denominator,MESH),F(-(-b.numerator//b.denominator),MESH)


def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def scale(a,x):return mul(a,(F(x),F(x)))
def square(a):
    lo=F(0) if a[0]<=0<=a[1] else min(a[0]*a[0],a[1]*a[1])
    return lo,max(a[0]*a[0],a[1]*a[1])
def distance_zero(a):return F(0) if a[0]<=0<=a[1] else min(abs(a[0]),abs(a[1]))


def root_upper(q,mesh=10**15):
    q=F(q)
    if q<0:raise ValueError('nonnegative root')
    k=isqrt(q.numerator*mesh*mesh//q.denominator)
    return F(k if F(k,mesh)**2==q else k+1,mesh)


@lru_cache(maxsize=4096)
def exp_negative(x):
    x=F(x)
    if not 0<=x<=128:raise ValueError('exponential budget 0<=x<=128')
    z=x;steps=0
    while z>1:z/=2;steps+=1
    term=value=F(1)
    for j in range(1,81):term*=z/j;value+=term
    tail=term*z/81/(1-z/82)
    lo,hi=outward(1/(value+tail),1/value)
    for _ in range(steps):lo,hi=outward(lo*lo,hi*hi)
    return lo,hi


@lru_cache(maxsize=4096)
def laplace_moment(m,y):
    """Integral_0^1 u^m exp(-y*u) du with signed Taylor remainder."""
    y=F(y)
    if type(m)is not int or not 0<=m<=128 or not 0<=y<=16:
        raise ValueError('moment order 0..128, 0<=y<=16')
    term=F(1);lower=F(0)
    for j in range(160):
        if j:term*=-y/j
        lower+=term/F(m+j+1)
    upper=lower+(-term*y/160)/F(m+161)
    # e^-yu is nonnegative even if an intentionally coarse lower were negative.
    return outward(max(F(0),lower),min(F(1,m+1),upper))


@lru_cache(maxsize=4096)
def bessel_sum(z,k,y=None,shift=0):
    """Sum z^l/(l!(l+k)!), optionally weighted by Laplace moment k+l."""
    z=F(z)
    if type(k)is not int or not 0<=k<=16 or not 0<=z<=256 or type(shift)is not int or not 0<=shift<=12:
        raise ValueError('Bessel order 0..16 and 0<=z<=256')
    y=None if y is None else F(y)
    lo=hi=F(0);weight=F(1,factorial(k))
    for l in range(65):
        if l:weight*=z/F(l*(l+k))
        if y is None:lo+=weight;hi+=weight
        else:
            a,b=laplace_moment(k+shift+l,y);lo+=weight*a;hi+=weight*b
    nxt=weight*z/F(65*(65+k));ratio=z/F(66*(66+k))
    require(ratio<1,'Bessel geometric tail ratio')
    hi+=nxt/(1-ratio)*(F(1) if y is None else F(1,k+shift+66))
    return outward(lo,hi)


@lru_cache(maxsize=8192)
def disk_gaussian(D,R2,t):
    D,R2,t=map(F,(D,R2,t))
    if D<0 or R2<=0 or t<=0:raise ValueError('nonnegative distance, positive smear/time')
    x,y=D/(4*t),R2/(4*t)
    return outward(*mul(exp_negative(x),bessel_sum(x*y,0,y)))


def cmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def cpow(a,k):
    out=(F(1),F(0))
    for _ in range(k):out=cmul(out,a)
    return out
def rotate(a):return -a[1],a[0]
def norm(a):return a[0]*a[0]+a[1]*a[1]


def quarter_probe(a,b,t):
    a,b=tuple(map(F,a)),tuple(map(F,b));t=F(t)
    if len(a)!=2 or len(b)!=2 or norm(a)!=norm(b) or norm(a)<=0 or t<=0:
        raise ValueError('two nonzero equal-radius rational centers and t>0')
    centers=[];weights=[]
    for c,w in [(a,1),(b,-1)]:
        for _ in range(4):centers.append(c);weights.append(w);c=rotate(c)
    if len(set(centers))!=8:raise ValueError('probe orbits must be distinct')
    return {'centers':centers,'weights':weights,'t':t,'center_squared_radius':norm(a)}


def standard_probe(t=F(1,16),radius=F(1,2)):
    r=F(radius)
    return quarter_probe((r,0),(r*F(3,5),r*F(4,5)),t)


def center_power(probe,k):
    terms=[(w,cpow(a,k)) for a,w in zip(probe['centers'],probe['weights'])]
    return tuple(sum((w*z[j] for w,z in terms),F(0)) for j in range(2))


def gaussian_gram(a,b,ti,tj):
    a,b=tuple(map(F,a)),tuple(map(F,b));ti,tj=F(ti),F(tj)
    if len(a)!=2 or len(b)!=2 or min(ti,tj)<=0:raise ValueError('two-dimensional centers and positive Gaussian times')
    d=sum(((x-y)**2 for x,y in zip(a,b)),F(0));s=ti+tj;x=d/(4*s)
    return scale(exp_negative(x),2*ti*tj/(s*s)*(1-x))


def probe_gram(left,right):
    result=(F(0),F(0))
    for a,wa in zip(left['centers'],left['weights']):
        for b,wb in zip(right['centers'],right['weights']):
            result=add(result,scale(gaussian_gram(a,b,left['t'],right['t']),wa*wb))
    return outward(*result)


def measurement(points,probe):
    points=[tuple(map(F,p)) for p in points]
    n=len(points)
    if n<2 or any(len(p)!=2 for p in points):raise ValueError('at least two two-dimensional input points')
    result=(F(0),F(0))
    for p in points:
        for a,w in zip(probe['centers'],probe['weights']):
            d=sum(((F(x)-y)**2 for x,y in zip(p,a)),F(0))
            result=add(result,scale(disk_gaussian(d,F(1,n),probe['t']),w))
    return outward(*result)


def lipschitz_upper(probe):
    return sum(map(abs,probe['weights']))*root_upper(1/(2*probe['t']))
