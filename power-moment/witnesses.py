"""Exact centered/normalized point witnesses for moment-region nonemptiness."""
from fractions import Fraction as F
from complex_interval import C,I,as_i,root_bounds


def normalized(raw,mass):
    raw=list(raw);n=len(raw);mass=F(mass)
    if n<2 or mass<=0 or len(set(raw))!=n:raise ValueError('distinct points and positive target mass')
    average=sum(raw,C(0))/n;raw=[z-average for z in raw]
    q=sum((z.norm() for z in raw),F(0))
    if not q:raise ValueError('nonzero centered second moment')
    ratio=mass/q;sqrt=I(root_bounds(ratio));powers=[as_i(n),as_i(0)]
    for k in range(2,n+1):
        value=as_i(sum((z**k for z in raw),C(0)))*ratio**(k//2)
        if k%2:value=value*sqrt
        powers.append(value)
    product=F(1)
    for i in range(n):
        for j in range(i):product*=(raw[i]-raw[j]).norm()*ratio
    return raw,ratio,powers,product
