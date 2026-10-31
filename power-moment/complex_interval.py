"""Gaussian-rational numbers and outward rectangular complex intervals."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt


def root_bounds(value,mesh=10**24):
    q=F(value)
    if q<0 or type(mesh)is not int or mesh<1:raise ValueError('nonnegative square and integer mesh')
    k=isqrt(q.numerator*mesh*mesh//q.denominator);lo=F(k,mesh)
    return lo,lo if lo*lo==q else F(k+1,mesh)


@dataclass(frozen=True)
class C:
    re:F=F(0)
    im:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'re',F(self.re));object.__setattr__(self,'im',F(self.im))
    def __add__(self,b):
        b=as_c(b);return C(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self):return C(-self.re,-self.im)
    def __sub__(self,b):return self+-as_c(b)
    def __rsub__(self,b):return as_c(b)+-self
    def __mul__(self,b):
        b=as_c(b);return C(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=as_c(b);q=b.norm()
        if not q:raise ZeroDivisionError('zero Gaussian rational')
        return self*b.conjugate()*F(1,q)
    def __rtruediv__(self,b):return as_c(b)/self
    def __pow__(self,k):
        if type(k)is not int or k<0:raise ValueError('nonnegative integer exponent')
        out=C(1);base=self
        while k:
            if k%2:out=out*base
            base=base*base;k//=2
        return out
    def conjugate(self):return C(self.re,-self.im)
    def norm(self):return self.re*self.re+self.im*self.im
    def __bool__(self):return bool(self.re or self.im)


def as_c(x):return x if isinstance(x,C) else C(x)


MESH=10**24
def outward(lo,hi):
    lo,hi=F(lo),F(hi)
    if lo>hi:raise ValueError('ordered real interval')
    a,b=lo*MESH,hi*MESH
    return F(a.numerator//a.denominator,MESH),F(-(-b.numerator//b.denominator),MESH)


def radd(a,b):return outward(a[0]+b[0],a[1]+b[1])
def rmul(a,b):
    values=[x*y for x in a for y in b];return outward(min(values),max(values))
def rneg(a):return -a[1],-a[0]


@dataclass(frozen=True)
class I:
    re:tuple
    im:tuple=(F(0),F(0))
    def __post_init__(self):
        object.__setattr__(self,'re',tuple(map(F,self.re)));object.__setattr__(self,'im',tuple(map(F,self.im)))
        if len(self.re)!=2 or len(self.im)!=2 or self.re[0]>self.re[1] or self.im[0]>self.im[1]:raise ValueError('ordered rectangular interval')
    def __add__(self,b):
        b=as_i(b);return I(radd(self.re,b.re),radd(self.im,b.im))
    __radd__=__add__
    def __neg__(self):return I(rneg(self.re),rneg(self.im))
    def __sub__(self,b):return self+-as_i(b)
    def __rsub__(self,b):return as_i(b)+-self
    def __mul__(self,b):
        b=as_i(b)
        return I(radd(rmul(self.re,b.re),rneg(rmul(self.im,b.im))),
                 radd(rmul(self.re,b.im),rmul(self.im,b.re)))
    __rmul__=__mul__
    def __truediv__(self,b):
        if isinstance(b,C) and b.im==0:b=b.re
        b=F(b)
        if b==0:raise ZeroDivisionError('zero real scalar')
        return self*F(1,b)
    def __pow__(self,k):
        if type(k)is not int or k<0:raise ValueError('nonnegative integer exponent')
        out=as_i(1);base=self
        while k:
            if k%2:out=out*base
            base=base*base;k//=2
        return out
    def abs1_upper(self):return max(map(abs,self.re))+max(map(abs,self.im))
    def contains(self,z):
        z=as_c(z);return self.re[0]<=z.re<=self.re[1] and self.im[0]<=z.im<=self.im[1]


def as_i(x):
    if isinstance(x,I):return x
    x=as_c(x);return I((x.re,x.re),(x.im,x.im))


def inverse_det(matrix):
    n=len(matrix)
    if not n or any(len(row)!=n for row in matrix):raise ValueError('nonempty square matrix')
    A=[[as_c(x) for x in row]+[C(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    det=C(1)
    for j in range(n):
        k=next((i for i in range(j,n) if A[i][j]),None)
        if k is None:raise ZeroDivisionError('singular exact center')
        if k!=j:A[k],A[j]=A[j],A[k];det=-det
        pivot=A[j][j];det=det*pivot;A[j]=[v/pivot for v in A[j]]
        for i in range(n):
            if i==j:continue
            factor=A[i][j]
            if factor:A[i]=[x-factor*y for x,y in zip(A[i],A[j])]
    return [row[n:] for row in A],det
