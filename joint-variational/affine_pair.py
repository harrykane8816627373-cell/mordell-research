"""Global rational affine minorants of the exact unordered pair loss."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'layer-circumference'))
from pair_integral import sqrt_bounds,acos_lower,pair_integral_lower,radial_h_lower
sys.path.insert(0,str(ROOT/'collective-stability'))
from ring_population import pi_interval


def require(ok,label):
    if not ok:raise ArithmeticError(label)


@lru_cache(maxsize=1)
def pi_bounds():
    lo,hi=pi_interval();a=F(3141592653589793,10**15);b=F(3141592653589794,10**15)
    require(a<lo<=hi<b,'Machin fixed pi bounds')
    return a,b


def acos_upper(x,terms=13,mesh=10**12):
    x=F(x)
    if not 0<=x<=1 or type(terms)is not int or terms<1 or terms%2!=1:
        raise ValueError('0<=x<=1 and positive odd atan length')
    t=sqrt_bounds((1-x)/(1+x),mesh)[1]
    require(t<=1,'atan alternating domain')
    return 2*sum(((-1)**j*t**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))


def tail_interval(d):
    """P'(d)=-2*tau(d)/d; tau is the two-disk convolution tail."""
    d=F(d)
    if d<=0 or d*d<2 or d>2:raise ValueError('sqrt(2)<=d<=2')
    al,au=acos_lower(d/2,12,10**12),acos_upper(d/2)
    bl,bu=sqrt_bounds(1-d*d/4,10**12)
    # The coefficient 1-d^2 is negative: reverse the angle bounds.
    vl=(1-d*d)*au+(d/2)*(1+d*d/2)*bl
    vu=(1-d*d)*al+(d/2)*(1+d*d/2)*bu
    pl,pu=pi_bounds()
    return max(F(0),2*vl/pu),min(F(1),2*vu/pl)


@lru_cache(maxsize=16)
def support_line(left,right,slope,bins=4096):
    left,right,beta=F(left),F(right),F(slope)
    if left<=0 or left*left<2 or not left<right<2 or beta<=0:
        raise ValueError('sqrt(2)<=left<right<2 and beta>0')
    tl=tail_interval(left)[0];tu=tail_interval(right)[1]
    require(2*tl/left>beta>2*tu/right,'bracket contains the minimizing derivative')
    raw=pair_integral_lower(right,bins)+beta*left
    alpha=F((raw*10**12).numerator//(raw*10**12).denominator,10**12)
    return {'left':left,'right':right,'beta':beta,'alpha':alpha,
            'left_derivative_magnitude_lower':2*tl/left,
            'right_derivative_magnitude_upper':2*tu/right,
            'raw_intercept_lower':raw,
            'scope':'P(d)>=alpha-beta*d for every d>0, with P(d)=0 for d>=2. Each unordered pair includes both ordered terms.'}


def pair_integral_upper(c,bins=4096):
    c=F(c)
    if c<=0 or c*c<2 or c>=2 or type(bins)is not int or not 1<=bins<=16384:
        raise ValueError('sqrt(2)<=c<2 and 1..16384 bins')
    step=(2-c)/bins;total=F(0);rounding_mesh=10**18
    for j in range(bins):
        l=c-1+j*step;u=l+step
        cosine=(c*c+u*u-1)/(2*c*u)
        q=1-l*l
        h=radial_h_lower(l)+q**25/(50*(1-q))
        value=step*u*acos_upper(cosine)*h
        scaled=value*rounding_mesh
        total+=F(-(-scaled.numerator//scaled.denominator),rounding_mesh)
    return 4*total/pi_bounds()[0]


def fixed_line():return support_line(F(17901,10000),F(8951,5000),F(7,500))
