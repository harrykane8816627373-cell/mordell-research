"""Effective microscopic covering of interior disks for global optimizers.

Proof: MICROSCOPIC_COVERING.md. No point-set search or unknown threshold.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,isqrt


def exp_interval(x,terms=24):
    x=F(x)
    if x<0 or type(terms)is not int or terms<1:
        raise ValueError('nonnegative rational and positive series length required')
    reductions=0
    while x>1:
        x/=2
        reductions+=1
    lower=sum((x**k/F(factorial(k)) for k in range(terms+1)),F(0))
    tail=x**(terms+1)/factorial(terms+1)/(1-x/F(terms+2))
    upper=lower+tail
    for _ in range(reductions):
        lower,upper=lower*lower,upper*upper
    return lower,upper


def sqrt_interval(x,scale=10**12):
    x=F(x)
    if x<0 or type(scale)is not int or scale<1:
        raise ValueError('nonnegative rational and positive root scale required')
    numerator,denominator=x.numerator*scale*scale,x.denominator
    root=isqrt(numerator//denominator)
    return F(root,scale),F(root+int(root*root*denominator<numerator),scale)


def ceil_fraction(x):
    x=F(x)
    return -((-x.numerator)//x.denominator)


@lru_cache(maxsize=32)
def parameters(bulk_squared=F(1),tail_ratio=F(5,6),threshold=64):
    """Validate one uniform bulk radius/tail/starting-count certificate.

    Every n>=threshold is covered. Generic choices are accepted only after
    exact exponent, monotonicity and tail-margin validation.
    """
    bulk_squared,q=F(bulk_squared),F(tail_ratio)
    if not 0<bulk_squared<2 or not 0<q<1 or type(threshold)is not int or threshold<4 or threshold%2:
        raise ValueError('0<bulk radius squared<2, 0<q<1, and even threshold>=4 required')
    beta=bulk_squared/2
    chernoff_upper=beta*exp_interval(1-beta)[1]
    if chernoff_upper>q:
        raise ValueError('tail ratio does not bound beta*exp(1-beta)')
    m=threshold//2
    if q*(2*m+3)>2*m+1:
        raise ValueError('tail envelope not proved decreasing above threshold')
    tail=(2*m+1)*q**m+q**(2*m)
    mesh=10**12
    tail_upper=F(ceil_fraction(tail*mesh),mesh)
    if tail_upper>=F(1,4):
        raise ValueError('Poisson-tail margin too large')
    gaussian_upper=8/exp_interval(F(23,8),12)[0]
    if gaussian_upper>=F(1,2):
        raise ArithmeticError('packing Gaussian tail margin failed')
    return {'bulk_squared':bulk_squared,'tail_ratio':q,'threshold':threshold,
            'chernoff_base_upper':chernoff_upper,
            'poisson_error_upper':tail_upper,'gaussian_sum_upper':gaussian_upper,
            'combined_error_upper':tail_upper+gaussian_upper,
            'scaled_hole_radius':3,
            'scope':'All global optimizers for every n above threshold and every center in the stated bulk disk.'}


TABLE={F(1):(F(5,6),64),F(3,2):(F(27,28),512),F(7,4):(F(127,128),4096)}


def covering_certificate(n,bulk_squared=F(1),scale=10**12):
    if type(n)is not int or n<2:
        raise ValueError('integer n>=2 required')
    bulk_squared=F(bulk_squared)
    if bulk_squared not in TABLE:
        raise ValueError('choose a fixed certified bulk, or validate generic parameters separately')
    q,start=TABLE[bulk_squared]
    data=parameters(bulk_squared,q,start)
    if n<start:
        raise ValueError('n is below this bulk certificate threshold')
    root_lower=sqrt_interval(n-1,scale)[0]
    radius_upper=6/root_lower
    return {'n':n,'bulk_squared':bulk_squared,'starting_count':start,
            'covering_radius_strict_upper':radius_upper,
            'exact_covering_radius_squared_strict_upper':F(36,n-1),
            'parameters':data,
            'scope':'For every |u|^2<=bulk_squared some point satisfies |z_i-u|<6/sqrt(n-1), in original moment-n scale.'}


def disk_population(n,center,radius,bulk_squared=F(1),scale=10**12):
    """Microscopic positive-density interval for an entire disk in the bulk.

    Covering gives a lower count, strict global separation an upper count.
    This does not claim sharp uniform-disk population discrepancy.
    """
    cert=covering_certificate(n,bulk_squared,scale)
    if len(center)!=2:
        raise ValueError('rational center pair required')
    x,y=map(F,center)
    radius=F(radius)
    if radius<0:
        raise ValueError('nonnegative disk radius required')
    center_q=x*x+y*y
    remainder=F(bulk_squared)-center_q-radius*radius
    if remainder<0 or remainder*remainder<4*center_q*radius*radius:
        raise ValueError('whole disk must lie in the certified bulk')
    h=cert['covering_radius_strict_upper']
    lower_real=max(F(0),radius/h-1)**2
    root_upper=sqrt_interval(n-1,scale)[1]
    upper_real=(1+radius*root_upper)**2
    return {'n':n,'center':(x,y),'radius':radius,'bulk_squared':F(bulk_squared),
            'lower':min(n,ceil_fraction(lower_real)),
            'upper':min(n,upper_real.numerator//upper_real.denominator),
            'covering_lower_real':lower_real,'packing_upper_real':upper_real,
            'scope':'Necessary count in the closed disk, under global optimality; constants are coarse and do not determine the point arrangement.'}
