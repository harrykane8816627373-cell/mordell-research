"""Exact finite moment and Coulomb-deficit certificates; see README.md.

The universal radius is only for global maximizers. Arbitrary shapes get
their own rigorously enclosing radius. All returns use Fraction arithmetic.
"""
from fractions import Fraction as F
from math import isqrt,prod

def root_up(a,scale=10**12):
    a=F(a)
    if a<0:raise ValueError('nonnegative radicand required')
    m=isqrt(a.numerator*scale*scale//a.denominator)
    if m*m*a.denominator<a.numerator*scale*scale:m+=1
    return F(m,scale)

def global_radius_squared(n):
    if type(n)is not int or n<2:raise ValueError('n>=2 integer required')
    return 2*(F(n-1,n)+root_up(F(1,n)))**2

def universal_discount(n):
    """Positive rational c with Mn<=Un exp(-c); not a sharp value."""
    qr=global_radius_squared(n);t=F(n,n-1)*qr
    return F(4,9)/sum((t**k/F(k)for k in range(2,n+1)),F(0))

def cmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def norm(a):return a[0]*a[0]+a[1]*a[1]
def centered(points):
    v=[tuple(map(F,p))for p in points];n=len(v)
    if n<2 or len(set(v))!=n:raise ValueError('at least 2 distinct points')
    av=tuple(sum(p[j]for p in v)/n for j in range(2))
    return [(p[0]-av[0],p[1]-av[1])for p in v]

def power_sums(points,K):
    out=[];powers=[(F(1),F(0))for _ in points]
    for k in range(1,K+1):
        powers=[cmul(a,b)for a,b in zip(powers,points)]
        out.append(tuple(sum(p[j]for p in powers)for j in range(2)))
    return out

def shape_certificate(points,K=None):
    v=centered(points);n=len(v);N=n*(n-1)//2;Q=sum(map(norm,v))
    K=n if K is None else K
    if type(K)is not int or K<1:raise ValueError('positive degree')
    scale2=F(n-1)/Q
    qs=[scale2*norm(p)for p in v]
    qr=max(F(2),(root_up(max(qs))+root_up(F(2,n)))**2)
    moments=power_sums(v,K)
    angular=sum((scale2**k*norm(p)/(k*qr**k)for k,p in enumerate(moments,1)),F(0))
    delta=prod(norm((a[0]-b[0],a[1]-b[1]))for i,a in enumerate(v)for b in v[i+1:])
    return {'n':n,'Q':Q,'R2':qr,'angular_delta_lower':angular,
            'coulomb_radii_squared':qs,'normalized_product':(F(n)/Q)**N*delta,
            'power_sums_raw':moments}

def strip_global_penalty(n,width_squared):
    w2=F(width_squared)
    if w2<0:raise ValueError('nonnegative squared half-width')
    s=max(F(0),n-1-2*n*w2)
    return s*s/(2*global_radius_squared(n)**2)

def ln_unit_interval(q,terms=60):
    """0<=ln q bracket, q in [1,2], exact positive atanh tail."""
    q=F(q)
    if not 1<=q<=2:raise ValueError('argument outside [1,2]')
    t=(q-1)/(q+1);t2=t*t;power=t;lo=F(0)
    for j in range(terms):lo+=2*power/(2*j+1);power*=t2
    tail=2*power/((2*terms+1)*(1-t2))
    return lo,lo+tail

def ln_int(v):
    if type(v)is not int or v<1:raise ValueError('positive integer')
    k=v.bit_length()-1;scale=1<<70
    a=v*scale//(1<<k)
    b=a+int(a*(1<<k)!=v*scale)
    l2,h2=ln_unit_interval(F(2))
    lo=ln_unit_interval(F(a,scale))[0]
    hi=ln_unit_interval(F(b,scale))[1]
    return k*l2+lo,k*h2+hi

def log_interval(q):
    q=F(q)
    if q<=0:raise ValueError('positive rational')
    a,b=ln_int(q.numerator);c,d=ln_int(q.denominator)
    return a-d,b-c

def deficit_interval(n,attainable):
    N=n*(n-1)//2
    a,b=log_interval(n);c,d=log_interval(F(2*n,n-1));e,f=log_interval(attainable)
    return F(n,2)*a+N*c-F(N,2)-f,F(n,2)*b+N*d-F(N,2)-e

def conservative_decimal(q,digits=8,upper=False):
    """Printable outward decimal, never used to establish a certificate."""
    q=F(q);s=10**digits;a=q.numerator*s//q.denominator
    if upper and a*q.denominator<q.numerator*s:a+=1
    return ('-'if a<0 else '')+str(abs(a)//s)+'.'+str(abs(a)%s).zfill(digits)

if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser();p.add_argument('n',type=int);args=p.parse_args()
    n=args.n;c=universal_discount(n)
    print(json.dumps({'n':n,'R_squared_upper':str(global_radius_squared(n)),
        'positive_rational_exponent_discount':str(c),'decimal_lower':conservative_decimal(c,16),
        'scope':'Mn <= Un exp(-c); c need not improve other previously known upper bounds.'},indent=2))
