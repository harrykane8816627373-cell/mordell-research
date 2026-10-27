"""Standard-library exact checks for the effective second-order grid theorem.

The general analytic theorem is proved in GRID_PROOF.md. These checks
exercise exact normalization, force cancellation, log intervals, and constants.
"""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
import json
checks=0
def require(v,name):
    global checks
    if not v:raise ArithmeticError(name)
    checks+=1
def add(z,w):return(z[0]+w[0],z[1]+w[1])
def neg(z):return(-z[0],-z[1])
def sub(z,w):return add(z,neg(w))
def scale(z,a):return(z[0]*a,z[1]*a)
def mul(z,w):return(z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def norm(z):return z[0]**2+z[1]**2
def div(z,w):return scale(mul(z,(w[0],-w[1])),1/norm(w))
def total(v):
    out=(F(0),F(0))
    for a in v:out=add(out,a)
    return out
def center(v):
    avg=scale(total(v),F(1,len(v)))
    return[sub(z,avg)for z in v]
def delta(v):
    out=F(1)
    for i,j in combinations(range(len(v)),2):out*=norm(sub(v[i],v[j]))
    return out
def log_interval(x,order=16):
    require(x>0,'positive logarithm argument')
    y=(x-1)/(x+1);p=y;acc=F(0)
    for j in range(order):
        acc+=2*p/(2*j+1);p*=y*y
    tail=2*abs(p)/((2*order+1)*(1-y*y))
    return acc-tail,acc+tail
def ceil_sqrt(x):
    q=isqrt(x.numerator//x.denominator)
    return q if q*q>=x else q+1

for n in range(3,101):
    nn=F(n);lam=F(n-1,2);NN=F(n*(n-1),2)
    require(NN==nn*lam,'homogeneous scale identity')
    for eta in [F(1,6),F(1,3),F(1),F(2)]:
        sigma=6*eta/(n+6*eta)
        require(n*sigma/(6*(1-sigma))==eta,'exact eta remainder constant')
        threshold=F(n-1,2)/sigma**2
        require(threshold==F((n-1)*(n+6*eta)**2,72*eta**2),'eta q threshold')
    q=ceil_sqrt(F((n-1)*(n+1)**2,2))
    sigma=F(1,n+1)
    require(2*q*q*sigma*sigma>=n-1,'rational t enclosure')
    require(1+n*sigma/(6*(1-sigma))==F(7,6),'7/6 constant')
    for eps in[F(1,10),F(1,2),F(2)]:
        K=F(7*n*(n-1),12)
        qq=ceil_sqrt(max(F((n-1)*(n+1)**2,2),K*(1+eps)/eps))
        rr=K/(qq*qq)
        require(0<rr<1,'positive rational loss below one')
        require(1/(1-rr)<=1+eps,'exact rational stopping envelope')

for n in range(3,11):
    raw=[(F((i*i+3*i)%13-6,37),F((i*i*i+2*i)%11-5,41))for i in range(n)]
    ec=center(raw);E=sum(norm(z)for z in ec)
    pair=sum(norm(sub(ec[i],ec[j]))for i,j in combinations(range(n),2))
    require(pair==n*E,'centered pair variance identity')
    require(E<=sum(norm(z)for z in raw),'centering reduces error norm')

# A physical four-point square, with exact Gaussian-rational perturbations.
z=[(F(1),F(0)),(F(0),F(1)),(F(-1),F(0)),(F(0),F(-1))]
require(delta(z)==256,'square exact Vandermonde value')
for seed in range(1,21):
    raw=[(F((seed*(i+1)+i*i)%11-5,100),F((seed*(i+2)+3*i)%13-6,100))for i in range(4)]
    e=center(raw);E=sum(norm(a)for a in e);w=[add(a,b)for a,b in zip(z,e)]
    pp=[div(sub(e[i],e[j]),sub(z[i],z[j]))for i,j in combinations(range(4),2)]
    first=total(pp)[0]
    dot=sum(mul(ei,(zi[0],-zi[1]))[0]for ei,zi in zip(e,z))
    require(first==F(3,2)*dot,'exact force linear cancellation')
    a2=total([mul(v,v)for v in pp])[0]
    require(abs(a2)<=F(3,2)*E,'second order real quadratic bound')
    sigma=F(1,4)
    require(all(norm(v)<sigma*sigma for v in pp),'strict perturbation disk')
    Q=sum(norm(v)for v in w)
    ratio=F(256)* (Q/4)**6/delta(w)
    lo,hi=log_interval(ratio)
    bound=3*(1+4*sigma/(6*(1-sigma)))*E
    require(ratio>=1,'square normalized loss is nonnegative')
    require(hi<=bound,'rigorous log interval within second order bound')

report={'passed':True,'exact_checks':checks,'square_perturbation_instances':20,'arithmetic':'fractions.Fraction; no floating point','scope':'Finite implementation and constant checks; the all-n analytic proof is GRID_PROOF.md.'}
(Path(__file__).parent/'second_order_grid_verification.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps(report,indent=2))
