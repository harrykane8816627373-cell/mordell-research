"""Independent static Gaussian-probe bridge; no local-smoothing theorem input.

All acceptance arithmetic is Fraction/integer, including actual disk means.
The input is a normalized rational example, not an asserted optimizer.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
MESH=10**30


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def outward(lo,hi):
    a=lo*MESH;b=hi*MESH
    return F(a.numerator//a.denominator,MESH),F(-(-b.numerator//b.denominator),MESH)


@lru_cache(maxsize=256)
def exp_negative(x):
    x=F(x)
    if x<0 or x>10:raise ValueError('0<=x<=10 required by this example')
    value=F(1);term=F(1)
    for k in range(1,81):term*=x/k;value+=term
    nxt=term*x/81;tail=nxt/(1-x/82)
    return outward(1/(value+tail),1/value)


@lru_cache(maxsize=256)
def disk_gaussian(distance_squared,smear_squared,t):
    D,R2,t=map(F,(distance_squared,smear_squared,t));x=D/(4*t);y=R2/(4*t)
    if D<0 or R2<=0 or t<=0 or y>1:raise ValueError('positive data and smear_squared<=4t')
    z=x*y;lo=hi=F(0)
    for k in range(25):
        weight=z**k/F(factorial(k)**2)
        # Odd and even Taylor polynomials for exp(-y*u), integrated exactly.
        lower=sum(((-y)**j/F(factorial(j)*(k+j+1)) for j in range(20)),F(0))
        upper=lower+y**20/F(factorial(20)*(k+21))
        lo+=weight*lower;hi+=weight*upper
    ratio=z/F(26**2)
    require(ratio<1,'Bessel tail geometric ratio')
    hi+=z**25/F(factorial(25)**2*26)/(1-ratio)
    el,eu=exp_negative(x)
    return outward(el*lo,eu*hi)


def gaussian_gram(a,b,ti,tj):
    r2=sum(((F(x)-F(y))**2 for x,y in zip(a,b)),F(0));s=ti+tj;x=r2/(4*s)
    factor=2*ti*tj/(s*s)*(1-x);el,eu=exp_negative(x)
    return min(factor*el,factor*eu),max(factor*el,factor*eu)


def log_interval(x,terms=96):
    x=F(x)
    if x<=0:raise ValueError('positive log')
    power=0
    while x>=2:x/=2;power+=1
    while x<1:x*=2;power-=1
    def atanh_log(y):
        y=F(y)
        u=(y-1)/(y+1)
        low=2*sum((u**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
        return low,low+2*u**(2*terms+1)/(F(2*terms+1)*(1-u*u))
    a,b=atanh_log(x);l2,u2=atanh_log(2)
    return (a+power*l2,b+power*u2) if power>=0 else (a+power*u2,b+power*l2)


def rotate(a):return -a[1],a[0]


def run():
    points=[(F(0),F(0)),(F(7,9),0),(F(-7,9),0),(0,F(7,9)),(0,F(-7,9))]
    points += [(F(s*4,9),F(r*4,9)) for s in (-1,1) for r in (-1,1)]
    points=[tuple(map(F,p)) for p in points];n=len(points)
    require(n==9 and sum(x*x+y*y for x,y in points)==4,'actual q normalization')
    require(max(x*x+y*y for x,y in points)<=F(64,81),'listed global radius condition')
    require(min((x-u)**2+(y-v)**2 for i,(x,y) in enumerate(points) for u,v in points[i+1:])>=F(2,9),'listed global separation condition')
    def cmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
    powers=[]
    for k in (1,2,3):
        total=[F(0),F(0)]
        for p in points:
            z=(F(1),F(0))
            for _ in range(k):z=cmul(z,p)
            total[0]+=z[0];total[1]+=z[1]
        require(total==[0,0],'vanishing low angular power sum');powers.append(total)
    centers=[];weights=[]
    for center,weight in [((F(1,2),F(0)),1),((F(3,10),F(2,5)),-1)]:
        for _ in range(4):centers.append(center);weights.append(weight);center=rotate(center)
    require(all(x*x+y*y==F(1,4) for x,y in centers) and sum(weights)==0,'zero radial component')
    require(all(weights[centers.index(rotate(a))]==w for a,w in zip(centers,weights)),'quarter-turn invariance')
    t=F(1,16);R2=F(1,9);ml=mu=F(0)
    for p in points:
        for a,w in zip(centers,weights):
            D=sum(((x-y)**2 for x,y in zip(p,a)),F(0))
            lo,hi=disk_gaussian(D,R2,t)
            ml+=min(w*lo,w*hi);mu+=max(w*lo,w*hi)
    hl=hu=F(0)
    for i,a in enumerate(centers):
        for j,b in enumerate(centers):
            lo,hi=gaussian_gram(a,b,t,t);w=weights[i]*weights[j]
            hl+=min(w*lo,w*hi);hu+=max(w*lo,w*hi)
    require(hl>0 and (ml>0 or mu<0),'nonzero certified localized angular signal')
    min_measurement=min(abs(ml),abs(mu));energy=min_measurement**2/hu
    require(energy>0,'strict new angular-probe energy')
    require(energy>F(611,12500),'simple strict angular-energy lower')
    # Independent true product gap: z=(3/2)q, sum |z|^2=9.
    logdl=logdu=F(0)
    for i,(x,y) in enumerate(points):
        for u,v in points[i+1:]:
            d=F(9,4)*((x-u)**2+(y-v)**2);l,h=log_interval(d);logdl+=l;logdu+=h
    l9,u9=log_interval(9);lr,ur=log_interval(F(9,4))
    gap_lower=F(9,2)*l9+36*lr-18-logdu
    gap_upper=F(9,2)*u9+36*ur-18-logdl
    require(energy<gap_lower,'probe lower is below independently computed actual total gap')
    return {'passed':True,'q_points':points,'original_z_scale':'3/2',
        'q_power_sums_1_2_3':powers,'centers':centers,'weights':weights,
        'gaussian_time':t,'actual_smear_squared':R2,
        'measurement_interval':outward(ml,mu),'Dirichlet_norm_interval':outward(hl,hu),
        'angular_energy_lower':outward(energy,energy)[0],
        'simple_strict_angular_energy_lower':'611/12500',
        'passes_listed_radius_and_separation':True,
        'independent_actual_total_gap_interval':outward(gap_lower,gap_upper),
        'orthogonal_modes':'All radial tests and all angular k=2,3 tests; quarter-turn invariance leaves only nonzero multiples of 4.',
        'source_scope':'The example is normalized and has vanishing low angular moments, but is not asserted stationary or optimal. No local-smoothing theorem is used in this static certificate.'}


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x


if __name__=='__main__':
    report=encode(run())
    (HERE/'gaussian_bridge_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS gaussian_bridge_verification.json')
    for key in ['measurement_interval','Dirichlet_norm_interval','angular_energy_lower','independent_actual_total_gap_interval']:
        val=report[key];print(key,[float(F(v)) for v in val] if isinstance(val,list) else float(F(val)))
