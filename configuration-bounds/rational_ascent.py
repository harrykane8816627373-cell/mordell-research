"""Exact constructive ascent for rational imaginary-quadratic point shapes.

A point (a,b) represents a+i*sqrt(d)*b, d a positive integer. All squared
norms and pair products are rational. This routine proves an improvement
when the projected force residual is nonzero; residual zero is not a
global optimality certificate.
"""
from fractions import Fraction as F
from itertools import combinations
from math import prod


def add(z,w):return (z[0]+w[0],z[1]+w[1])
def sub(z,w):return (z[0]-w[0],z[1]-w[1])
def scale(z,a):return (z[0]*a,z[1]*a)
def conjugate(z):return (z[0],-z[1])
def mul(z,w,d):return (z[0]*w[0]-d*z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def norm(z,d):return z[0]*z[0]+d*z[1]*z[1]
def divide(z,w,d):return scale(mul(z,conjugate(w),d),1/norm(w,d))


def centered(points,d=1):
    if type(d) is not int or d<1 or len(points)<2:
        raise ValueError('positive integer d and at least two points required')
    v=[(F(a),F(b))for a,b in points]
    if len(set(v))!=len(v):raise ValueError('distinct points required')
    n=len(v)
    mean=(sum(p[0]for p in v)/n,sum(p[1]for p in v)/n)
    return tuple(sub(p,mean)for p in v)


def value(points,d=1):
    v=centered(points,d);n=len(v);N=n*(n-1)//2
    Q=sum(norm(p,d)for p in v)
    return (F(n)/Q)**N*prod(norm(sub(a,b),d)for a,b in combinations(v,2))


def residual_data(points,d=1):
    v=centered(points,d);n=len(v);N=n*(n-1)//2
    Q=sum(norm(p,d)for p in v)
    force=[]
    for i,p in enumerate(v):
        f=(F(0),F(0))
        for j,q in enumerate(v):
            if i!=j:f=add(f,divide((F(1),F(0)),sub(p,q),d))
        force.append(f)
    residual=tuple(sub(conjugate(f),scale(p,F(N)/Q))for f,p in zip(force,v))
    R=sum(norm(r,d)for r in residual)
    ratios=[divide(sub(residual[i],residual[j]),sub(v[i],v[j]),d)
            for i,j in combinations(range(n),2)]
    A=sum(norm(w,d)for w in ratios)
    # |a+i sqrt(d)b| <= |a|+d|b| for integer d>=1.
    C=max(abs(a)+d*abs(b)for a,b in ratios)
    H=2*A+F(N)*R/Q
    return {'points':v,'d':d,'Q':Q,'N':N,'force':tuple(force),
            'residual':residual,'R':R,'normalized_residual':Q*R,
            'ratios':tuple(ratios),'A':A,'C':C,'H':H}


def ascent_step(points,d=1):
    data=residual_data(points,d)
    before=value(data['points'],d)
    R,C,H=data['R'],data['C'],data['H']
    if R==0:
        return {**data,'stationary':True,'step':F(0),'log_gain_lower':F(0),
                'before':before,'after':before,'new_points':data['points']}
    denominator=H+2*C*R
    t=R/denominator
    g=t*R
    out=tuple(add(p,scale(r,t))for p,r in zip(data['points'],data['residual']))
    after=value(out,d)
    if not 0<t*C<F(1,2) or not after>before or not after/before>1+g:
        raise ArithmeticError('exact ascent certificate failed')
    return {**data,'stationary':False,'step':t,'log_gain_lower':g,
            'before':before,'after':after,'new_points':out}


def stopping_budget(points,epsilon,d=1):
    """A theoretical finite update budget for Q*R<=epsilon, not execution.

    The bound is deliberately coarse and can be extremely large.
    It does not guarantee global objective accuracy.
    """
    v=centered(points,d);n=len(v);N=n*(n-1)//2;epsilon=F(epsilon)
    if epsilon<=0:raise ValueError('positive rational tolerance required')
    F0=value(v,d)
    a=F0/F(n**N*2**(N-1))
    B=F(N)+(2*n+2*n*(n-1)*(d+1))/a
    U=F(2*n,n-1)**N
    limit=B*(U/F0-1)/epsilon
    budget=(limit.numerator+limit.denominator-1)//limit.denominator
    return {'F0':F0,'normalized_separation_floor':a,'gain_denominator':B,
            'global_upper':U,'residual_tolerance':epsilon,'update_budget':budget,
            'search_executed':False,
            'scope':'Stationarity tolerance only, not global optimality.'}


def iterate(points,epsilon,d=1,max_steps=100):
    """Run at most max_steps certified updates, without claiming completion.

    The status explicitly distinguishes tolerance attainment from a budget
    stop. Every returned update is a strictly improved rational bound.
    """
    epsilon=F(epsilon)
    if epsilon<=0 or type(max_steps)is not int or max_steps<0:
        raise ValueError('positive tolerance and nonnegative integer budget required')
    v=centered(points,d);rows=[]
    for k in range(max_steps+1):
        data=residual_data(v,d)
        if data['normalized_residual']<=epsilon:
            return {'status':'stationarity_tolerance_met','points':v,'rows':rows,
                    'normalized_residual':data['normalized_residual'],'value':value(v,d)}
        if k==max_steps:break
        row=ascent_step(v,d);rows.append(row);v=row['new_points']
    return {'status':'update_budget_exhausted','points':v,'rows':rows,
            'normalized_residual':data['normalized_residual'],'value':value(v,d)}
