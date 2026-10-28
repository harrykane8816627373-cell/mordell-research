"""Exact checks of the constructive rational ascent theorem."""
from fractions import Fraction as F
from itertools import combinations
from math import lcm, prod
from pathlib import Path
import hashlib
import json

from rational_ascent import (add,sub,scale,conjugate,mul,norm,centered,value,
                             residual_data,ascent_step,stopping_budget,iterate)


def require(ok,name):
    if not ok:raise ArithmeticError(name)


def independent_integer_value(points,d):
    n=len(points);N=n*(n-1)//2
    common=lcm(*(q.denominator for p in points for q in p))
    raw=[(int(a*common),int(b*common))for a,b in points]
    sx=sum(x for x,y in raw);sy=sum(y for x,y in raw)
    D=n*sum(x*x+d*y*y for x,y in raw)-sx*sx-d*sy*sy
    distances=[(x-u)**2+d*(y-v)**2 for (x,y),(u,v) in combinations(raw,2)]
    return F(n*n,D)**N*prod(distances)


def log_lower(x,terms=24):
    require(x>=1,'positive logarithm ratio')
    y=(x-1)/(x+1)
    return 2*sum((y**(2*j+1)/(2*j+1)for j in range(terms)),F(0))


def first_points(n,triangular=False):
    q=lambda p:p[0]*p[0]+(p[0]*p[1]if triangular else 0)+p[1]*p[1]
    points=sorted(((a,b)for a in range(-n,n+1)for b in range(-n,n+1)),key=lambda p:(q(p),p))[:n]
    if triangular:return[(F(a)+F(b,2),F(b,2))for a,b in points]
    return[(F(a),F(b))for a,b in points]


def digest(q):
    def ih(x):return hashlib.sha256(x.to_bytes(max(1,(x.bit_length()+7)//8),'big')).hexdigest()
    return {'numerator_sha256':ih(q.numerator),'denominator_sha256':ih(q.denominator)}


def main():
    stationary=[];improved=[];rows=[];checks=0
    for d,triangular in ((1,False),(3,True)):
        for n in range(2,21):
            v=first_points(n,triangular)
            row=ascent_step(v,d)
            points=row['points'];res=row['residual'];Q=row['Q'];R=row['R']
            require(sum(p[0]for p in points)==sum(p[1]for p in points)==0,'centered seed')
            require(sum(r[0]for r in res)==sum(r[1]for r in res)==0,'zero-sum residual')
            inner=(F(0),F(0))
            force_contract=(F(0),F(0))
            for p,r,f in zip(points,res,row['force']):
                inner=add(inner,mul(conjugate(p),r,d))
                force_contract=add(force_contract,mul(f,r,d))
            require(inner==(0,0),'complex Hermitian tangent identity')
            require(force_contract==(R,0),'linear force gain identity')
            require(R==sum(norm(f,d)for f in row['force'])-F(row['N']**2)/Q,'force residual norm identity')
            require(row['before']==independent_integer_value(points,d),'independent seed product')
            require(row['after']==independent_integer_value(row['new_points'],d),'independent updated product')
            if n<=7:
                sharp={2:F(4),3:F(27),4:F(256),5:F(3125),
                       6:F(6**15,5**10),7:F(7**21,6**15)}[n]
                require(row['after']<=sharp,'known small-n sharp ceiling')
            scaled=residual_data([scale(p,F(7,3))for p in points],d)
            require(scaled['normalized_residual']==Q*R,'scale invariant normalized residual')
            if row['stationary']:
                stationary.append((n,d))
                require(R==0 and row['before']==row['after'],'stationary branch')
            else:
                t=row['step'];g=row['log_gain_lower']
                require(sum(norm(p,d)for p in row['new_points'])==Q+t*t*R,'exact updated norm')
                require(t*row['C']<F(1,2),'collision-free step disk')
                require(log_lower(row['after']/row['before'])>=g,'independent exact logarithm lower comparison')
                plan=stopping_budget(points,F(1,10),d)
                require(g>=Q*R/plan['gain_denominator'],'global level-set gain estimate')
                require(plan['normalized_separation_floor']<=min(norm(sub(p,q),d)for p,q in combinations(points,2))/Q,'normalized separation floor')
                improved.append((n,d))
            rows.append({'n':n,'d':d,'stationary':row['stationary'],
                         'before':digest(row['before']),'after':digest(row['after']),
                         'log_gain_lower':str(row['log_gain_lower'])if n<=6 else digest(row['log_gain_lower'])})
            checks+=1
    # A saddle and a strict nonglobal local maximum both have zero residual.
    line=ascent_step([(-1,0),(0,0),(1,0)])
    require(line['stationary'] and line['before']==F(27,2),'collinear stationary counterexample')
    square=ascent_step([(0,0),(1,0),(0,1),(-1,0),(0,-1)])
    require(square['stationary'] and square['before']==F(5**10,4**6)<5**5,'nonglobal stationary center square')
    two=iterate([(-1,0),(1,0)],F(1,100),max_steps=0)
    require(two['status']=='stationarity_tolerance_met' and not two['rows'],'zero-update true tolerance')
    truncated=iterate([(0,0),(2,0),(0,1)],F(1,10**10),max_steps=0)
    require(truncated['status']=='update_budget_exhausted','budget scope is explicit')
    follow=iterate([(0,0),(2,0),(0,1)],F(1,10**10),max_steps=2)
    require(len(follow['rows'])==2 and all(r['after']>r['before']for r in follow['rows']),'two certified steps')
    report={'passed':True,'point_shapes_checked':checks,
            'strictly_improved_cases':[list(x)for x in improved],
            'stationary_seed_cases':[list(x)for x in stationary],
            'arithmetic':'integers, Fraction and exact logarithm lower series; no floating optimizer',
            'scope':'Finite identities and source seeds; general theorem is ascent_audit.md. Residual tolerance is not global objective tolerance.',
            'rows':rows}
    (Path(__file__).parent/'rational_ascent_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items()if k!='rows'},indent=2))


if __name__=='__main__':main()
