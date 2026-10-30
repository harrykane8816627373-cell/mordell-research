"""Exact interpolation, exponential/tail and effective covering certificates."""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import json
from microscopic_covering import (exp_interval,sqrt_interval,parameters,TABLE,
                                  covering_certificate,disk_population)


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def cadd(z,w):return z[0]+w[0],z[1]+w[1]
def csub(z,w):return z[0]-w[0],z[1]-w[1]
def cmul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def cdiv(z,w):
    d=w[0]*w[0]+w[1]*w[1]
    return (z[0]*w[0]+z[1]*w[1])/d,(z[1]*w[0]-z[0]*w[1])/d


def run():
    interpolation=0
    for n in range(2,9):
        nodes=[(F(j-2),F((j*j)%5-2,3)) for j in range(n)]
        for u in [(F(1,3),F(2,5)),(F(-2,7),F(3,8))]:
            def poly(z):
                ans,power=(F(0),F(0)),(F(1),F(0))
                for j in range(n):
                    ans=cadd(ans,cmul((F(j+1),F((-1)**j,2)),power))
                    power=cmul(power,z)
                return ans
            total=(F(0),F(0))
            for i,node in enumerate(nodes):
                cardinal=(F(1),F(0))
                for j,other in enumerate(nodes):
                    if i!=j:cardinal=cmul(cardinal,cdiv(csub(u,other),csub(node,other)))
                total=cadd(total,cmul(poly(node),cardinal))
            require(total==poly(u),'independent exact polynomial interpolation')
            interpolation+=1
    poisson=0
    for n in (4,8,16,32,64):
        for fraction in (F(1,4),F(1,3),F(1,2)):
            s=n*fraction
            lowexp,highexp=exp_interval(s,40)
            partial=sum((s**k/F(factorial(k)) for k in range(n)),F(0))
            tau_upper=1-partial/highexp
            q=F(5,6)
            require(0<=1-partial/lowexp<=tau_upper<=q**n,'independent Poisson tail enclosure')
            poisson+=1
    rows=[]
    monotonic=0
    for bulk,(q,start) in TABLE.items():
        data=parameters(bulk,q,start)
        require(data['combined_error_upper']<F(3,4),'strict covering contradiction margin')
        for m in range(start//2,start//2+20):
            require(q*(2*m+3)<=2*m+1,'tail envelope decreases for all later m')
            monotonic+=1
        rows.append(covering_certificate(start,bulk))
    populations=[disk_population(10**6,(0,0),F(1,10)),
                 disk_population(10**6,(F(1,2),F(1,4)),F(1,5)),
                 disk_population(10000,(0,0),F(1,2))]
    require(populations[0]['lower']==246 and populations[0]['upper']==10200,'million-point microscopic population example')
    invalid=0
    for call in (lambda:covering_certificate(63),lambda:parameters(2,F(5,6),64),
                 lambda:parameters(1,F(1,2),64),lambda:parameters(1,F(5,6),4),
                 lambda:disk_population(100,(1,1),1),lambda:exp_interval(-1)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid covering input accepted')
    return {'passed':True,'independent_interpolation_identities':interpolation,
        'independent_poisson_tail_checks':poisson,'tail_monotonicity_checks':monotonic,
        'invalid_inputs_rejected':invalid,'fixed_bulk_certificates':rows,
        'population_examples':populations,
        'theorem':'For every global optimizer, n>=64 and every |u|<=1, min_i|z_i-u|<6/sqrt(n-1).',
        'scope':'Finite explicit microscopic bulk covering and coarse local population bounds; no boundary covering, optimal constants, exact optimizer or crystallization claimed.'}


if __name__=='__main__':
    report=run()
    Path(__file__).with_name('microscopic_covering_verification.json').write_text(
        json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,ensure_ascii=False,indent=2))
