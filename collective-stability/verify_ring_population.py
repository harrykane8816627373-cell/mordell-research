"""Exact population budgets, complete small integer compositions and constants."""
from fractions import Fraction as F
from pathlib import Path
import json
from ring_population import pi_interval,population_cost,minimum_cost,minimum_layers,population_certificate


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def compositions(n):
    if n==0:
        yield ()
    for m in range(1,n+1):
        for tail in compositions(n-m):
            yield (m,)+tail


def run():
    lo,hi=pi_interval()
    require(F(314159265358979323846,10**20)<lo<hi<F(314159265358979323847,10**20),
            'rational Machin enclosure')
    # Exact tangent arithmetic proves tan(4atan(1/5)-atan(1/239))=1.
    def tangent_add(a,b):return (a+b)/(1-a*b)
    t=tangent_add(F(1,5),F(1,5))
    t=tangent_add(t,t)
    require(tangent_add(t,F(-1,239))==1,'Machin tangent identity')
    count=0
    for n in range(2,17):
        true={}
        for pops in compositions(n):
            cost=population_cost(pops)
            k=len(pops)
            true[k]=min(true.get(k,cost),cost)
            count+=1
        costs=[]
        for k in range(1,n+1):
            require(minimum_cost(n,k)['cost']==true[k],'independent complete composition minimum')
            costs.append(true[k])
        require(costs==sorted(costs,reverse=True),'cost decreases on splitting layers')
        necessary=minimum_layers(n)['minimum_layers']
        require(necessary==next(k for k in range(1,n+1) if true[k]<hi**2*n*(n-1)),
                'binary search matches full K sweep')
    # Verify the singleton monotonicity derivative numerator exactly.
    derivatives=0
    for n in range(3,80):
        for k in range(1,(2*n)//3+1):
            for s in range(k):
                require(2*n-3*k+s>=0,'singleton derivative sign')
                require(F((n-s)**3,(k-s)**2)>=F(n**3,k*k),'continuous cost floor')
                derivatives+=1
    rows=[]
    for n in (2,8,9,30,40,80,100,1000,10000,10**6,10**8):
        cert=minimum_layers(n)
        if cert['last_excluded_layers']:
            require(cert['last_excluded_cost']['cost']>=cert['cost_cap_upper'],'last K impossible')
        require(minimum_cost(n,cert['minimum_layers'])['cost']<cert['cost_cap_upper'],'next K not excluded')
        rows.append(cert)
    require(population_certificate((80,20))['excluded_from_global_maxima'],'dominant hundred-point circle excluded')
    require(not population_certificate((2,8,14,20,26,30))['excluded_from_global_maxima'],
            'previous constructive population is not spuriously excluded')
    invalid=0
    for action in (lambda:minimum_cost(0,1),lambda:minimum_cost(10,11),
                   lambda:minimum_cost(100,20,max_cases=1),lambda:population_cost((True,2))):
        try:action()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid input accepted')
    return {'passed':True,'independent_ordered_compositions':count,
        'singleton_monotonicity_checks':derivatives,'invalid_inputs_rejected':invalid,
        'general_theorem':'Every global maximizer has K > n/(pi sqrt(n-1)) > sqrt(n)/pi, n>=2.',
        'pi_interval':[lo,hi],'rows':rows,
        'scope':'Exact necessary population budgets and general proof; finite composition checks validate the minimization code only.'}


if __name__=='__main__':
    report=run()
    Path(__file__).with_name('ring_population_verification.json').write_text(
        json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','pi_interval')},ensure_ascii=False,indent=2))
