"""Combinatorial ledger, actual separated inputs, and exact linear corollaries."""
from fractions import Fraction as F
from pathlib import Path
from math import prod,isqrt
import json
import sys
from packing_constants import KAPPA,PAIR_COST,HOLE_COST
from linear_gap import linear_gap_lower,ledger_lower,matching_plan,gamma_upper_interval
from subset_moment import log_interval

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'stationary-screening/imported-angular'))
import angular_certificate as old


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def hex_points(m):
    return [(2*a+b,b) for a in range(-m,m+1) for b in range(-m,m+1)
            if max(abs(a),abs(b),abs(a+b))<=m]


def run():
    counts=0
    for n in range(2,151):
        # Any boundary count B and actual quiet count satisfying Q>=n-B-32M.
        for B in (0,n//3,n):
            for m in range(n//2+1):
                q=max(0,n-B-32*m)
                if q>n-2*m:continue
                result=ledger_lower(n,m,q)
                require(result['energy_lower']>=KAPPA*max(0,n-B),'matching plus quiet-hole accounting')
                counts+=1
    examples=[]
    for n in (122,144,400,10000,10**6,10**8,10**12):
        lower=linear_gap_lower(n)
        require(lower>=0,'effective linear gap nonnegative')
        if n in (144,400,10000,10**6,10**8,10**12):require(lower>0,'positive large finite correction')
        examples.append({'n':n,'linear_gap_lower':lower})
    plans=[]
    for m in (2,3,4,5):
        points=hex_points(m)
        plan=matching_plan(points,imaginary_unit_squared=3)
        require(not plan['matching'] and plan['quiet_population']>0 and plan['selected_quiet_centers'],
                'real separated quiet configurations exercise hole branch')
        if m==2:
            n=len(points);Q=sum(x*x+3*y*y for x,y in points);N=n*(n-1)//2
            product=F(n,Q)**N*prod((x-a)**2+3*(y-b)**2
                for i,(x,y) in enumerate(points) for a,b in points[i+1:])
            require(old.deficit_interval(n,product)[0]>plan['conditional_energy_lower'],
                    'physical conditional energy below independently computed true deficit')
        plans.append(plan)
    stored=json.loads((HERE.parent/'symmetry-breaking/integer_witnesses.json').read_text(encoding='utf-8'))
    for n in (10,30,100):
        pts=next(r['integer_points'] for r in stored['rows'] if r['n']==n)
        plan=matching_plan(pts)
        require(plan['matching'] and plan['conditional_energy_lower']>0,'actual prior input exercises close-pair branch')
        plans.append(plan)
    scale=10**12;q=F(2,143)
    square_root_lower=F(isqrt(q.numerator*scale*scale//q.denominator),scale)
    old_B_lower=(1+square_root_lower)**2
    old_discount_upper=144*log_interval(2)[1]**2/old_B_lower**144
    require(old_discount_upper<F(1,10**10)<linear_gap_lower(144),'new finite correction strictly exceeds the old Newton correction')
    gamma=gamma_upper_interval()
    require(gamma[1]<F('0.403426350115382569900758939272'),'outward conditional gamma endpoint')
    require(KAPPA==F(1,16777216),'effective asymptotic gap')
    invalid=0
    for call in (lambda:linear_gap_lower(True),lambda:ledger_lower(8,5,0),
                 lambda:matching_plan([(0,0)]*4),lambda:matching_plan(hex_points(2),max_points=10),
                 lambda:matching_plan(hex_points(2),imaginary_unit_squared=0)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid deficit input accepted')
    return {'passed':True,'independent_integer_ledgers':counts,'effective_examples':examples,
        'conditional_geometric_plans':plans,'gamma_upper_outward_if_existing_expansion_used':gamma,
        'n144_old_Newton_discount_strict_upper':'1/10000000000',
        'n144_new_linear_discount_lower':linear_gap_lower(144),
        'invalid_inputs_rejected':invalid,'universal_correction':'2^-24*(n-11sqrt(n))_+',
        'scope':'Finite all-n theorem is analytic and uses inherited global separation/radius. Plans certify only their geometric and energy hypotheses, not optimum status. Gamma corollary retains the previously cited expansion input.'}


if __name__=='__main__':
    out=run()
    (HERE/'defect_ledger_verification.json').write_text(json.dumps(out,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='conditional_geometric_plans'},default=str,ensure_ascii=False,indent=2))
