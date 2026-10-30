"""Validated parameter families and independent integer vertex-cover ledgers."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from parameter_certificate import validate,FAMILY,family_certificate,gap_lower,gamma_upper_interval

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'local-defect'))
from linear_gap import linear_gap_lower


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    families={name:family_certificate(name) for name in FAMILY}
    ledgers=0
    for cert in families.values():
        A,G,k=[cert[k] for k in ('endpoint_cap','quiet_cap','linear_coefficient')]
        for n in range(2,101):
            for boundary in (0,n//3,n):
                for edges in range(n+1):
                    quiet=max(0,n-boundary-A*edges)
                    holes=(quiet+G-1)//G
                    energy=cert['pair_cost_used']*edges+cert['hole_cost']*holes
                    require(energy>=k*max(0,n-boundary),'independent all-integer cover/quiet ledger')
                    ledgers+=1
    rows=[]
    for n in (122,144,400,10000,10**6,10**8,10**12):
        old=linear_gap_lower(n)
        require(gap_lower(n,'simple32')==32*old and gap_lower(n,'strong37')==37*old,'exact finite improvement factor')
        rows.append({'n':n,'strong_gap_lower':gap_lower(n),'old_gap_lower':old})
    base=list(FAMILY['strong37']);invalid=0
    calls=[]
    for pos,bad in ((0,F(7,5)),(1,F(1)),(2,F(0)),(3,1),(4,1),(5,F(1,1000))):
        args=base.copy();args[pos]=bad
        calls.append(lambda args=args:validate(*args))
    calls.extend((lambda:validate(*base,boundary=1),lambda:validate(*base,cosine_lower=F(99,100)),
                  lambda:gap_lower(1),lambda:family_certificate('unlisted')))
    for call in calls:
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid parameter certificate accepted')
    return {'passed':True,'fixed_parameter_certificates':families,'integer_ledger_checks':ledgers,
        'finite_gap_examples':rows,'conditional_gamma_upper_interval':gamma_upper_interval(),
        'invalid_inputs_rejected':invalid,
        'scope':'All-n sufficiency for the supplied parameter families; no optimization over all parameter choices or point configurations is claimed.'}


if __name__=='__main__':
    result=run()
    (HERE/'parameters_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='fixed_parameter_certificates'},default=str,ensure_ascii=False,indent=2))
