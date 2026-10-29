"""Exact acceptance of the explicit finite-n layer corollary."""
from fractions import Fraction as F
from pathlib import Path
import json
from effective_layers import sqrt_interval,effective_layers


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    budget=F(1,100)+F(77,600)+F(35,40000)
    require(budget==F(16705,120000) and budget<F(7,50),'effective polar budget for all n>=10000')
    require(F(7,5)**2<2 and F(8,3)**2>7,'coarse certified roots')
    coefficient_upper=(F(5,4)+F(8,30))/(2*F(7,5))
    require(coefficient_upper==F(13,24)<F(11,20),'uniform strict finite error coefficient')
    root_checks=0
    for num in range(201):
        for den in (1,3,7):
            q=F(num,den)
            lo,hi=sqrt_interval(q,1000)
            require(lo*lo<=q<=hi*hi and 0<=hi-lo<=F(1,1000),'outward rational square root')
            root_checks+=1
    rows=[]
    for n in (10000,100000,10**6,10**8,10**12):
        row=effective_layers(n)
        lo,hi=row['sqrt_interval']
        require(lo*lo<=n<=hi*hi,'sqrt source enclosed')
        require(row['fourth_root_upper']**4>=n,'fourth root direction')
        require(row['minimum_layers']>row['strict_layer_lower_endpoint'],'integer strict endpoint')
        rows.append(row)
    require(rows[0]['minimum_layers']==45,'all 10000-point global optimizers have at least45 radii')
    require(rows[2]['minimum_layers']==483,'million-point minimum count')
    require(rows[3]['minimum_layers']==4946,'hundred-million strict integer boundary')
    require(rows[4]['minimum_layers']==499451,'trillion strict integer boundary')
    bad=0
    for action in (lambda:effective_layers(9999),lambda:effective_layers(True),
                   lambda:sqrt_interval(-1),lambda:sqrt_interval(2,0)):
        try:action()
        except ValueError:bad+=1
        else:raise ArithmeticError('invalid finite corollary input accepted')
    return {'passed':True,'root_enclosures_checked':root_checks,'invalid_inputs_rejected':bad,
        'old_effective_deficit_rate_upper_at_10000':budget,
        'error_coefficient_upper':coefficient_upper,
        'theorem':'Every global optimizer has K>sqrt(n)/2-(11/20)*n**(1/4) for all n>=10000.',
        'rows':rows,
        'scope':'Arithmetic verification of the closed finite corollary; universal bridges are proved in EFFECTIVE_LAYERS.md and the transport/capacity documents.'}


if __name__=='__main__':
    report=run()
    Path(__file__).with_name('effective_layers_verification.json').write_text(
        json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,ensure_ascii=False,indent=2))
