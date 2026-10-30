"""Exact local geometric proof constants and independent polynomial identities."""
from fractions import Fraction as F
from pathlib import Path
import json
from packing_constants import (validate,packing_less_than,CLOSE,QUERY,HOLE,INFLATED,
                              NEAR,GREEDY,KAPPA,PAIR_COST,HOLE_COST)

HERE=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    constants=validate();distances=arcs=0
    # Every rectangle endpoint bound used by the separate-convexity argument.
    for i in range(21):
        x=CLOSE+(NEAR-CLOSE)*F(i,20)
        for j in range(21):
            y=CLOSE+(NEAR-CLOSE)*F(j,20)
            d2=x*x+y*y-2*x*y*F(31,50)
            require(d2<=constants['seven_neighbor_corner_distance_squared_upper']<CLOSE*CLOSE,
                    'independent rational annulus distance bound')
            distances+=1
    beta=constants['cap_cosine_lower']
    for j in range(21):
        d=CLOSE+(NEAR-CLOSE)*F(j,20)
        current=(d*d+QUERY*QUERY-INFLATED*INFLATED)/(2*d*QUERY)
        require(current>=beta and 4*current*current>3,'inflated arc is shorter than pi/3')
        arcs+=1
    # h(1-u)=-log(1-u)-u+u^2/2 >=u^2 from positive terms.
    pair=0
    for j in range(1,20):
        u=F(j,40)
        h_lower=sum((u**k/k for k in range(1,30)),F(0))-u+u*u/2
        require(h_lower>=u*u,'one-disk logarithmic gap from independent positive series')
        pair+=1
    require((1-CLOSE/2)**4/8==PAIR_COST==32*KAPPA,'ordered pair coefficient')
    require(HOLE_COST/17>=KAPPA,'local holes compensate all quiet population')
    invalid=0
    for call in (lambda:packing_less_than(-1,17),lambda:packing_less_than(1,0)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid geometric bound accepted')
    return {'passed':True,'constants':constants,'rational_annulus_checks':distances,
        'inflated_cap_checks':arcs,'independent_log_gap_checks':pair,
        'invalid_inputs_rejected':invalid,
        'scope':'Finite exact constant acceptance. The general local alternative uses angular pigeonhole, connected cap lengths and packing, proved in LOCAL_ALTERNATIVE.md.'}


if __name__=='__main__':
    out=run()
    (HERE/'geometry_verification.json').write_text(json.dumps(out,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,default=str,ensure_ascii=False,indent=2))
