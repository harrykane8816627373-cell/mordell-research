"""Independent parabolic-domain and polynomial-integration acceptance."""
from fractions import Fraction as F
from pathlib import Path
import json
from overlap_cost import pair_cost_lower,pair_cost_squared_lower,certified_pi_upper

HERE=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    # Direct integrate (a-s)^2 on [-a,a], then the normalized y polynomial.
    inner=sum((F(v,degree+1)*(1-(-1)**(degree+1)) for degree,v in enumerate([1,-2,1])),F(0))
    outer=sum((F(v,2*j+1) for j,v in enumerate([1,-3,3,-1])),F(0))
    require(inner==F(8,3) and outer==F(16,35),'independent parabolic integral factors')
    require(2*inner*2*outer==F(512,105),'two ordered terms and two y halves')
    beta_integral=2*(F(1,3)-F(2,5)+F(1,7))
    require(beta_integral==F(16,105) and 4*beta_integral==F(64,105),'independent near-contact leading integral')
    domain=derivatives=0
    for d in (F(3,2),F(8,5),F(7,4),F(19,10),F(1913,1000)):
        eps=1-d/2;D=2*(d-1)
        for i in range(9):
            y2=D*eps*F(i,8)**2;a=eps-y2/D
            for j in range(-8,9):
                s=a*F(j,8)
                q1=(d/2-s)**2+y2;q2=(d/2+s)**2+y2
                upper1=d/2-s+y2/D;upper2=d/2+s+y2/D
                require(0<=upper1<=1 and 0<=upper2<=1 and upper1*upper1>=q1 and upper2*upper2>=q2,
                        'independent full parabolic domain containment')
                require(1-upper2==a-s,'asymmetric gap retained, not replaced by |s|')
                domain+=1
        direct=2*eps**7-7*(d-1)*eps**6
        require(direct==eps**6*(9-8*d)<0,'monotone squared pair bound')
        derivatives+=1
    rows=[]
    for c in (F(8,5),F(7,4),F(19,10),F(1913,1000)):
        new=pair_cost_lower(c);old=(1-c/2)**4/8
        require(new>old,'parabolic region improves old inscribed-circle bound')
        rows.append({'cutoff':c,'new_pair_cost_lower':new,'old_pair_cost_lower':old,'improvement_ratio_lower':new/old})
    require(pair_cost_lower(F(8,5))>F(3,500),'simple short-edge coefficient')
    invalid=0
    for call in (lambda:pair_cost_lower(1),lambda:pair_cost_lower(2),lambda:pair_cost_lower(F(8,5),0)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid pair cost accepted')
    return {'passed':True,'normalized_inner_integral':inner,'normalized_outer_integral':outer,
        'parabolic_domain_checks':domain,'monotonicity_checks':derivatives,'pi_upper':certified_pi_upper(),
        'near_contact_leading_constant_times_pi':F(64,105),'near_contact_power':'7/2',
        'pair_cost_examples':rows,'invalid_inputs_rejected':invalid,
        'scope':'Finite integration/domain checks; PARABOLIC_PAIR.md proves the general overlap lower bound for sqrt(2)<=distance<2.'}


if __name__=='__main__':
    result=run()
    (HERE/'overlap_cost_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,default=str,ensure_ascii=False,indent=2))
