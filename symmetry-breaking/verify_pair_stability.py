"""Exact algebra, independent complex evaluation and whole-box exclusion."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from pair_stability import pair_slack,pair_box_upper

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'configuration-bounds'))
from rational_ascent import add,sub,scale,mul,conjugate,norm,divide


def require(ok,label):
    if not ok:
        raise ArithmeticError(label)


def run():
    identities=0
    # Values need not be stationary: this checks the exact algebraic identity.
    points=[(F(a,3),F(b,5)) for a,b in [(0,0),(1,1),(2,-3),(-5,2),(4,7)]]
    for n in range(2,16):
        lam=F(n-1,2)
        for z in points:
            for w in points:
                if z==w:
                    continue
                a,b,s=norm(z,1),norm(w,1),norm(sub(z,w),1)
                d=sub(z,w)
                first=divide((1/lam,F(0)),mul(d,d,1),1)
                second=scale(conjugate(mul(z,w,1)),F(1,n))
                analytic=(1-F(1,n)-a/n)*(1-F(1,n)-b/n)-norm(add(first,second),1)
                require(pair_slack(n,a,b,s)==lam*n*s*s*analytic,'independent complex identity')
                identities+=1
    # Known rational-coordinate stationary square, and center+square.
    square=[(F(1),F(0)),(F(0),F(1)),(F(-1),F(0)),(F(0),F(-1))]
    stable_checks=0
    for n,points,q in [(4,square,F(1)),(5,[(F(0),F(0))]+square,F(5,4))]:
        for i,z in enumerate(points):
            for w in points[i+1:]:
                require(pair_slack(n,q*norm(z,1),q*norm(w,1),q*norm(sub(z,w),1))>=0,
                        'known local maximum pair passes')
                stable_checks+=1
    box=((F(49,100),F(51,100)),(F(149,100),F(151,100)),(F(1,3),F(17,50)))
    upper=pair_box_upper(8,*box)
    require(upper<0,'entire thick radial-distance box excluded')
    require(pair_slack(8,F(1,2),F(3,2),F(1,3))==F(-925,1008),'exact central pair example')
    # Independent corners and internal points are coding checks, not a coverage proof.
    sample_checks=0
    for i in range(5):
        for j in range(5):
            for k in range(5):
                a=box[0][0]+F(i,4)*(box[0][1]-box[0][0])
                b=box[1][0]+F(j,4)*(box[1][1]-box[1][0])
                s=box[2][0]+F(k,4)*(box[2][1]-box[2][0])
                require(pair_slack(8,a,b,s)<=upper,'box upper covers evaluation')
                sample_checks+=1
    bad=0
    for action in (lambda:pair_slack(1,1,1,1),lambda:pair_slack(8,-1,1,1),
                   lambda:pair_box_upper(8,(1,0),(1,2),(1,2)),
                   lambda:pair_box_upper(8,(0,1),(1,2),(F(1,10),1))):
        try:action()
        except ValueError:bad+=1
        else:raise ArithmeticError('invalid input accepted')
    return {'passed':True,'independent_complex_identities':identities,
        'known_local_maximum_pairs':stable_checks,'box_code_checks':sample_checks,
        'invalid_inputs_rejected':bad,'excluded_n8_box':box,'box_slack_upper':upper,
        'exact_example_slack':F(-925,1008),
        'scope':'The analytic projection and interval-product proofs cover arbitrary local maxima and the stated entire box. Finite evaluations only cross-check implementation; no full eight-point exclusion is claimed.'}


if __name__=='__main__':
    report=run()
    Path(__file__).with_name('pair_stability_verification.json').write_text(
        json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,ensure_ascii=False,indent=2))
