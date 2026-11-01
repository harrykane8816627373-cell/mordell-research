"""Joint fourth-mode Hessian, old-ellipse and normalized radial budget."""
from fractions import Fraction as F
from math import comb
import refinement as R
import mass_dual as D

E=R.E
# g(t)=56t-1587/32-(-120t^4+288t^3-162t^2).
COEFFICIENTS=[F(-1587,32),F(56),F(162),F(-288),F(120)]


def evaluate(t):
    result=F(0)
    for c in reversed(COEFFICIENTS):
        result=result*t+c
    return result


def bernstein(interval):
    """Exact degree-four Bernstein coefficients and linear mass coefficients."""
    lo,hi=map(F,interval)
    if lo>hi:
        raise ValueError('ordered exact radial interval')
    width=hi-lo
    power=[sum(COEFFICIENTS[j]*comb(j,k)*lo**(j-k)*width**k for j in range(k,5)) for k in range(5)]
    beta=[sum(power[j]*F(comb(k,j),comb(4,j)) for j in range(k+1)) for k in range(5)]
    mass=[lo+F(k,4)*width for k in range(5)]
    return beta,mass


def mass_lower(intervals):
    if len(intervals)!=8:
        raise ValueError('eight radial intervals')
    intervals=[tuple(map(F,p)) for p in intervals]
    if sum(a for a,b in intervals)>7 or sum(b for a,b in intervals)<7:
        return {'mass_infeasible':True}
    records=[bernstein(p) for p in intervals]
    hulls=[]
    candidates={F(0)}
    for beta,mass in records:
        lines={}
        for b,m in zip(beta,mass):
            lines[-m]=max(lines.get(-m,-b),-b)
        hull=D.upper_hull(lines)
        hulls.append(hull)
        candidates.update(start for slope,intercept,start in hull[1:])
    def lower(eta):
        return 7*eta-sum(D.hull_value(hull,eta) for hull in hulls)
    chosen=max(candidates,key=lambda t:(lower(t),-t))
    return {'mass_infeasible':False,'lower':lower(chosen),'multiplier':chosen,'candidate_count':len(candidates)}


def exclusion(qbox):
    norms=[E.cnorm(q) for q in qbox]
    if all(4*t[0]>=3*E.S for t in norms):
        return 'Hessian_inner_particle',{'minimum_radial_lower_raw':min(t[0] for t in norms)}
    intervals=[(F(max(0,a),E.S),F(b,E.S)) for a,b in norms]
    bound=mass_lower(intervals)
    if bound['mass_infeasible']:
        return 'Hessian_radial_mass',{}
    s3,s4=E.ci(0),E.ci(0)
    for q in qbox:
        q2=E.cm(q,q)
        s3=E.ca(s3,E.cm(q2,q))
        s4=E.ca(s4,E.cm(q2,q2))
    joint=bound['lower']+F(E.cnorm(s4)[0],E.S)+F(216*E.cnorm(s3)[0],175*E.S)
    if joint>=F(-3,4):
        return 'Hessian_radial_budget',{'joint_lower':str(joint),'mass_multiplier':str(bound['multiplier']),
            'parameter_candidates':bound['candidate_count']}
    return None,None
