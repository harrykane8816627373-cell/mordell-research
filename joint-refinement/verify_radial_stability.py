"""Independent identities, SOS and Bernstein checks for the fourth-mode budget."""
from fractions import Fraction as F
from math import comb
import json
import radial_stability as S
from verify import C
from verify_ellipse_dual import physical_witness

E=S.E


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def add(p,q):
    return [(p[i] if i<len(p) else F(0))+(q[i] if i<len(q) else F(0)) for i in range(max(len(p),len(q)))]


def multiply(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            out[i+j]+=a*b
    return out


def scale(p,a):
    return [a*v for v in p]


def identity_checks():
    # Q=120(t-3/4)(t-1)^2+132(t-191/176)^2+469/704.
    q=add(add(scale(multiply([F(-3,4),1],multiply([-1,1],[-1,1])),120),
        scale(multiply([F(-191,176),1],[F(-191,176),1]),132)),[F(469,704)])
    need(q==[F(529,8),F(27,2),F(-198),F(120)],'independent exact SOS expansion')
    need(multiply([F(-3,4),1],q)==S.COEFFICIENTS,'g has its certified nonnegative factorization on t>=3/4')
    need(56*7-8*F(1587,32)==F(-19,4),'normalization leaves the exact -19/4 offset')
    radial_eighth_upper=F(7,4)**4+7*F(3,4)**4
    need(radial_eighth_upper==F(371,32),'annular mass simplex gives an exact eighth-moment upper bound')
    singular_squared_lower=16+F(3,2)/radial_eighth_upper
    need(singular_squared_lower==F(5984,371) and singular_squared_lower>F(257,64)**2,
        'annular fourth-mode gap gives full real Hessian minimum strictly below -1/64')
    cases=[physical_witness(),[C(1),C(-1),C(0,1),C(0,-1),C(1,F(1,2)),
        C(-1,-F(1,2)),C(F(1,2)),C(-F(1,2))]]
    for points in cases:
        need(sum(points,C())==C() and sum(z.norm() for z in points)==7,'actual centered normalized identity input')
        t=[z.norm() for z in points]
        s2=sum((z**2 for z in points),C())
        s4=sum((z**4 for z in points),C())
        Q2=sum((r*z**2 for r,z in zip(t,points)),C())
        W2=8*Q2-13*s2
        action=[(32*r-36)*z**2-2*s2 for r,z in zip(t,points)]
        M44=64*sum(r**4 for r in t)-8*s4.norm()-sum(z.norm() for z in action)
        f=sum(-120*r**4+288*r**3-162*r*r for r in t)
        need(M44/8==f+4*s2.norm()-s4.norm()+2*(s2.conjugate()*W2).re,
            'full Schur diagonal is independently expanded with the nonstationary Ward residual retained')
        need(f==F(-19,4)-sum(S.evaluate(r) for r in t),'actual radial mass identity')
    return {'actual_mass_seven_identity_cases':len(cases),'nonstationary_residual_retained':True,
        'exact_SOS_coefficients_checked':True,'inner_particle_threshold':'3/4',
        'required_s2_squared_upper':'19/16','old_ellipse_implies_inner_particle':True,
        'annular_eighth_moment_upper':'371/32','annular_Hessian_min_strict_upper':'-1/64'}


def bernstein_checks():
    intervals=[(F(0),F(49,32)),(F(3,4),F(1)),(F(7,8),F(7,8)),(F(1,2),F(3,2)),
        (F(13,16),F(15,16))]
    queries=0
    for interval in intervals:
        beta,mass=S.bernstein(interval)
        lo,hi=interval
        for k in range(21):
            u=F(k,20)
            t=lo+(hi-lo)*u
            basis=[comb(4,j)*u**j*(1-u)**(4-j) for j in range(5)]
            need(sum(basis)==1 and all(v>=0 for v in basis),'Bernstein basis is a nonnegative partition of one')
            need(sum(b*v for b,v in zip(beta,basis))==S.evaluate(t),'direct polynomial agrees with independent Bernstein expansion')
            for eta in [F(-7,3),F(0),F(101,7)]:
                value=S.evaluate(t)-eta*t
                need(min(b-eta*m for b,m in zip(beta,mass))<=value,
                    'each mass-adjusted lower coefficient is a genuine continuous lower bound')
            queries+=1
    singleton=[(F(7,8),F(7,8))]*8
    bound=S.mass_lower(singleton)
    need(not bound['mass_infeasible'] and bound['lower']==8*S.evaluate(F(7,8)),
        'normalized equal-radius singleton budget is exact')
    broader=[(F(3,4),F(1))]*8
    bound=S.mass_lower(broader)
    need(bound['lower']<=8*S.evaluate(F(7,8)),'mass-dual lower bound retains an actual normalized radial witness')
    # Independently scan additional parameters; the selected finite-candidate lower bound dominates them.
    records=[S.bernstein(p) for p in broader]
    for k in range(-100,101):
        eta=F(k,7)
        value=7*eta+sum(min(b-eta*m for b,m in zip(beta,mass)) for beta,mass in records)
        need(value<=bound['lower'],'exact finite dual candidates dominate independent parameter scan')
    return {'independent_Bernstein_queries':queries,'additional_mass_parameters_checked':201,
        'floating_point_acceptance':False}


def annular_witness():
    base=[C(F(2,3)),C(-F(2,3)),C(0,F(5,6)),C(0,-F(5,6)),
        C(1),C(-1),C(0,F(7,6)),C(0,-F(7,6))]
    ints=[(935,0),(661,661),(0,935),(-661,661),(-935,0),(-661,-661),(0,-935),(661,-661)]
    target=[C(F(a,1000),F(b,1000)) for a,b in ints]
    direction=[a-b for a,b in zip(base,target)]
    numerator=sum((a.conjugate()*b for a,b in zip(direction,base)),C()).re
    denominator=sum(z.norm() for z in direction)
    q=[a-2*numerator/denominator*b for a,b in zip(base,direction)]
    norms=[z.norm() for z in q]
    s2=sum((z**2 for z in q),C())
    s3=sum((z**3 for z in q),C())
    s4=sum((z**4 for z in q),C())
    x,y=sum(t*t for t in norms),sum(t**3 for t in norms)
    A=196-24*x-s2.norm()
    B=312*x-128*y-1183-2*s3.norm()
    mixed=sum((t*t*z for t,z in zip(norms,q)),C())
    cross=40*mixed+s2.conjugate()*s3
    need(sum(q,C())==C() and sum(norms)==7 and q[0].im==0 and len(set(q))==8,
        'actual annular witness has exact center, mass and coordinate gauge')
    need(min(norms)>F(3,4) and s2 and s3 and s2.norm()+F(54,175)*s3.norm()<1,
        'nonempty physical annular region has both low moments nonzero and satisfies the old ellipse')
    need(A>0 and B>0 and A*B-2*cross.norm()>0,'old algebraic low-mode inequalities alone do not reject this actual witness')
    upper=sum(-120*t**4+288*t**3-162*t*t for t in norms)+4*s2.norm()-s4.norm()
    need(upper<0,'the full fourth-mode radial necessary inequality rejects it')
    E.setup(40)
    qbox=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in q]
    reason,evidence=S.exclusion(qbox)
    need(reason=='Hessian_inner_particle','exact interval implementation rejects the actual annular point')
    return {'actual_centered_mass_seven':True,'inside_nonzero_moment_ellipse':True,
        'all_squared_radii_above_three_quarters':True,'old_A_B_cross_inequalities_strict':True,
        'fourth_mode_radial_rejection':True,'witness_is_claimed_stationary':False}


if __name__=='__main__':
    print(json.dumps({'passed':True,'fourth_mode_identity':identity_checks(),
        'continuous_radial_bounds':bernstein_checks(),'actual_annular_witness':annular_witness(),
        'full_eight_point_exclusion_proved':False},indent=2))
