"""Exact identities and constants for the old-ellipse Schur radius bound."""
from fractions import Fraction as Q
import json
from verify import C,need
from verify_schwarz2_radius import multiply
import ellipse_schur_radius as R


def check():
    upper=Q(231,160);threshold=Q(23,16);a=Q(329,500)
    x=Q(8,49)*upper
    need(x==Q(33,140)<Q(1,4),'the old radius lies on the increasing side of x(1-x)^3')
    # Its derivative is (1-x)^2(1-4x), verified coefficient by coefficient.
    inner=[C(),C(1),C(-3),C(3),C(-1)]
    derivative=[(i+1)*inner[i+1] for i in range(4)]
    factored=multiply(multiply([C(1),C(-1)],[C(1),C(-1)]),[C(1),C(-4)])
    need(derivative==factored,
        'monotonicity derivative factorization coefficients')
    alpha2=(Q(64,49)*upper)**7*(7-Q(8,7)*upper)**21/Q(7**7*6**15)
    need(alpha2<a*a,'exact improved constant-term bound below (329/500)^2')
    b=a/(1-a*a)*(Q(49,128)/threshold+Q(7,16))
    need(b==Q(10651375,13041828)<1,'uniform Schur constant-term upper below one')
    gap=(2+(1-b)/(1+b))*(1-a)/(1+a)
    need(gap==Q(8511842889,19641665287)>Q(3,7),
        'the Schur boundary gap strictly exceeds the high-radius derivative upper')
    need(gap-Q(3,7)==Q(657904362,137491657009),'independent rational difference')
    need(7-Q(32,7)*threshold==Q(3,7),'derivative upper at squared radius 23/16')
    need(upper-threshold==Q(1,160),'improvement over the universal second-order radius')
    offsets=[C(1),C(Q(-1,2),Q(1,3)),C(Q(1,5),Q(-2,7)),
        C(Q(-3,4),Q(1,6)),C(Q(2,3),Q(-1,5)),C(Q(-1,7),Q(2,9))]
    offsets.append(-sum(offsets,C()))
    products=0
    for point in [C(Q(1,2)),C(Q(2,3),Q(-1,4)),C(Q(-7,5),Q(2,9))]:
        points=[y-point/Q(7) for y in offsets]+[point]
        s2=sum((z**2 for z in points),C())
        p2=sum((y**2 for y in offsets),C())
        need(sum(points,C())==C() and p2==s2-Q(8,7)*point**2,
            's2 is the global power sum and the seven-point centered power sum is s2-8b squared/7')
        c=-point*Q(8,7);h=[C(1)]
        for y in offsets:h=multiply(h,[c,y])
        need(h[2]==-c**5*p2/2,'independent exact quadratic coefficient of H')
        need(sum(h,C())!=C(),'the algebra check has a nonzero normalization at t=1')
        alpha=h[0]/sum(h,C());f2=h[2]/sum(h,C())
        if alpha.norm()!=1:
            phase=(C(1)-alpha.conjugate())/(C(1)-alpha)
            beta=phase*f2/(1-alpha.norm())
            need(beta.norm()==alpha.norm()*p2.norm()/(4*c.norm()**2*(1-alpha.norm())**2),
                'exact squared modulus of the second Schur constant term')
        products+=1
    examples=[(C(),C()),(C(Q(1,3),Q(1,4)),C(Q(1,2))),
        (C(Q(-1,4),Q(1,5)),C(Q(1,4),Q(-1,3)))]
    for alpha,beta in examples:
        u=(C(1)-alpha)/(C(1)-alpha.conjugate())
        v=(C(1)-beta)/(C(1)-beta.conjugate())
        h1=(beta+v)/(C(1)+beta.conjugate()*v)
        hp1=v*(1-beta.norm())/(C(1)+beta.conjugate()*v)**2
        fp1=u*(1-alpha.norm())/(C(1)+alpha.conjugate()*u)**2*(2+hp1)
        expected=(C(1)-alpha).norm()/(1-alpha.norm())*(2+(C(1)-beta).norm()/(1-beta.norm()))
        need(h1==C(1) and fp1==C(expected),'exact two-stage Schur boundary-derivative identity')
    R.E.setup(40);settings=R.E.model('global_maximum')
    strong=R.model(settings,old_ellipse=True)
    need(strong['radius2']==R.E.rat(23,16)[1]
        and strong['radius']**2>=strong['radius2']*R.E.S,
        'the new radius is rounded outward')
    for settings_arg,kwargs in [(settings,{}),(settings,{'old_ellipse':False}),
        (R.E.model('user_physical_system'),{'old_ellipse':True})]:
        try:R.model(settings_arg,**kwargs)
        except ValueError:pass
        else:raise ArithmeticError('ellipse-Schur radius escaped one of its two required hypotheses')
    return {'passed':True,'centered_power_sum_products_checked':products,
        'two_stage_boundary_examples':len(examples),'constant_term_modulus_strict_upper':'329/500',
        'second_Schur_constant_term_modulus_upper':str(b),'boundary_derivative_strict_lower':str(gap),
        'high_radius_derivative_upper':'3/7','squared_radius_strict_upper':'23/16',
        'gain_over_universal_second_order_radius':'1/160',
        'requires_global_maximizer_and_old_moment_ellipse':True,
        'old_ellipse_hypothesis_must_be_explicit_in_model':True,'scope_guards_checked':3,
        's2_is_the_global_complex_power_sum':True,
        'analytic_Schur_bridge_requires_written_proof':True,
        'seven_point_sharp_bound_is_inherited_external_input':True,
        'seven_point_full_classification_rebuilt_here':False,'floating_point_acceptance':False,
        'full_eight_point_exclusion_proved':False}


if __name__=='__main__':print(json.dumps(check(),indent=2))
