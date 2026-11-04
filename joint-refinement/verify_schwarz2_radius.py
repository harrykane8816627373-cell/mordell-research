"""Exact centering and constants for the second-order Schwarz radius theorem."""
from fractions import Fraction as Q
import json
from verify import C,need
import schwarz2_radius as R


def multiply(a,b):
    out=[C()]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=out[i+j]+x*y
    return out


def check():
    # Independent exact products, including b=0: no division by the constant term.
    offsets=[C(1),C(Q(-1,2),Q(1,3)),C(Q(1,5),Q(-2,7)),
        C(Q(-3,4),Q(1,6)),C(Q(2,3),Q(-1,5)),C(Q(-1,7),Q(2,9))]
    offsets.append(-sum(offsets,C()))
    need(sum(offsets,C())==C(),'the seven offsets are exactly centered')
    products=0
    for b in [C(),C(1),C(Q(2,3),Q(-1,4)),C(Q(-7,5),Q(2,9)),C(Q(1,1000))]:
        c=-b*Q(8,7)
        h=[C(1)]
        for y in offsets:h=multiply(h,[c,y])
        need(len(h)==8 and h[0]==c**7 and h[1]==C(),
            'all factors expanded: H(0)=(-8b/7)^7 and H prime at zero is exactly zero')
        if b==C():need(all(x==C() for x in h[:7]),'the zero-center case is a seventh-order monomial')
        products+=1
    bound=Q(3**6*7**14,2**50)
    need(bound<Q(4,9),'inherited exact uniform constant-term bound is strictly below 4/9')
    examples=[C(),C(Q(1,2)),C(Q(1,3),Q(1,4)),C(Q(-1,4),Q(1,5)),C(Q(2,3)-Q(1,10**6))]
    for a in examples:
        # These rational disk maps have f'(0)=0 and realize the second-order
        # boundary derivative factor; the universal statement is proved in prose.
        u=(C(1)-a)/(C(1)-a.conjugate())
        numerator=[a,C(),u]
        denominator=[C(1),C(),a.conjugate()*u]
        need(numerator[1]*denominator[0]-numerator[0]*denominator[1]==C(),
            'the sharp rational example has derivative zero at the origin')
        f1=(a+u)/(C(1)+a.conjugate()*u)
        derivative=2*u*(1-a.norm())/(C(1)+a.conjugate()*u)**2
        expected=2*(C(1)-a).norm()/(1-a.norm())
        need(f1==C(1) and derivative==C(expected),
            'second-order sharp boundary derivative identity')
        need(expected>Q(2,5),'the checked constant terms give the strict second-order gap')
    qlower=(21+Q(2,5))/4
    radius=Q(7,8)*(7-qlower)
    need(qlower==Q(107,20) and radius==Q(231,160),'exact mass and squared-radius constants')
    need(Q(119,80)-radius==Q(7,160) and radius/Q(119,80)==Q(33,34),
        'strict improvement over the accepted first-order Schwarz radius')
    need(Q(49,32)-radius==Q(7,80) and radius/Q(49,32)==Q(33,35),
        'comparison with the earlier unquantified radius')
    R.E.setup(40)
    settings=R.E.model('global_maximum');strong=R.model(settings)
    need(strong['radius2']==R.E.rat(231,160)[1]
        and strong['radius']**2>=strong['radius2']*R.E.S,
        'squared radius and radius are rounded outward on the integer mesh')
    need(settings['radius2']!=strong['radius2'],'the model copy does not overwrite the old global model')
    try:R.model(R.E.model('user_physical_system'))
    except ValueError:pass
    else:raise ArithmeticError('second-order Schwarz theorem escaped its global scope')
    return {'passed':True,'centered_polynomial_products_checked':products,
        'zero_constant_term_case_checked':True,'conditional_polynomial_derivative_at_zero':'0',
        'strict_constant_term_modulus_upper':'2/3','sharp_second_order_boundary_examples':len(examples),
        'strict_boundary_derivative_lower':'2/5','seven_subset_centered_mass_strict_lower':'107/20',
        'squared_radius_strict_upper':'231/160','gain_over_first_order_radius':'7/160',
        'analytic_second_order_Schwarz_bridge_requires_written_proof':True,
        'seven_point_sharp_bound_is_inherited_external_input':True,
        'seven_point_full_classification_rebuilt_here':False,'global_scope_guard_checked':True,
        'floating_point_acceptance':False,'full_eight_point_exclusion_proved':False}


if __name__=='__main__':print(json.dumps(check(),indent=2))
