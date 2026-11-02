"""Exact constants and algebra for the seven-point/target-value/Schwarz radius theorem."""
from fractions import Fraction as Q
import json
from verify import C,need
import target_radius as R


def add(a,b):
    n=max(len(a),len(b));return [(a[i] if i<len(a) else Q())+(b[i] if i<len(b) else Q()) for i in range(n)]


def multiply(a,b):
    out=[Q()]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def check():
    # 27/256-x(1-x)^3=(x-1/4)^2((x-5/4)^2+1/8).
    left=[Q(27,256),Q(-1),Q(3),Q(-3),Q(1)]
    right=multiply(multiply([Q(-1,4),Q(1)],[Q(-1,4),Q(1)]),
        add(multiply([Q(-5,4),Q(1)],[Q(-5,4),Q(1)]),[Q(1,8)]))
    need(left==right,'all coefficients of the global quartic maximum certificate')
    coefficient=Q(2**6*7**14,3**15)
    bound=coefficient*Q(27,256)**7
    need(bound==Q(3**6*7**14,2**50)<Q(4,9),'exact uniform constant-term modulus bound below 2/3')
    for t in [Q(0),Q(1),Q(119,80),Q(49,32),Q(49,8)]:
        x=Q(8,49)*t;mass=7-Q(8,7)*t
        lhs=Q(64,49)**7*t**7*mass**21/Q(7**7*6**15)
        rhs=coefficient*(x*(1-x)**3)**7
        need(lhs==rhs<=bound,'independent exact scale reduction including both endpoints')
    # Check the disk automorphism derivative factor and its sharp examples independently.
    examples=[C(),C(Q(1,2)),C(Q(1,3),Q(1,4)),C(Q(-1,4),Q(1,5)),C(Q(2,3)-Q(1,10**6))]
    for a in examples:
        u=(C(1)-a)/(C(1)-a.conjugate())
        f1=(a+u)/(C(1)+a.conjugate()*u)
        derivative=u*(1-a.norm())/(C(1)+a.conjugate()*u)**2
        expected=(C(1)-a).norm()/(1-a.norm())
        need(f1==C(1) and derivative==C(expected),'sharp boundary-derivative disk-automorphism identity')
        need(expected>Q(1,5),'each checked constant term below 2/3 gives the required strict derivative gap')
    qlower=(21+Q(1,5))/4
    radius=Q(7,8)*(7-qlower)
    need(qlower==Q(53,10) and radius==Q(119,80),'exact final mass and squared-radius constants')
    need(Q(49,32)-radius==Q(7,160) and radius/Q(49,32)==Q(34,35),'comparison with the inherited unquantified radius')
    R.E.setup(40);settings=R.E.model('global_maximum');strong=R.model(settings)
    need(strong['radius2']==R.E.rat(119,80)[1] and strong['radius']**2>=strong['radius2']*R.E.S,
        'integer implementation rounds the new squared radius and radius outward')
    try:R.model(R.E.model('user_physical_system'))
    except ValueError:pass
    else:raise ArithmeticError('global target-radius theorem escaped its scope')
    return {'passed':True,'quartic_nonnegative_factorization_coefficients_checked':5,
        'uniform_constant_term_modulus_squared_upper':str(bound),'strict_modulus_upper':'2/3',
        'sharp_boundary_automorphism_examples':len(examples),'strict_boundary_derivative_lower':'1/5',
        'seven_subset_centered_mass_strict_lower':'53/10','squared_radius_strict_upper':'119/80',
        'squared_radius_gain_over_49_over_32':'7/160','new_to_old_squared_radius_ratio':'34/35',
        'seven_point_sharp_bound_is_inherited_external_input':True,
        'seven_point_full_classification_rebuilt_here':False,'analytic_Schwarz_bridge_requires_written_proof':True,
        'global_scope_guard_checked':True,'floating_point_acceptance':False,
        'full_eight_point_exclusion_proved':False}


if __name__=='__main__':print(json.dumps(check(),indent=2))
