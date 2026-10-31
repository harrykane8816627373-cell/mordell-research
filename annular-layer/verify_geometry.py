"""Independent nontrivial ring and convex-mass checks for the m-1 lemma."""
from fractions import Fraction as F
import sys
from acceptance import ROOT,require,write_report
sys.path.insert(0,str(ROOT/'joint-variational'))
from affine_pair import sqrt_bounds,pi_bounds,pair_integral_upper
sys.path.insert(0,str(ROOT/'layer-circumference'))
from pair_integral import pair_integral_lower
from affine_pair import fixed_line
from separated_dual import strengthened_certificate
from shared_budget import radial_packet


def hinge(d):return max(F(0),2-d)**2


def run():
    pl,pu=pi_bounds()
    # A 24-gon is a real circle, not an abstract independently variable edge list.
    # cos(pi/12)=sqrt(2+sqrt(3))/2; adjacent chord at R=23/4 is exact.
    s3=sqrt_bounds(3)
    c_lo=sqrt_bounds(2+s3[0])[0]/2;c_hi=sqrt_bounds(2+s3[1])[1]/2
    R=F(23,4)
    chord=sqrt_bounds(2*R*R*(1-c_lo))[1]
    require(chord*chord>=2*R*R*(1-c_hi) and chord*chord>2,'actual chord upper in lens domain')
    actual_lower=24*pair_integral_lower(chord)
    new_lower=23*pair_integral_lower(pu/2)
    new_upper=23*pair_integral_upper(pl/2)
    old_upper=22*pair_integral_upper(23*pl/44)
    require(actual_lower>new_upper,'actual regular-ring cost dominates the new convex artificial mass')
    require(new_lower>old_upper,'new surrogate is strictly stronger than the old m-2 ring surrogate')
    # Irregular rational points on one circle; test the geometry with a different
    # convex decreasing kernel, so this check is independent of P quadrature.
    circle=[(F(1),F(0)),(F(3,5),F(4,5)),(F(0),F(1)),
            (F(-4,5),F(3,5)),(F(-1),F(0)),(F(-3,5),F(-4,5)),
            (F(0),F(-1)),(F(4,5),F(-3,5))]
    checks=0
    for m in range(2,9):
        for radius in [F(1,4),F(1,2),F(3,4),F(1)]:
            points=circle[:m]
            edges=[(0,1)] if m==2 else [(j,(j+1)%m) for j in range(m)]
            cost=sum((hinge(sqrt_bounds(radius*radius*((points[i][0]-points[j][0])**2+
                  (points[i][1]-points[j][1])**2))[1]) for i,j in edges),F(0))
            lower=(m-1)*hinge(2*pu*radius/(m-1))
            require(cost>=lower,'irregular circle convex mass inequality')
            checks+=1
    # Capacity gives the restricted domain exactly, without a large-n limit.
    for m in range(2,101):
        radius=F(m-1,4)
        require(2*pl*radius/(m-1)==pl/2,'capacity-to-restricted-domain identity')
    line=fixed_line();alpha,beta=line['alpha'],line['beta']
    r0_upper=alpha/(2*pl*beta)
    require(2*r0_upper<sqrt_bounds(2)[0],'two unpaid singleton centers contradict global separation')
    require(alpha<sqrt_bounds(2)[0]*(pl-1)*beta,'two-point layer is paid at full mass')
    rejected=0;packet=radial_packet(100,F(5,16),F(11,10),[0])
    for action in [lambda:strengthened_certificate(packet,4,dict(line,alpha=10*beta)),
                   lambda:strengthened_certificate(packet,4,line,force_radius_lower=F(1,100))]:
        try:action()
        except ValueError:rejected+=1
    require(rejected==2,'invalid separation-line or inner-radius hypotheses rejected')
    return {'passed':True,'regular_polygon_population':24,'regular_polygon_radius_force':R,
        'actual_adjacent_chord_upper':chord,'actual_cycle_cost_lower':actual_lower,
        'new_m_minus_one_surrogate_lower':new_lower,'new_m_minus_one_surrogate_upper':new_upper,
        'old_m_minus_two_surrogate_upper':old_upper,
        'strict_surrogate_gain_lower':new_lower-old_upper,
        'irregular_rational_circle_checks':checks,'exact_capacity_identities':99,
        'unpaid_singleton_force_radius_upper':r0_upper,
        'separation_adapter_invalid_inputs_rejected':rejected,
        'scope':'Independent ring geometry and exact convex-mass checks. The 24-gon is used to test the subfunctional and is not claimed an actual normalized global optimizer.'}


if __name__=='__main__':write_report('geometry_verification.json',run())
