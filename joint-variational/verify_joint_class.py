"""A concrete n100 angular-energy/four-layer exclusion beyond separate budgets."""
from fractions import Fraction as F
from affine_pair import fixed_line,sqrt_bounds,pi_bounds,pair_integral_upper,ROOT
from shared_budget import radial_packet,add_centered_holes,exact_certificate,interval_certificate,geometry_constants
from acceptance import require,write_report
import json,sys
sys.path.insert(0,str(ROOT/'stationary-screening'))
from stationary_moments import low_mode_coefficient


def run():
    n,K,a,T=100,4,F(5,16),F(11,10)
    packet=radial_packet(n,a,T,[0]);line=fixed_line()
    result=exact_certificate(packet,K,line)
    E=F(413923829261,500000000000);A=F(151,200);remaining=E-A
    prior=json.loads((ROOT/'overlap-budget/graph_budget_verification.json').read_text(encoding='utf-8'))
    require(F(prior['short_edge_capacity_examples'][1]['remaining_gap'])==E,'actual inherited global n100 budget')
    require(result['shared_gap_lower']>E-A,'new shared certificate excludes the whole class')
    require(A<E,'orthogonal energy alone does not exclude')
    # The strengthened separate relaxation even uses the matched-moment Schur residual.
    L,V=geometry_constants(n,K,a);h=result['residual_norm_squared']
    d_lower=2*pi_bounds()[0]*(V+sqrt_bounds(h*remaining,10**12)[0]/4)/L
    d_test=F(8999,5000)
    require(d_lower>d_test,'separate-relaxation mean distance direction')
    pair_upper=L*pair_integral_upper(d_test)
    require(pair_upper<remaining,'separate radial and pair caps remain simultaneously feasible')
    S=F(11881,5000);c2=low_mode_coefficient(100,2,S);c3=low_mode_coefficient(100,3,S)
    require(c2*F(27,10)**2>A,'concrete original-z second-power-sum threshold')
    require(c2*F(27,10)**2<E,'previous angular ellipse permits this threshold')
    box_packet=radial_packet(100,a,T,[0,0])
    whole_box=interval_certificate(box_packet,4,line,[(0,0),(F(-1,100),F(1,100))],orthogonal_energy_lower=A)
    require(whole_box['shared_gap_lower']>E,'joint exclusion survives the entire declared radial moment box')
    holes=add_centered_holes(packet,[F(1,10)])
    hole_result=exact_certificate(holes,4,line)
    hole_A=F(51,100)
    require(hole_A+hole_result['shared_gap_lower']>E,'angular plus central-hole plus layer shared exclusion')
    require(hole_A+hole_result['projection_energy']<E,'angular plus hole Gram alone passes')
    require(c2*F(9,4)**2>hole_A and c2*F(9,4)**2+hole_result['projection_energy']<E,'concrete hole p2 threshold and old projection comparison')
    return {'passed':True,'n':100,'hypothesized_maximum_layers':4,'a':a,'cover_radius_q':T,
            'certified_pair_line':line,'joint_radial_pair_certificate':result,
            'orthogonal_angular_energy_class_lower':A,'global_gap_upper':E,
            'joint_total_lower':A+result['shared_gap_lower'],
            'strict_exclusion_margin':A+result['shared_gap_lower']-E,
            'separate_residual_budget':remaining,'separate_mean_distance_lower':d_lower,
            'separate_pair_cost_upper':pair_upper,
            'original_z_p2_threshold':'27/10','n100_p2_coefficient':c2,'n100_p3_coefficient':c3,
            'conditional_ellipse':'K<=4 implies C2*|P2|^2+C3*|P3|^2<151/200; hence |P2|<27/10',
            'whole_radial_measurement_box':whole_box,
            'centered_energy_hole_radius_q':'1/10','centered_hole_shared_certificate':hole_result,
            'centered_hole_angular_class_lower':hole_A,
            'centered_hole_joint_total_lower':hole_A+hole_result['shared_gap_lower'],
            'centered_hole_conditional_ellipse':'K<=4 and no smeared mass in q-radius 1/10 imply C2*|P2|^2+C3*|P3|^2<51/100; hence |P2|<9/4',
            'sufficient_original_z_hole_squared_radius':'8/99',
            'scope':'Universal conditional n100 exclusion. The comparison proves that the specified separate energy relaxation permits the data; it does not claim a physical optimizer or that every previous necessary condition permits the class.'}


if __name__=='__main__':write_report('joint_class_verification.json',run())
