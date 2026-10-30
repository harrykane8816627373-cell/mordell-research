"""Finite n100 local angular-density contrasts coupled to layer and moment budgets."""
from fractions import Fraction as F
from pathlib import Path
import sys,json
from gaussian_energy import standard_probe,probe_gram
from acceptance import require,write_report

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'joint-variational'))
sys.path.insert(0,str(ROOT/'stationary-screening'))
from affine_pair import fixed_line
from shared_budget import radial_packet,add_centered_holes,exact_certificate
from stationary_moments import low_mode_coefficient


def run():
    n=100;probe=standard_probe();norm=probe_gram(probe,probe);h=F(487,2000)
    require(0<norm[0]<=norm[1]<h,'certified fixed Gaussian norm cap')
    E=F(413923829261,500000000000)
    prior=json.loads((ROOT/'joint-variational/joint_class_verification.json').read_text(encoding='utf-8'))
    require(F(prior['global_gap_upper'])==E,'inherited actual n100 global gap')
    packet=radial_packet(100,F(5,16),F(11,10),[0]);line=fixed_line()
    radial=exact_certificate(packet,4,line)['shared_gap_lower']
    with_hole=exact_certificate(add_centered_holes(packet,[F(1,10)]),4,line)['shared_gap_lower']
    cases=[('four_layers',F(0),F(857,2000),radial,F(151,200)),
           ('four_layers_and_low_modes',F(7,10),F(57,500),radial,F(151,200)),
           ('four_layers_and_centered_hole',F(0),F(351,1000),with_hole,F(51,100))]
    results=[]
    for name,A,contrast,geometry,old_cap in cases:
        angular=A+contrast*contrast/h;total=geometry+angular
        require(angular<old_cap<E,'displayed previous angular cap permits the compared data')
        require(total>E,'new local-density class excluded by precise shared budget')
        results.append({'name':name,'old_modes_energy_lower':A,'excluded_contrast_modulus_lower':contrast,
            'Gaussian_energy_lower':contrast*contrast/h,'joint_gap_lower':total,
            'global_gap_upper':E,'strict_margin':total-E,'displayed_previous_angular_cap':old_cap})
    S=F(11881,5000)
    return {'passed':True,'n':100,'probe':probe,'Dirichlet_norm_interval':norm,
            'accepted_norm_upper':h,'q_coordinates':'q=sqrt(99/200)*z','actual_smear_squared':'1/100',
            'legacy_p2_coefficient':low_mode_coefficient(100,2,S),
            'legacy_p3_coefficient':low_mode_coefficient(100,3,S),
            'radial_layer_shared_lower':radial,'centered_hole_layer_shared_lower':with_hole,
            'conditional_contrast_bounds':results,
            'universal_normalized_contrast_bound':'sup_(gamma>0,theta) |M_(gamma,theta)|/n <= sqrt((487/2000)*E)/n',
            'scope':'Necessary bounds for global n100 optimizers, valid at every self-similar scale and rotation by exact norm invariance. K<=4 is required; the second case also requires C2|P2|^2+C3|P3|^2>=7/10; the third requires every original |z_i|^2>=8/99. Does not prove n100 needs five layers.'}


if __name__=='__main__':write_report('joint_budget_verification.json',run())
