"""n100 consequences of the stronger n-K layer mass in the shifted dual."""
from fractions import Fraction as F
import sys
from acceptance import ROOT,require,write_report
sys.path.insert(0,str(ROOT/'joint-variational'))
from shared_budget import radial_packet,add_centered_holes,exact_certificate
from affine_pair import fixed_line
sys.path.insert(0,str(ROOT/'stationary-screening'))
from stationary_moments import low_mode_coefficient
from separated_dual import strengthened_certificate
import json


def run():
    E=F(413923829261,500000000000);line=fixed_line()
    prior=json.loads((ROOT/'gaussian-angular/joint_budget_verification.json').read_text(encoding='utf-8'))
    require(all(F(row['global_gap_upper'])==E for row in prior['conditional_contrast_bounds']),'actual inherited global gap upper')
    packet=radial_packet(100,F(5,16),F(11,10),[0])
    base_cert=strengthened_certificate(packet,4,line)
    hole_cert=strengthened_certificate(add_centered_holes(packet,[F(1,10)]),4,line,force_radius_lower=F(2))
    base=base_cert['shared_gap_lower'];hole=hole_cert['shared_gap_lower']
    c2=low_mode_coefficient(100,2,F(11881,5000))
    c3=low_mode_coefficient(100,3,F(11881,5000))
    require(base+F(143,250)>E and hole+F(149,500)>E,'new angular-energy thresholds')
    require(c2*F(47,20)**2>F(143,250),'n100 |P2|<47/20')
    require(c2*F(17,10)**2>F(149,500),'center-hole n100 |P2|<17/10')
    contrasts=[]
    for label,radial,A,bound in [('K<=4',base,F(0),F(3,8)),
            ('K<=4 and low-mode energy>=11/20',base,F(11,20),F(73,1000)),
            ('K<=4 and all |z_i|^2>=8/99',hole,F(0),F(27,100))]:
        total=radial+A+F(2000,487)*bound*bound
        require(total>E,'all-scale Gaussian contrast exclusion '+label)
        contrasts.append({'hypothesis':label,'strict_contrast_upper':bound,
            'joint_lower_at_excluded_threshold':total,'strict_margin':total-E})
    return {'passed':True,'n':100,'maximum_layers':4,'global_gap_upper':E,
        'pair_line':{k:line[k] for k in ('alpha','beta','left','right')},
        'new_base_certificate':base_cert,'new_centered_hole_certificate':hole_cert,
        'new_base_shared_lower':base,'new_centered_hole_shared_lower':hole,
        'angular_energy_strict_upper':'143/250','centered_hole_angular_energy_strict_upper':'149/500',
        'original_z_p2_strict_upper':'47/20','centered_hole_original_z_p2_strict_upper':'17/10',
        'n100_original_z_c2':c2,'n100_original_z_c3':c3,
        'all_scale_gaussian_contrast_bounds':contrasts,
        'scope':'Conditional necessary inequalities for actual global optimizers, using h=1 and full prior radial Gram. Does not prove five layers or a sharp n100 value. Gaussian contrasts use actual default disk means and are orthogonal to radial and k=2,3 tests.'}


if __name__=='__main__':write_report('joint_verification.json',run())
