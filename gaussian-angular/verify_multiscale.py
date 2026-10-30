"""Actual multiscale gain over the full old k=4 pair of complex tests, and a whole family."""
from fractions import Fraction as F
import json
from gaussian_energy import (standard_probe,measurement,probe_gram,add,sub,scale,
    square,root_upper,lipschitz_upper)
from legacy_projection import old_mode,packet_certificate,assembled_packet,dual_lower
from acceptance import require,write_report,nine_points,HERE


def run():
    points=nine_points();S=F(5,4);old=old_mode(points,4,S)
    probes=[standard_probe(F(1,24),F(3,4)),standard_probe()]
    values=[measurement(points,p) for p in probes]
    single=[packet_certificate(old,[p],[v]) for p,v in zip(probes,values)]
    joint=packet_certificate(old,probes,values)
    require(joint['energy_lower']>F(1288,10000),'fixed center joint lower')
    require(joint['energy_lower']>F(1059,1000)*old['old_projection'],'strict center improvement exceeds 5.9 percent')
    require(joint['energy_lower']>max(x['energy_lower'] for x in single),'two scales stronger than either accepted single scale')
    original=json.loads((HERE/'imported-review/gaussian_bridge_verification.json').read_text(encoding='utf-8'))
    require(joint['energy_lower']<F(original['independent_actual_total_gap_interval'][0]),'accepted projection below independently computed true product gap')
    # An exact normalized one-parameter family: a^2+2b^2=1 identically.
    eps=F(1,10**7);interval=(F(1,4)-eps,F(1,4)+eps)
    def a(s):return (1-2*s*s)/(1+2*s*s)
    def b(s):return 2*s/(1+2*s*s)
    abox=(a(interval[1]),a(interval[0]));bbox=(b(interval[0]),b(interval[1]))
    require(interval[0]>0 and 2*interval[1]**2<1,'monotone rational parameterization branches')
    qmax=max(abox[1]**2,2*bbox[1]**2);margin=S-qmax-F(1,9)
    require(margin>0 and margin*margin>4*qmax/9,'entire family smear inside old/new common support')
    a4=(abox[0]**4,abox[1]**4);b4=(bbox[0]**4,bbox[1]**4)
    a6=(abox[0]**6,abox[1]**6);b6=(bbox[0]**6,bbox[1]**6)
    p4=sub(scale(a4,4),scale(b4,16));Q4=sub(scale(a6,4),scale(b6,32))
    bubble=sub(scale(p4,S-F(5,18)),Q4)
    da=max(abs(v-F(7,9)) for v in abox);db=max(abs(v-F(4,9)) for v in bbox)
    displacement_sum=4*da+4*root_upper(2)*db
    family_values=[]
    for probe,v in zip(probes,values):
        error=lipschitz_upper(probe)*displacement_sum
        family_values.append((v[0]-error,v[1]+error))
    H,_=assembled_packet(old,probes,values)
    m=[p4,(F(0),F(0)),bubble,(F(0),F(0))]+family_values
    family_lower=dual_lower(H,m,joint['fixed_dual_vector'])
    old_upper=sum((square(v)[1]/h for v,h in zip(m[:4],old['diagonal'])),F(0))
    require(old_upper<F(61,500) and family_lower>F(16,125),'entire family improves the displayed class upper bound')
    require(family_lower>F(21,20)*old_upper,'entire family projection gain exceeds five percent')
    return {'passed':True,'common_support_squared_q':S,'old_full_k4_projection':old,
            'single_scale_energy_lowers':[x['energy_lower'] for x in single],
            'joint_certificate':joint,'center_improvement_factor_lower':joint['energy_lower']/old['old_projection'],
            'family_parameter_interval':interval,'family_a_interval':abox,'family_b_interval':bbox,
            'family_old_projection_upper':old_upper,'family_joint_projection_lower':family_lower,
            'family_fixed_dual_vector':joint['fixed_dual_vector'],'family_measurement_intervals':m,
            'scope':'A reflection/quarter-turn-preserving one-parameter normalized family, not an independent box in all coordinates, not asserted stationary or optimal. Comparison is against the specified complete old k=4 holomorphic/first-mixed projection, not every prior criterion.'}


if __name__=='__main__':write_report('multiscale_verification.json',run())
