"""Accept the universal tail certificate, graph counts and asymptotic bridge."""
from fractions import Fraction as F
from math import isqrt
from pair_integral import fixed_certificates,FIXED_COSTS,certified_pi_upper,ROOT
from circumference import layer_exclusion,elementary_gap_upper
from reporting import require,write_report
import json
import sys
sys.path.insert(0,str(ROOT/'stationary-screening'))
sys.path.insert(0,str(ROOT/'local-defect'))
from stationary_moments import low_mode_coefficient
from radial_hierarchy import eight_radial_coefficient


def run():
    fixed_certificates()
    n0=10**16;k=F(5007,10000);tmax=F(1,10000);rootmin=10**8
    emax=F(10001,10**6)
    require(elementary_gap_upper(n0)/n0<emax,'all-n decreasing elementary normalized budget')
    require(emax*(1+tmax*tmax)/2<F(3,40)**2,'Dirichlet test error less than 3t/40')
    require(1+F(3,40)+F(1,8)==F(6,5),'radial-test loss at most 6t/5')
    dist=2*certified_pi_upper()*(k-F(1,4)+(k/4+F(3,10))*tmax)/(1-2*k/rootmin)
    cost=FIXED_COSTS[F(3151,2000)]
    require(dist<F(3151,2000),'all-n adjacent distance cutoff')
    lower=(1-2*k/rootmin)*cost
    require(lower>emax,'all-n pair energy contradicts elementary budget')
    examples=[]
    for n in (10**16,2*10**16,10**18,10**20,10**24):
        K=isqrt(n*5007**2)//10000
        fourth=isqrt(isqrt(n));fourth+=fourth**4<n
        row=layer_exclusion(n,K,elementary_gap_upper(n),F(1,fourth))
        require(row['entire_at_most_K_class_excluded'],'finite generic accepting example')
        examples.append({'n':n,'excluded_maximum_layers':K,'guaranteed_minimum_layers':K+1,
                         'generic_certificate':row})
    # Ordinary finite graph statistics use the new accurate integral, not the asymptotic input.
    gap8=F(523,5000);gap100=F(413923829261,500000000000);edgecost=FIXED_COSTS[F(8,5)]
    inherited=json.loads((ROOT/'overlap-budget/graph_budget_verification.json').read_text(encoding='utf-8'))
    require(F(inherited['short_edge_capacity_examples'][1]['remaining_gap'])==gap100,'inherited actual n100 gap')
    caps=[int(gap8/edgecost),int(gap100/edgecost)]
    require(caps==[12,101],'improved global short-edge populations')
    cr=eight_radial_coefficient();c2=low_mode_coefficient(8,2,F(151,50))
    base=cr*(F(101,10)-F(64,7))**2+c2*F(4,5)**2
    require(base<gap8<base+3*edgecost,'weaker eighth-point moment threshold plus three edges excluded')
    # Asymptotic result with the separately cited, inherited external input.
    require(FIXED_COSTS[F(44,25)]>F(1,800),'asymptotic budget separation')
    asym_low=F(1,4)+F(22,25)/certified_pi_upper()
    effective_asym_low=F(1,4)+F(3151,4000)/certified_pi_upper()
    require(asym_low>F(53,100) and effective_asym_low>k,'pi bounds on two asymptotic constants')
    return {'passed':True,'finite_all_n_threshold':n0,'finite_layer_coefficient':k,
            'tail_t_upper':tmax,'tail_normalized_gap_upper':emax,
            'tail_adjacent_distance_upper':dist,'tail_pair_energy_per_n_lower':lower,
            'finite_examples':examples,'force_short_edge_cutoff':'8/5',
            'n8_global_short_edge_cap':caps[0],'n100_global_short_edge_cap':caps[1],
            'n8_conditional_fourth_moment_lower':'101/10','n8_conditional_p2_modulus_lower':'4/5',
            'n8_conditional_short_edge_cap':2,'n8_joint_three_edge_lower':base+3*edgecost,
            'elementary_asymptotic_bound':'1/4+3151/(4000*pi)',
            'elementary_asymptotic_rational_lower':effective_asym_low,
            'external_input_asymptotic_bound':'1/4+22/(25*pi)',
            'external_input_asymptotic_rational_lower':asym_low,
            'external_input':'The previously cited linear energy theorem plus triangular trial bound limsup(delta_n/n)<=1/800; no effective threshold is inferred for this input.',
            'scope':'Finite n>=10^16 conclusion uses only the explicit elementary gap bound, author stability/radius, and exact smearing. The >0.53 asymptotic conclusion retains its external energy input.'}


if __name__=='__main__':write_report('layers_verification.json',run())
