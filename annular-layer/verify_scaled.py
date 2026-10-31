"""Accept variable-smearing all-n tail and finite consequences."""
from fractions import Fraction as F
from math import isqrt
import sys
from acceptance import ROOT,require,write_report
from scaled_circumference import certificate,smearing_penalty_upper,pi_bounds
sys.path.insert(0,str(ROOT/'layer-circumference'))
from pair_integral import pair_integral_lower
from circumference import elementary_gap_upper


def run():
    h=F(109,100);c=F(3,2);cost=F(179,10000)
    raw=pair_integral_lower(c)
    require(raw>cost,'exact P(3/2)>179/10000')
    n0=10**14;k=F(51,100);tmax=F(1,3000);rootmin=10**7
    require(n0>=3000**4,'all-tail fourth-root direction')
    e=F(143,8000)
    require(elementary_gap_upper(n0)/n0+smearing_penalty_upper(h)<e,'all-tail decreasing corrected energy budget')
    require(e*(1+tmax*tmax)/2<F(19,200)**2,'test transport root bound')
    require(1+F(19,200)+h*h/8<F(5,4),'all-tail radial loss bounded by 5t/4')
    dist=2*pi_bounds()[1]*(k-F(1,4)+(k/4+F(5,16))*tmax)/(h*(1-k/rootmin))
    lower=(1-k/rootmin)*cost
    require(dist<c and lower>e,'one rational comparison covers every n>=10^14')
    examples=[]
    for n in (10**14,2*10**14,10**16,10**20):
        K=isqrt(n*51**2)//100
        fourth=isqrt(isqrt(n));fourth+=fourth**4<n
        result=certificate(n,K,elementary_gap_upper(n),F(1,fourth),h,c,cost)
        require(result['whole_class_excluded'],'generic finite certificate')
        examples.append(result)
    rejected=0
    for args in [(1,1,1,F(1,2),1,c,cost),(8,8,1,F(1,2),1,c,cost),
                 (8,2,-1,F(1,2),1,c,cost),(8,2,1,0,1,c,cost),
                 (8,2,1,F(1,2),0,c,cost),(8,2,1,F(1,2),1,c,0)]:
        try:certificate(*args)
        except ValueError:rejected+=1
    require(rejected==6,'invalid scaled inputs rejected')
    return {'passed':True,'smearing_multiplier':h,'pair_cutoff':c,
        'pair_cost_lower':cost,'quadrature_exceeds_claim':raw>cost,
        'all_n_threshold':n0,'all_n_strict_layer_coefficient':k,
        'tail_corrected_normalized_budget_upper':e,'tail_scaled_distance_upper':dist,
        'tail_pair_energy_per_n_lower':lower,'tail_strict_energy_margin':lower-e,
        'finite_examples':examples,'invalid_inputs_rejected':rejected,
        'scope':'For every global optimizer and every n>=10^14, K>(51/100)sqrt(n), using only the inherited effective elementary gap and author capacity/radius inputs. Examples do not replace the all-tail proof.'}


if __name__=='__main__':write_report('scaled_verification.json',run())
