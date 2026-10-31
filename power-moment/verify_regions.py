"""Whole moment boxes; a separate local-rigidity witness; explicit incomplete cover."""
from fractions import Fraction as F
from itertools import combinations
import json
from acceptance import HERE,require,write_report
from complex_interval import C,I
from power_coordinates import objective_box
from witnesses import normalized


def abs2_upper(v):return max(map(abs,v.re))**2+max(map(abs,v.im))**2
def abs2_lower(v):
    def dist(a):return F(0) if a[0]<=0<=a[1] else min(map(abs,a))
    return dist(v.re)**2+dist(v.im)**2


def run():
    source=json.loads((HERE/'region_inputs.json').read_text(encoding='utf-8'))
    require(source['full_remaining_region_covered'] is False,'nonexhaustive seed list cannot be promoted to a global cover')
    target=4**28*7**7;results=[]
    for row in source['regions']:
        raw=[C(*v) for v in row['raw_points']];center=[C(*v) for v in row['power_center']];radius=F(row['radius'])
        centered,ratio,powers,product=normalized(raw,28)
        require(center[:2]==[C(8),C(0)],'centered moment coordinates')
        for k in range(2,9):
            require(center[k].re-radius<=powers[k].re[0]<=powers[k].re[1]<=center[k].re+radius and
                    center[k].im-radius<=powers[k].im[0]<=powers[k].im[1]<=center[k].im+radius,'actual normalized point witness lies in the entire box')
        for k in (2,3):
            box=I((center[k].re-radius,center[k].re+radius),(center[k].im-radius,center[k].im+radius))
            require(abs2_lower(box)>0,'entire box belongs to the nonzero p2,p3 branch')
        cert=objective_box(center,radius);upper=cert['objective_modulus_upper']
        require(product<=upper<target,'whole moment box excluded at the exact candidate threshold')
        results.append({'id':row['id'],'radius':radius,'objective_upper_over_candidate':upper/target,
            'strict_relative_margin':1-upper/target,'exact_witness_product_over_candidate':product/target,
            'preconditioned_row_error_sums':cert['preconditioned_row_error_sums']})
    # A physically normalized nonstationary input in the rigidity box passes
    # every inherited strict subset inequality and the displayed low ellipse.
    raw=[C(*v) for v in source['near_candidate_raw_points']]
    centered,ratio,powers,product=normalized(raw,28)
    eps=F(1,10)
    for k in range(2,9):
        target_value=C(896) if k==7 else C(0)
        require(abs2_upper(powers[k]-target_value)<=eps**2,'nonempty full complex-radius local moment neighborhood')
    require(abs2_lower(powers[2])>0 and abs2_lower(powers[3])>0,'nonzero low moments inside the local region')
    subset_checks=0;minimum_slack=None
    for m in range(2,8):
        for indices in combinations(range(8),m):
            values=[centered[i] for i in indices];mean=sum(values,C(0))/m
            internal=ratio*sum(((z-mean).norm() for z in values),F(0))
            slack=internal-F(m*(m-1),2)
            require(slack>0,'actual normalized input passes inherited strict subset condition')
            minimum_slack=slack if minimum_slack is None else min(minimum_slack,slack)
            subset_checks+=1
    ellipse=F(2,7)**2*abs2_upper(powers[2])+F(27,100)*F(2,7)**3*abs2_upper(powers[3])
    require(ellipse<F(64,49) and product<target,'nearby nonzero-moment input passes old ellipse and is genuinely below the candidate')
    return {'passed':True,'n':8,'force_scale_candidate_product':target,
        'entire_independent_moment_boxes_excluded':len(results),'boxes':results,
        'local_neighborhood_complex_power_radius':'1/10',
        'nonempty_local_nonzero_p2_p3_witness':True,
        'nearby_witness_all_strict_subsets_checked':subset_checks,
        'nearby_witness_minimum_subset_slack':minimum_slack,
        'nearby_witness_old_ellipse_upper':ellipse,
        'full_eight_point_exclusion_proved':False,
        'scope':'Eleven literal whole-box objective exclusions and one nonempty local stationary-rigidity neighborhood, with all rotations/relabelings. No claim that this union covers all normalized stationary moment vectors. Some distant seed configurations also fail older subset/separation tests; no superiority over all old constraints is asserted.'}


if __name__=='__main__':write_report('regions_verification.json',run())
