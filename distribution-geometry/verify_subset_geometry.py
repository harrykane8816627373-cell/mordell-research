"""Exact acceptance for thin-strip geometry and a full 200-point witness."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import sys

from subset_geometry import *

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'symmetry-breaking'))
from pair_stability import pair_slack


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def run():
    span_checks=0
    # Independent enumeration of which elementary gaps enter each skip.
    for m in range(2,25):
        gaps=[F(j*j+1,10) for j in range(1,m)]
        x=[F(0)]
        for gap in gaps:x.append(x[-1]+gap)
        L=x[-1]
        for k in range(1,m):
            direct=sum(x[i+k]-x[i] for i in range(m-k))
            require(direct<=min(k,m-k)*L,'all skip spans fit exact multiplicity bound')
            counts=[sum(i<=j<i+k for i in range(m-k)) for j in range(m-1)]
            require(max(counts)==min(k,m-k),'maximum elementary-gap multiplicity')
            require(direct==sum(a*b for a,b in zip(counts,gaps)),'independent span expansion')
            span_checks+=1
    require(odd_skip_sum(2,1)==1 and odd_skip_sum(3,1)==8 and odd_skip_sum(4,1)==28,
            'independent small line reciprocal sums')
    convexity_checks=0
    for W in (F(0),F(1,100),F(1,10),F(1)):
        for j in range(1,25):
            t=W+F(j,13)
            # F''=(6t^2-2W^2)/(t^2+W^2)^3.
            require(6*t*t-2*W*W>0,'convexity on permitted longitudinal differences')
            convexity_checks+=1
    phase_checks=0
    for dx in (F(1,5),F(1,2),F(1),F(3)):
        for dy in (F(-1,100),F(0),F(1,100)):
            W=F(1,100)
            d2=dx*dx+dy*dy
            delta=F(4,199)
            if d2<delta:continue
            eta=1-2*W*W/delta
            # Real part of -B = Re(1/(dx+i dy)^2).
            actual=(dx*dx-dy*dy)/(d2*d2)
            require(actual>=eta/d2,'direct complex reciprocal real-part lower')
            require(1/d2>=1/(dx*dx+W*W),'transverse rectangle distance bound')
            phase_checks+=1
    # Unequal horizontal gaps and nonzero transverse coordinates independently
    # test the Jensen/ordering bound before the operator inequality is imposed.
    geometric_checks=0
    for m in (4,7,12,20):
        W=F(1,100)
        xs=[F(i,5)+F(i*i,500) for i in range(m)]
        ys=[(-1)**i*W/2 for i in range(m)]
        L=xs[-1]-xs[0]
        direct=F(0)
        for i in range(m):
            for j in range(i+1,m):
                if (j-i)%2:
                    direct+=1/((xs[j]-xs[i])**2+(ys[j]-ys[i])**2)
        require(direct>=odd_skip_sum(m,L*L,W*W),'direct unequal-gap reciprocal sum dominates new lower')
        geometric_checks+=1

    n,m,L,W=200,20,F(27,10),F(1,100)
    local=strip_capacity(n,m,L,W)
    global_test=strip_capacity(n,m,L,W,True)
    microscopic=microscopic_strip_capacity(m,(n-1)*L*L,(n-1)*W*W)
    require(microscopic['squared_stability_slack']==global_test['squared_stability_slack']/(n-1)**2,
            'exact microscopic scaling removes all n dependence')
    require(local['excluded'] and global_test['excluded'],'entire rectangle excludes local and global maxima')
    require(local['cross_operator_lower']>1018 and local['cross_operator_upper_squared']==995**2,
            'simple exact strict strip contradiction')
    witness=explicit_two_hundred_points()
    p=witness['points'];f=witness['normalization_squared']
    require(f==F(3724000,3724007)<1,'exact physical normalization')
    require(witness['grid_pair_norm_sum']==4571 and len(witness['grid_swaps'])==21,
            'deterministic finite grid selection replay')
    require(sum(x for x,y in p)==sum(y for x,y in p)==0,'full centered witness')
    require(f*sum(x*x+y*y for x,y in p)==n,'full normalized moment')
    require(f*(max(x for x,y in p[:20])-min(x for x,y in p[:20]))**2<L*L,
            'selected twenty points lie strictly within rectangle length')
    require(all(y==0 for x,y in p[:20]),'selected twenty points are inside positive-width strip')
    radii=[f*(x*x+y*y) for x,y in p]
    min_distance=None;min_slack=None;pairs=0
    for i,j in combinations(range(n),2):
        s=f*((p[i][0]-p[j][0])**2+(p[i][1]-p[j][1])**2)
        t=pair_slack(n,radii[i],radii[j],s)
        min_distance=s if min_distance is None else min(min_distance,s)
        min_slack=t if min_slack is None else min(min_slack,t)
        pairs+=1
    require(min_distance>F(4,n-1) and min_slack>0,'full witness passes all old global distances and local pair slacks')
    require(max(radii)<F(2*(n-1),n),'full witness passes old global radius bound')
    subset_Q=f*sum(x*x+y*y for x,y in p[:20])
    complement_Q=f*sum(x*x+y*y for x,y in p[20:])
    require(subset_Q>F(m*(m-1),n-1) and complement_Q>F((n-m)*(n-m-1),n-1),
            'selected subset AND complement pass old strict internal-moment constraints')
    require(L*L>F(4*(m-1)**2,n-1),'old collinear pair packing allows stated count and length')
    stronger_line_examples=[]
    for M in (16,20,40,100):
        new=line_length_squared_lower(200,M)
        old=F(4*(M-1)**2,199)
        require(new>old,'recorded new finite line bound improves old global pair packing')
        stronger_line_examples.append({'population':M,'new_local_span_squared_lower':new,
                                       'old_global_pair_span_squared_lower':old})
    invalid=[lambda:strip_capacity(8,1,1),lambda:strip_capacity(8,9,1),
             lambda:strip_capacity(8,3,0),lambda:strip_capacity(8,3,1,1),
             lambda:strip_capacity(200,20,1,max_terms=1),lambda:odd_skip_sum(4,1,-1)]
    for call in invalid:
        try:call()
        except ValueError:pass
        else:raise ArithmeticError('invalid or incomplete strip certificate accepted')
    return {'passed':True,'arithmetic':'integers/Fraction; no floating-point acceptance',
        'independent_skip_span_checks':span_checks,'convexity_checks':convexity_checks,
        'complex_phase_real_checks':phase_checks,'direct_unequal_gap_examples':geometric_checks,
        'old_full_pair_slacks_checked':pairs,'rejected_invalid_cases':len(invalid),
        'entire_strip_local_certificate':local,'entire_strip_global_certificate':global_test,
        'universal_scaled_strip_certificate':microscopic,
        'nonempty_full_witness':{k:v for k,v in witness.items() if k!='points'}|{
            'point_count':n,'selected_strip_population':m,'minimum_squared_distance':min_distance,
            'maximum_squared_radius':max(radii),'minimum_old_pair_slack':min_slack,
            'selected_subset_internal_moment':subset_Q,'complement_internal_moment':complement_Q},
        'finite_line_improvements':stronger_line_examples,
        'scope':'Universal strip capacity follows from the original norm bound plus ordered-gap convexity. The full normalized witness passes the listed prior conditions, not every earlier global theorem. The new criterion is necessary geometry, not an energy or optimum certificate.'}


if __name__=='__main__':
    report=run()
    ROOT.joinpath('subset_geometry_verification.json').write_text(json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('nonempty_full_witness','finite_line_improvements')},default=str,ensure_ascii=False,indent=2))
