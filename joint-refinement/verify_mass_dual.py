"""Independent exact checks for normalization ranges and grid-dual envelopes."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import gzip
import json
import tempfile
import mass_dual as M
import refinement as R
from check import check

E=R.E
ROOT=Path(__file__).resolve().parent.parent


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def vertices(intervals):
    active=[i for i,(lo,hi) in enumerate(intervals) if lo!=hi]
    for free in active or [None]:
        other=[i for i in range(8) if i!=free]
        choices=[(intervals[i][0],) if intervals[i][0]==intervals[i][1] else intervals[i] for i in other]
        for values in product(*choices):
            x=[0]*8
            for i,value in zip(other,values):
                x[i]=value
            if free is not None:
                x[free]=-sum(values)
            if sum(x)==0 and all(lo<=v<=hi for v,(lo,hi) in zip(x,intervals)):
                yield x


def independent_minimum(intervals):
    choices=[('fixed',) if lo==hi else ('lo','hi','free') for lo,hi in intervals]
    best=None
    for labels in product(*choices):
        free=[i for i,label in enumerate(labels) if label=='free']
        x=[F(intervals[i][1] if label=='hi' else intervals[i][0]) for i,label in enumerate(labels)]
        if free:
            value=-sum(x[i] for i in range(8) if i not in free)/len(free)
            for i in free:
                x[i]=value
        if sum(x)==0 and all(lo<=v<=hi for v,(lo,hi) in zip(x,intervals)):
            norm=sum(v*v for v in x)
            best=norm if best is None else min(best,norm)
    return best


def projection_and_lines():
    E.setup(40)
    cases=[[(1,2)]+[(-1,0)]*7,
        [(10,100),(-1,0),(-100,-5)]+[(0,0)]*5,
        [(0,2),(-3,-1),(-1,2),(0,4)]+[(0,0)]*4,
        [(-2,-1),(1,3),(-1,1)]+[(0,0)]*5,
        [(-1,-1),(1,1)]+[(0,0)]*6,
        [(-1,1)]+[(0,0)]*7,
        [(1,2)]*8]
    records=0
    for index,small in enumerate(cases):
        intervals=[(lo*E.S,hi*E.S) for lo,hi in small]
        actual=M.projection(intervals)
        expected=independent_minimum(intervals)
        need((actual is None)==(expected is None),'centroid projection feasibility')
        if actual is None:
            continue
        need(actual['minimum_mass']==expected,'independent face enumeration gives the same exact minimum')
        need(sum(actual['point'])==0 and all(lo<=v<=hi for v,(lo,hi) in zip(actual['point'],intervals)),
            'exact continuous projection lies in the slice')
        weights=[(i+2*j+index+1)*E.S//23 for i,j in E.PAIRS]
        lines={}
        norms=[]
        for x in vertices(intervals):
            norm=sum(v*v for v in x)
            graph=sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
            lines[norm]=max(lines.get(norm,graph),graph)
            norms.append(norm)
            records+=1
            p=actual['point']
            need(sum(F(v)**2 for v in x)-actual['minimum_mass']>=sum((F(v)-u)**2 for v,u in zip(x,p)),
                'projection normal-cone inequality at independent vertices')
        generated=M.vertex_lines(intervals,weights)
        need(generated['lines']==lines and generated['feasible_vertex_records']==len(norms),
            'Gray-code graph/mass lines match direct Cartesian vertices')
        unweighted=M.vertex_lines(intervals)
        need(max(unweighted['lines'])==max(norms),'exact maximum centered mass')
    special=M.projection([(E.S,2*E.S)]+[(-E.S,0)]*7)
    need(special['minimum_mass']==F(8*E.S*E.S,7),'continuous minimum is not rounded to a grid point')
    upper_box=[((-E.S,0),(0,0)) for _ in range(8)]
    upper_box[0]=((E.S,E.S),(0,0))
    minimum,maximum=M.mass_range(upper_box)
    need(minimum==F(8*E.S*E.S,7) and maximum==2*E.S*E.S<7*E.S*E.S,
         'centered maximum mass rejects a box whose independent squared upper sum is eight')
    lower_box=[((-E.S,E.S),(0,0)) for _ in range(8)]
    lower_box[0]=((5*E.S//2,5*E.S//2),(0,0))
    minimum,maximum=M.mass_range(lower_box)
    need(minimum==F(50*E.S*E.S,7)>7*E.S*E.S,'centered minimum mass exceeds seven exactly')
    return {'independent_projection_cases':len(cases),'independent_vertex_records':records,
        'non_grid_projection_mass':'8/7','upper_mass_rejection_example':'2',
        'lower_mass_rejection_example':'50/7','floating_point_acceptance':False}


def envelope_checks():
    E.setup(40)
    families=[{0:0,2:8,4:5},{1:2,3:10,5:1},{1:1,2:2,3:3},{3:7},{4:3},
        {0:8,1:6,2:4,3:2,4:0}]
    queries=0
    for lines in families:
        hull=M.upper_hull(lines)
        points=[F(v) for v in range(-20,21)]
        for slope,intercept,start in hull[1:]:
            points += [start-F(1,1000),start,start+F(1,1000)]
        for point in points:
            need(M.hull_value(hull,point)==max(graph-point*mass for mass,graph in lines.items()),
                'upper envelope equals direct maximum at integer and rational parameters')
            queries+=1
    checked=0
    for left,right in [(families[0],families[1]),(families[3],families[4]),
                       (families[0],families[3]),(families[5],families[1])]:
        lx={mass*E.S*E.S:graph*E.S*E.S for mass,graph in left.items()}
        ly={mass*E.S*E.S:graph*E.S*E.S for mass,graph in right.items()}
        for limit in [-8,-1,0,3,10,50]:
            result=M.minimize_dual(lx,ly,limit)
            if max(left)+max(right)<7:
                need(result is None,'unbounded left dual is a normalization infeasibility case')
                continue
            def direct(t):
                return 7*t*E.S*E.S+max(g-t*m for m,g in lx.items())+max(g-t*m for m,g in ly.items())
            expected=min(direct(t) for t in range(-100,limit+1))
            need(result['value']==expected and result['multiplier']<=limit
                 and direct(result['multiplier'])==expected,'integer-grid optimum equals independent exhaustive parameter scan')
            checked+=1
    need(M.minimize_dual({E.S*E.S:0},{E.S*E.S:0},0) is None,'mass upper bound below seven is not a finite dual minimum')
    return {'independent_upper_envelope_queries':queries,'independent_grid_minimization_cases':checked,
        'parameters_checked_per_scan_at_least':93,'floating_point_acceptance':False}


def physical_check():
    from verify import C,exact_log_product
    E.setup(40)
    q=[C(v) for v in [F(6,5),F(-6,5),F(9,10),F(-9,10),F(1),F(-1),F(1,2),F(-1,2)]]
    values=[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]]
    center=[round(v*E.S) for v in values]
    radius=E.S//12
    box=[(v-radius,v+radius) for v in center]
    cp=E.positions([E.point(v) for v in center])
    qbox=[((z[0][0]-radius,z[0][0]+radius),(z[1][0]-radius,z[1][0]+radius)) for z in cp]
    qbox[0]=box[0],E.ZERO
    new=M.target_upper(box,qbox,E.model()['radius'],center)
    old=R.centered_target_upper(box,qbox,E.model()['radius'],center)
    lower,_=exact_log_product(q)
    need(new is not None and old is not None and lower<=F(new['log_Delta_over_target_upper'],E.S)
        and new['log_Delta_over_target_upper']<=old['log_Delta_over_target_upper']<0,
        'actual normalized point and whole box improve or match the previous two-parameter bound')
    point=[C(1),C(-1),C(0,1),C(0,-1),C(1,F(1,2)),C(-1,-F(1,2)),C(F(1,2)),C(-F(1,2))]
    values=[point[0].re]+[v for z in point[1:7] for v in [z.re,z.im]]
    exact_center=[int(v*E.S) for v in values]
    point_box=[E.point(v) for v in exact_center]
    tight=M.target_upper(point_box,E.positions(point_box),E.model()['radius'],exact_center)
    exact_lower,exact_upper=exact_log_product(point)
    need(tight is not None and exact_lower<=F(tight['log_Delta_over_target_upper'],E.S)<exact_upper+F(1,10**6),
         'self-reference point checks the exact normalization and rounding scale')
    return {'actual_center_and_mass_exact':True,'whole_box_log_upper':str(F(new['log_Delta_over_target_upper'],E.S)),
        'previous_two_parameter_log_upper':str(F(old['log_Delta_over_target_upper'],E.S)),
        'parameter_candidates':new['dual_candidate_count'],'self_reference_point_tight':True,
        'benchmark_is_a_remaining_stationary_branch':False}


def guards():
    base=json.loads(gzip.decompress((ROOT/'joint-closure/certificates/partial_cover.json.gz').read_bytes()))
    toy={'base_sha256':R.BASE_SHA256,'precision_bits':40,'settings':base['settings'],
        'roots':[[i,j] for j,(i,b) in enumerate(base['pending'])],
        'nodes':[{'kind':'unknown'} for _ in base['pending']],
        'pending':[[j,b] for j,(i,b) in enumerate(base['pending'])],
        'complete_cover':False,'full_eight_point_exclusion_proved':False}
    rejected=0
    with tempfile.TemporaryDirectory(prefix='mass-dual-guard-',dir=ROOT) as temp:
        path=Path(temp)/'trace.gz'
        for reason in ['centered_mass_range','mass_envelope_log']:
            data=json.loads(json.dumps(toy))
            data['pending'].pop(0)
            data['nodes'][0]={'kind':'excluded','reason':reason,'algorithm':5}
            path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
            try:
                check(path,False)
            except ArithmeticError:
                rejected+=1
            else:
                raise ArithmeticError('unjustified mass/dual exclusion accepted')
    need(rejected==2,'both false mass/dual exclusions rejected')
    return {'false_exclusions_rejected':rejected}


if __name__=='__main__':
    print(json.dumps({'passed':True,'projection_and_mass_lines':projection_and_lines(),
        'grid_envelope':envelope_checks(),'physical_box':physical_check(),'coverage_guards':guards(),
        'full_eight_point_exclusion_proved':False},indent=2))
