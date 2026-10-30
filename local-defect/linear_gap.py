"""Effective linear deficit from a maximal close-pair matching and local holes."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import sys
from packing_constants import CLOSE,NEAR,GREEDY,KAPPA,PAIR_COST,HOLE_COST

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'joint-bounds'))
from subset_moment import log_interval,outward_decimal


def require_n(n):
    if type(n)is not int or n<2:raise ValueError('integer n>=2 required')


def root_up(q,scale=10**12):
    q=F(q)
    if q<0 or type(scale)is not int or scale<1:raise ValueError('q>=0 and integer scale>=1')
    a,b=q.numerator*scale*scale,q.denominator
    k=isqrt(a//b)
    return F(k+int(k*k*b<a),scale)


def linear_gap_lower(n,scale=10**12):
    require_n(n)
    return KAPPA*max(F(0),n-11*root_up(n,scale))


def ledger_lower(n,matching_size,quiet_population):
    """Energy of a VERIFIED matching/quiet plan, not arbitrary inferred counts."""
    require_n(n)
    if type(matching_size)is not int or not 0<=2*matching_size<=n or type(quiet_population)is not int or not 0<=quiet_population<=n-2*matching_size:
        raise ValueError('valid integer matching size and quiet population required')
    holes=(quiet_population+16)//17
    return {'matching_size':matching_size,'quiet_population':quiet_population,
        'disjoint_holes_guaranteed':holes,
        'energy_lower':PAIR_COST*matching_size+HOLE_COST*holes,
        'scope':'Requires a maximal close-pair matching, quiet interior centers and the stated separated/radius-capped geometry. Counts are not sufficient without those hypotheses.'}


def matching_plan(points,max_points=300,imaginary_unit_squared=1):
    """Exact conditional geometric plan for centered/renormalized rational input.

    Global optimizers satisfy the verified geometric hypotheses, but passing
    them here does not certify optimality or stationarity of this input.
    Hole centers are furnished by the analytic local alternative, not sampled.
    """
    pts=[tuple(map(F,p)) for p in points];n=len(pts)
    D=F(imaginary_unit_squared)
    if D<=0:raise ValueError('positive quadratic-field metric required')
    require_n(n)
    if n<4 or n>max_points or any(len(p)!=2 for p in pts) or len(set(pts))!=n:
        raise ValueError('4..budget distinct rational point pairs required')
    mean=tuple(sum((p[j] for p in pts),F(0))/n for j in range(2))
    v=[(p[0]-mean[0],p[1]-mean[1]) for p in pts]
    Q=sum((x*x+D*y*y for x,y in v),F(0))
    f=F(n*(n-1),2)/Q
    radii=[f*(x*x+D*y*y) for x,y in v]
    distance={(i,j):f*((v[i][0]-v[j][0])**2+D*(v[i][1]-v[j][1])**2)
              for i in range(n) for j in range(i+1,n)}
    if min(distance.values())<2 or max(radii)>n:
        raise ValueError('input fails separated radius-capped geometry')
    def d(i,j):return F(0) if i==j else distance[min(i,j),max(i,j)]
    endpoints=set();matching=[]
    for (i,j),q in distance.items():
        if q<CLOSE*CLOSE and i not in endpoints and j not in endpoints:
            matching.append((i,j));endpoints.update((i,j))
    blocked={i for i in range(n) if any(d(i,j)<=NEAR*NEAR for j in endpoints)}
    if len(blocked)>32*len(matching):raise ArithmeticError('endpoint packing cap failed')
    boundary=set()
    a=GREEDY/2
    for i,q in enumerate(radii):
        margin=n+a*a-q
        if margin<0 or margin*margin<4*a*a*n:boundary.add(i)
    quiet=[i for i in range(n) if i not in blocked and i not in boundary]
    selected=[];removed=set()
    for i in quiet:
        if i not in removed:
            selected.append(i);removed.update(j for j in quiet if d(i,j)<=GREEDY*GREEDY)
    if 17*len(selected)<len(quiet):raise ArithmeticError('greedy packing cap failed')
    if 32*len(matching)+len(quiet)+len(boundary)<n:raise ArithmeticError('population accounting failed')
    remaining=[i for i in range(n) if i not in endpoints]
    if any(d(i,j)<CLOSE*CLOSE for i in remaining for j in remaining if i<j):
        raise ArithmeticError('matching not maximal')
    return {'n':n,'force_scale_squared':f,'imaginary_unit_squared':D,'minimum_squared_distance':min(distance.values()),
        'maximum_squared_radius':max(radii),'matching':matching,'boundary_population':len(boundary),
        'blocked_population':len(blocked),'quiet_population':len(quiet),
        'selected_quiet_centers':selected,
        'conditional_energy_lower':PAIR_COST*len(matching)+HOLE_COST*len(selected),
        'universal_linear_lower':linear_gap_lower(n),
        'scope':'Geometric necessary-input plan only. Empty-ball locations are existential from the proved local alternative; no global optimizer is inferred.'}


def gamma_upper_interval(places=30):
    lo,hi=log_interval(2,80)
    return outward_decimal((F(3,4)-hi/2-KAPPA,F(3,4)-lo/2-KAPPA),places)
