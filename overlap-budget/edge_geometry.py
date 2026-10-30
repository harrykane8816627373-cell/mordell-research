"""Exact vertex-cover plan using all short edges and a c-separated remainder."""
from fractions import Fraction as F
from parameter_certificate import family_certificate,gap_lower


def vertex_plan(points,name='strong37',imaginary_unit_squared=1,max_points=300):
    pts=[tuple(map(F,p)) for p in points];n=len(pts);D=F(imaginary_unit_squared)
    if type(max_points)is not int or not 4<=n<=max_points or D<=0 or any(len(p)!=2 for p in pts) or len(set(pts))!=n:
        raise ValueError('distinct rational quadratic-field point pairs, n>=4 within budget required')
    cert=family_certificate(name)
    mean=tuple(sum((p[j] for p in pts),F(0))/n for j in range(2))
    v=[(p[0]-mean[0],p[1]-mean[1]) for p in pts]
    Q=sum((x*x+D*y*y for x,y in v),F(0));f=F(n*(n-1),2)/Q
    radii=[f*(x*x+D*y*y) for x,y in v]
    pairs={(i,j):f*((v[i][0]-v[j][0])**2+D*(v[i][1]-v[j][1])**2)
           for i in range(n) for j in range(i+1,n)}
    if min(pairs.values())<2 or max(radii)>n:raise ValueError('input fails global separated/radius-capped geometry')
    def distance(i,j):return F(0) if i==j else pairs[min(i,j),max(i,j)]
    c,U,L=[cert[k] for k in ('close_cutoff','near_radius','greedy_radius')]
    edges=[ij for ij,d in pairs.items() if d<=c*c]
    # One endpoint for each edge, with union size <= number of edges.
    cover=sorted({i for i,j in edges})
    if len(cover)>len(edges) or any(i not in cover and j not in cover for i,j in edges):
        raise ArithmeticError('constructed vertex cover failed')
    blocked={i for i in range(n) if any(distance(i,j)<=U*U for j in cover)}
    if len(blocked)>cert['endpoint_cap']*len(cover):raise ArithmeticError('endpoint cap failed')
    a=cert['query_radius']+cert['hole_radius'];boundary=set()
    for i,r2 in enumerate(radii):
        margin=n+a*a-r2
        if margin<0 or margin*margin<4*a*a*n:boundary.add(i)
    quiet=[i for i in range(n) if i not in blocked and i not in boundary]
    if any(distance(i,j)<=c*c for i in quiet for j in quiet if i<j):raise ArithmeticError('quiet set not c-separated')
    selected=[];removed=set()
    for i in quiet:
        if i not in removed:
            selected.append(i);removed.update(j for j in quiet if distance(i,j)<=L*L)
    if cert['quiet_cap']*len(selected)<len(quiet):raise ArithmeticError('quiet greedy cap failed')
    if cert['endpoint_cap']*len(edges)+len(quiet)+len(boundary)<n:raise ArithmeticError('vertex-cover accounting failed')
    energy=cert['pair_cost_used']*len(edges)+cert['hole_cost']*len(selected)
    if energy<gap_lower(n,name):raise ArithmeticError('actual ledger below universal bound')
    return {'n':n,'certificate':name,'imaginary_unit_squared':D,'force_scale_squared':f,
        'short_edge_count':len(edges),'cover_size':len(cover),'blocked_population':len(blocked),
        'boundary_population':len(boundary),'quiet_population':len(quiet),
        'selected_quiet_centers':selected,'conditional_energy_lower':energy,
        'universal_gap_lower':gap_lower(n,name),
        'scope':'Input geometric checks and conditional energy plan only. No optimizer or stationary status is inferred; actual hole locations follow from the local geometric lemma.'}
