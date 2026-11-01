"""Exact normalization feasibility and grid-dual envelopes on centered boxes."""
from fractions import Fraction as F
import refinement as R

E=R.E


def projection(intervals):
    if len(intervals)!=8 or any(type(v) is not int for p in intervals for v in p):
        raise ValueError('eight exact endpoint pairs')
    if any(lo>hi for lo,hi in intervals):
        raise ValueError('ordered intervals')
    if sum(lo for lo,hi in intervals)>0 or sum(hi for lo,hi in intervals)<0:
        return None
    endpoints=sorted({v for p in intervals for v in p})
    values=[sum(min(hi,max(lo,t)) for lo,hi in intervals) for t in endpoints]
    tau=None
    for t,value in zip(endpoints,values):
        if value==0:
            tau=F(t)
            break
    if tau is None:
        for a,b,fa,fb in zip(endpoints,endpoints[1:],values,values[1:]):
            if fa<0<fb:
                middle=a+b
                free=[i for i,(lo,hi) in enumerate(intervals) if 2*lo<middle<2*hi]
                fixed=sum(lo if middle<=2*lo else hi for i,(lo,hi) in enumerate(intervals) if i not in free)
                tau=F(-fixed,len(free))
                break
    if tau is None:
        raise ArithmeticError('a feasible centroid slice has a projection threshold')
    point=[min(F(hi),max(F(lo),tau)) for lo,hi in intervals]
    if sum(point)!=0:
        raise ArithmeticError('exact projection centroid')
    return {'point':point,'tau':tau,'minimum_mass':sum(v*v for v in point)}


def vertex_lines(intervals,weights=None):
    if len(intervals)!=8 or any(type(v) is not int for p in intervals for v in p) or any(lo>hi for lo,hi in intervals):
        raise ValueError('eight ordered intervals')
    if sum(lo for lo,hi in intervals)>0 or sum(hi for lo,hi in intervals)<0:
        return None
    if weights is not None and (len(weights)!=28 or any(type(v) is not int or v<0 for v in weights)):
        raise ValueError('28 nonnegative exact weights')
    matrix=[[0]*8 for _ in range(8)]
    if weights is not None:
        for (i,j),w in zip(E.PAIRS,weights):
            matrix[i][i]+=w
            matrix[j][j]+=w
            matrix[i][j]=matrix[j][i]=-w
    active=[i for i,(lo,hi) in enumerate(intervals) if lo!=hi]
    lines={}
    tested,feasible=0,0

    def record(x,value,mass):
        nonlocal tested,feasible
        tested+=1
        if any(not lo<=v<=hi for v,(lo,hi) in zip(x,intervals)):
            return
        feasible+=1
        lines[mass]=max(lines.get(mass,value),value)

    if not active:
        x=[lo for lo,hi in intervals]
        value=0 if weights is None else sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
        record(x,value,sum(v*v for v in x))
    for free in active:
        other=[i for i in active if i!=free]
        x=[lo for lo,hi in intervals]
        x[free]=-sum(x[i] for i in range(8) if i!=free)
        state=[0]*8
        row=[sum(matrix[i][j]*x[j] for j in range(8)) for i in range(8)]
        value=0 if weights is None else sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
        mass=sum(v*v for v in x)
        record(x,value,mass)
        for step in range(1,1 << len(other)):
            k=other[(step & -step).bit_length()-1]
            nxt=1-state[k]
            d=intervals[k][nxt]-x[k]
            if weights is not None:
                value+=2*d*(row[k]-row[free])+d*d*(matrix[k][k]+matrix[free][free]-2*matrix[k][free])
            mass+=2*d*(x[k]-x[free])+2*d*d
            x[k]+=d
            x[free]-=d
            state[k]=nxt
            if weights is not None:
                for j in range(8):
                    row[j]+=d*(matrix[j][k]-matrix[j][free])
            record(x,value,mass)
    if not feasible:
        return None
    return {'lines':lines,'vertex_records_checked':tested,'feasible_vertex_records':feasible}


def mass_range(qbox):
    minimum=F(0)
    maximum=0
    for component in [0,1]:
        intervals=[z[component] for z in qbox]
        p=projection(intervals)
        v=vertex_lines(intervals)
        if p is None or v is None:
            return None
        minimum+=p['minimum_mass']
        maximum+=max(v['lines'])
    return minimum,maximum


def upper_hull(lines):
    # Each record denotes graph - t*mass. Slopes increase along the envelope.
    hull=[]
    for slope,intercept in sorted((-mass,graph) for mass,graph in lines.items()):
        start=None
        while hull:
            previous=hull[-1]
            start=F(previous[1]-intercept,slope-previous[0])
            if len(hull)>1 and start<=previous[2]:
                hull.pop()
            else:
                break
        hull.append((slope,intercept,start))
    return hull


def hull_value(hull,t):
    lo,hi=0,len(hull)
    while lo+1<hi:
        mid=(lo+hi)//2
        if hull[mid][2]<=t:
            lo=mid
        else:
            hi=mid
    slope,intercept,_=hull[lo]
    return slope*t+intercept


def minimize_dual(lines_x,lines_y,limit):
    if type(limit) is not int:
        raise ValueError('integer parameter limit')
    if max(lines_x)+max(lines_y)<7*E.S*E.S:
        return None
    hx,hy=upper_hull(lines_x),upper_hull(lines_y)
    candidates={limit,min(4*E.S,limit)}
    for h in [hx,hy]:
        for slope,intercept,start in h[1:]:
            floor=start.numerator//start.denominator
            ceiling=-(-start.numerator//start.denominator)
            candidates.update(t for t in [floor,ceiling] if t<=limit)
    chosen=min(candidates,key=lambda t:(7*t*E.S*E.S+hull_value(hx,t)+hull_value(hy,t),t))
    value=7*chosen*E.S*E.S+hull_value(hx,chosen)+hull_value(hy,chosen)
    return {'value':value,'multiplier':chosen,'candidate_count':len(candidates),
        'hull_segments':[len(hx),len(hy)]}


def target_upper(box,qbox,radius,midpoint=None):
    if midpoint is None:
        midpoint,_=E.slab_center(box,radius,qbox[7])
    if len(midpoint)!=13 or not all(type(v) is int for v in midpoint):
        raise ValueError('exact rational reference center')
    center=E.positions([E.point(v) for v in midpoint])
    product,weights=1,[]
    for i,j in E.PAIRS:
        dx=center[i][0][0]-center[j][0][0]
        dy=center[i][1][0]-center[j][1][0]
        d=dx*dx+dy*dy
        if d==0:
            return None
        product*=d
        weights.append(E.rat(E.S*E.S,d)[1])
    rows=[0]*8
    graph=[[0]*8 for _ in range(8)]
    for (i,j),w in zip(E.PAIRS,weights):
        rows[i]+=w
        rows[j]+=w
        graph[i][j]=graph[j][i]=w
    pairs=set()
    for component in [0,1]:
        active=[i for i in range(8) if qbox[i][component][0]!=qbox[i][component][1]]
        pairs.update((i,j) for i in active for j in active if i<j)
    limit=min((rows[i]+rows[j]+2*graph[i][j])//2 for i,j in pairs) if pairs else 4*E.S
    records=[vertex_lines([q[component] for q in qbox],weights) for component in [0,1]]
    if any(v is None for v in records):
        return None
    result=minimize_dual(records[0]['lines'],records[1]['lines'],limit)
    if result is None:
        return None
    if any(rows[i]+rows[j]+2*graph[i][j]-2*result['multiplier']<0 for i,j in pairs):
        raise ArithmeticError('chosen grid-dual multiplier must preserve exchange convexity')
    upper=R.log_ratio(product,7**7*E.S**56)[1]-28*E.S-(-result['value']//(E.S*E.S))
    return {'log_Delta_over_target_upper':upper,'mass_multiplier':result['multiplier'],
        'exchange_convexity_limit':limit,'dual_candidate_count':result['candidate_count'],
        'dual_hull_segments':result['hull_segments'],
        'vertex_records_checked':sum(v['vertex_records_checked'] for v in records),
        'feasible_vertex_records':sum(v['feasible_vertex_records'] for v in records)}
