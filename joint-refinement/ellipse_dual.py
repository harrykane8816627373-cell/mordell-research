"""Exact joint use of the old moment ellipse and centered graph dual bounds."""
from fractions import Fraction as F
import refinement as R
import mass_dual as D

E=R.E


def mass_split(qbox):
    ranges=[]
    for component in [0,1]:
        intervals=[q[component] for q in qbox]
        minimum=D.projection(intervals)
        maximum=D.vertex_lines(intervals)
        if minimum is None or maximum is None:
            return None
        ranges.append((minimum['minimum_mass'],max(maximum['lines'])))
    dot=E.ZERO
    s3=E.ci(0)
    for q in qbox:
        dot=E.add(dot,E.mul(q[0],q[1]))
        s3=E.ca(s3,E.cm(E.cm(q,q),q))
    residual_upper=E.S-54*E.cnorm(s3)[0]//175-4*E.square(dot)[0]
    if residual_upper<0:
        return {'empty':True,'reason':'negative_ellipse_budget'}
    root=E.sqrt_hi(residual_upper)
    lo=max(ranges[0][0],7*E.S**2-ranges[1][1],F((7*E.S-root)*E.S,2),3*E.S**2)
    hi=min(ranges[0][1],7*E.S**2-ranges[1][0],F((7*E.S+root)*E.S,2),4*E.S**2)
    return {'empty':lo>hi,'reason':'disjoint_component_mass_ranges' if lo>hi else None,
        'lower':lo,'upper':hi,'component_ranges':ranges,'ellipse_residual_upper':residual_upper}


def break_candidates(hulls,limit):
    out={limit}
    for hull in hulls:
        for slope,intercept,start in hull[1:]:
            a=start.numerator//start.denominator
            b=-(-start.numerator//start.denominator)
            out.update(t for t in [a,b] if t<=limit)
    return out


def one_minimum(hull,mass,limit):
    candidates=break_candidates([hull],limit)
    if mass>-hull[0][0]:
        raise ValueError('mass exceeds maximum recorded vertex mass')
    return min(candidates,key=lambda t:(D.hull_value(hull,t)+mass*t,t))


def ordered_minimum(hx,hy,mx,my,tx,ty):
    """Minimize Fx(a)+Fy(b), a<=tx, b<=ty, b<=a, on the integer grid."""
    anchor=one_minimum(hy,my,ty)
    candidates=break_candidates([hx,hy],tx)
    candidates.add(min(tx,anchor))
    def value(a):
        b=min(a,anchor)
        return D.hull_value(hx,a)+mx*a+D.hull_value(hy,b)+my*b
    a=min(candidates,key=lambda t:(value(t),t))
    b=min(a,anchor)
    return {'value':value(a),'parameters':[a,b],'candidates':len(candidates)}


def minimize_split(lines_x,lines_y,limit_x,limit_y,lower,upper):
    if any(type(t) is not int for t in [limit_x,limit_y]):
        raise ValueError('exact integer parameter limits')
    if lower>upper or lower<0 or upper>7*E.S**2:
        raise ValueError('ordered physical mass-split interval')
    if upper>max(lines_x) or 7*E.S**2-lower>max(lines_y):
        raise ValueError('component mass bounds must respect maximum vertex masses')
    hx,hy=D.upper_hull(lines_x),D.upper_hull(lines_y)
    # a>=b uses the upper mass endpoint; a<=b uses the lower one.
    first=ordered_minimum(hx,hy,upper,7*E.S**2-upper,limit_x,limit_y)
    swapped=ordered_minimum(hy,hx,7*E.S**2-lower,lower,limit_y,limit_x)
    second=dict(swapped,parameters=list(reversed(swapped['parameters'])))
    chosen=min([first,second],key=lambda result:(result['value'],result['parameters']))
    return dict(chosen,total_candidates=first['candidates']+second['candidates'],
        hull_segments=[len(hx),len(hy)])


def target_upper(box,qbox,radius,midpoint=None,split=None):
    if split is None:
        split=mass_split(qbox)
    if split is None or split['empty']:
        return None
    if midpoint is None:
        midpoint,_=E.slab_center(box,radius,qbox[7])
    if len(midpoint)!=13 or any(type(v) is not int for v in midpoint):
        raise ValueError('exact thirteen-coordinate reference')
    center=E.positions([E.point(v) for v in midpoint])
    product,weights=1,[]
    for i,j in E.PAIRS:
        dx=center[i][0][0]-center[j][0][0]
        dy=center[i][1][0]-center[j][1][0]
        distance=dx*dx+dy*dy
        if distance==0:
            return None
        product*=distance
        weights.append(E.rat(E.S**2,distance)[1])
    rows=[0]*8
    graph=[[0]*8 for _ in range(8)]
    for (i,j),w in zip(E.PAIRS,weights):
        rows[i]+=w
        rows[j]+=w
        graph[i][j]=graph[j][i]=w
    pairs=[]
    limits=[]
    records=[]
    for component in [0,1]:
        active=[i for i in range(8) if qbox[i][component][0]!=qbox[i][component][1]]
        cp=[(i,j) for i in active for j in active if i<j]
        pairs.append(cp)
        limits.append(min((rows[i]+rows[j]+2*graph[i][j])//2 for i,j in cp) if cp else 4*E.S)
        records.append(D.vertex_lines([q[component] for q in qbox],weights))
    result=minimize_split(records[0]['lines'],records[1]['lines'],*limits,split['lower'],split['upper'])
    for component in [0,1]:
        if any(rows[i]+rows[j]+2*graph[i][j]-2*result['parameters'][component]<0 for i,j in pairs[component]):
            raise ArithmeticError('each component must retain exchange convexity')
    raw=F(result['value'],E.S**2)
    upper=R.log_ratio(product,7**7*E.S**56)[1]-28*E.S-(-raw.numerator//raw.denominator)
    return {'log_Delta_over_target_upper':upper,'mass_parameters':result['parameters'],
        'component_exchange_limits':limits,'dual_candidate_count':result['total_candidates'],
        'mass_split_lower_raw':str(split['lower']),'mass_split_upper_raw':str(split['upper']),
        'dual_hull_segments':result['hull_segments']}


def contract(box,settings,passes=5):
    box,reason,aux=R.contract(box,settings,passes,algorithm=5)
    if box is None:
        return box,reason,aux
    from radial_stability import exclusion
    reason,evidence=exclusion(aux[0])
    if reason is not None:
        return None,reason,evidence
    split=mass_split(aux[0])
    if split is not None and split['empty']:
        return None,'ellipse_mass_split',{'split_reason':split['reason']}
    candidates=[('midpoint',None)]
    midpoint=[(a+b)//2 for a,b in box]
    nearest=sorted(enumerate(R.REFERENCES['centers']),
        key=lambda item:sum((a-b)**2 for a,b in zip(midpoint,item[1])))[:2]
    candidates += [('reference_'+str(i),c) for i,c in nearest]
    for source,center in candidates:
        upper=target_upper(box,aux[0],settings['radius'],center,split)
        if upper is not None and upper['log_Delta_over_target_upper']<0:
            upper['center_source']=source
            return None,'ellipse_split_log',upper
    return box,None,aux
