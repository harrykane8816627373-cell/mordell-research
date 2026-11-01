"""Integer acceptance of mass-adjusted logarithmic Taylor exclusion."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'joint-closure'))
import interval_engine as E

BASE_SHA256 = '9b40704b8be9b15a219c89b64989bd71caddf74564386570900164494f857b58'
REFERENCES = json.loads((Path(__file__).resolve().parent/'reference_centers.json').read_text())


def _log_reduced(a, b):
    """Enclose log(a/b) for 1 <= a/b <= 2 by a positive atanh series."""
    if not (0 < b <= a <= 2*b):
        raise ValueError('reduced positive ratio')
    if a == b:
        return E.ZERO
    x = E.rat(a-b, a+b)
    x2, power, result = E.mul(x, x), x, E.ZERO
    terms = 40
    for j in range(terms):
        result = E.add(result, E.divide(E.scale(power, 2), 2*j+1))
        power = E.mul(power, x2)
    # x <= 1/3 exactly. Positive tail <= 9/[4(2m+1)3^(2m+1)].
    tail = E.rat(9, 4*(2*terms+1)*3**(2*terms+1))[1]
    return result[0], result[1]+tail


def log_ratio(a, b):
    """Integer fixed-point interval for log(a/b), for arbitrary positive integers."""
    if a <= 0 or b <= 0:
        raise ValueError('positive logarithm arguments')
    exponent = a.bit_length()-b.bit_length()
    num, den = (a, b << exponent) if exponent >= 0 else (a << -exponent, b)
    if num < den:
        exponent -= 1
        num *= 2
    if not den <= num < 2*den:
        raise ArithmeticError('binary logarithm reduction')
    return E.add(_log_reduced(num, den), E.scale(_log_reduced(2, 1), exponent))


def target_log_upper(box, qbox, radius, midpoint=None):
    """Upper log(Delta/7^7) for every mass-seven actual point in these slabs.

    None means the Taylor disks were too wide; it never means exclusion.
    """
    if midpoint is None:
        midpoint, _ = E.slab_center(box, radius, qbox[7])
    if len(midpoint) != 13 or not all(type(v) is int for v in midpoint):
        raise ValueError('exact rational center in 13 centered-coordinate variables')
    center = E.positions([E.point(x) for x in midpoint])
    raw_norms = [z[0][0]**2+z[1][0]**2 for z in center]
    product = 1
    remainder_upper = 0
    for i, j in E.PAIRS:
        dx = center[i][0][0]-center[j][0][0]
        dy = center[i][1][0]-center[j][1][0]
        raw_distance = dx*dx+dy*dy
        if not raw_distance:
            return None
        product *= raw_distance
        distance_lower = E.sqrt_lo(raw_distance//E.S)
        if distance_lower <= 0:
            return None
        displacement = E.cs(E.cs(qbox[i], qbox[j]), E.cs(center[i], center[j]))
        displacement_upper = E.sqrt_hi(E.cnorm(displacement)[1])
        rho = -(-displacement_upper*E.S//distance_lower)
        if rho >= E.S:
            return None
        remainder_upper += E.mul(E.square((0, rho)),
            E.reciprocal((E.S-rho, E.S-rho)))[1]
    inverse, _ = E.invert_edges(center)
    gradient, _ = E.force_gradient(center, inverse)
    gradient = E.reduced_gradient(gradient)
    linear_upper = 0
    for g, interval, c in zip(gradient, box, midpoint):
        linear_upper += E.scale(E.mul(E.neg(g), E.sub(interval, E.point(c))), 2)[1]
    mass_correction = E.scale(E.sub(E.rat(7), E.rat(sum(raw_norms), E.S**2)), 4)
    center_log = log_ratio(product, 7**7*E.S**56)
    upper = center_log[1]+mass_correction[1]+linear_upper+remainder_upper
    return {'log_Delta_over_target_upper': upper,
        'center_log_upper': center_log[1], 'mass_correction_upper': mass_correction[1],
        'linear_residual_upper': linear_upper, 'Taylor_remainder_upper': remainder_upper}


def _quadratic_max(d, b, interval):
    """Upper of d*x*x+b*x on an integer-grid real interval; all b in a box."""
    lo, hi = interval
    best = None
    for slope in b:
        for x in [lo, hi]:
            value = -(-(d*x*x+slope*x*E.S)//(E.S*E.S))
            best = value if best is None else max(best, value)
        if d < 0:
            numerator, denominator = slope*E.S, -2*d
            if lo*denominator <= numerator <= hi*denominator:
                value = -(-(slope*slope)//(-4*d))
                best = max(best, value)
    return best


def quadratic_target_upper(box, qbox, radius, midpoint=None):
    """Global mass-preserving log majorant; no small-displacement condition."""
    if midpoint is None:
        midpoint, _ = E.slab_center(box, radius, qbox[7])
    if len(midpoint) != 13 or not all(type(v) is int for v in midpoint):
        raise ValueError('exact rational center')
    center = E.positions([E.point(x) for x in midpoint])
    delta = [E.cs(z, c) for z, c in zip(qbox, center)]
    product, weights = 1, []
    for i, j in E.PAIRS:
        dx = center[i][0][0]-center[j][0][0]
        dy = center[i][1][0]-center[j][1][0]
        squared = dx*dx+dy*dy
        if not squared:
            return None
        product *= squared
        weights.append(E.rat(E.S**2, squared))
    try:
        inverse, _ = E.invert_edges(center)
    except ZeroDivisionError:
        return None
    gradient, _ = E.force_gradient(center, inverse)
    reduced = E.reduced_gradient(gradient)
    linear = sum(E.scale(E.mul(E.neg(g), E.sub(a, E.point(c))), 2)[1]
        for g, a, c in zip(reduced, box, midpoint))
    mass = sum(z[0][0]**2+z[1][0]**2 for z in center)
    constant = log_ratio(product, 7**7*E.S**56)[1]
    constant += E.scale(E.sub(E.rat(7), E.rat(mass, E.S**2)), 4)[1]
    adaptive = [max(1, E.sqrt_hi(E.cnorm(z)[1])) for z in delta]
    best = None
    for name, a in [('equal', [1]*8), ('radius', adaptive)]:
        diagonal = [-4*E.S]*8
        for (i, j), w in zip(E.PAIRS, weights):
            diagonal[i] += E.mul(w, E.rat(a[i]+a[j], a[i]))[1]
            diagonal[j] += E.mul(w, E.rat(a[i]+a[j], a[j]))[1]
        norm_bound = 0
        separable = 0
        for i in range(8):
            norm_bound += E.mul(E.point(diagonal[i]), E.cnorm(delta[i]))[1]
            for component in [0, 1]:
                separable += _quadratic_max(diagonal[i], E.scale(E.neg(gradient[i][component]), 2),
                    delta[i][component])
        variants = [('centered_linear', linear+norm_bound), ('separable', separable)]
        for mode, bound in variants:
            result = {'log_Delta_over_target_upper': constant+bound,
                'weight_scheme': name, 'quadratic_box_bound': mode,
                'negative_diagonal_coefficients': sum(d < 0 for d in diagonal)}
            if best is None or result['log_Delta_over_target_upper'] < best['log_Delta_over_target_upper']:
                best = result
    return best


def _contract_v1(box, settings, passes=5):
    box, reason, aux = E.contract(box, settings, passes=passes)
    if box is None:
        return box, reason, aux
    upper = target_log_upper(box, aux[0], settings['radius'])
    if upper is not None and upper['log_Delta_over_target_upper'] < 0:
        upper['center_source'] = 'midpoint'
        return None, 'mass_log_Taylor', upper
    if REFERENCES['precision_bits'] != E.BITS:
        raise ArithmeticError('reference grid matches interval grid')
    midpoint = [(a+b)//2 for a, b in box]
    nearest = sorted(enumerate(REFERENCES['centers']),
        key=lambda item: sum((a-b)**2 for a, b in zip(midpoint, item[1])))[:2]
    for index, center in nearest:
        upper = target_log_upper(box, aux[0], settings['radius'], center)
        if upper is not None and upper['log_Delta_over_target_upper'] < 0:
            upper['center_source'] = 'reference_'+str(index)
            return None, 'mass_log_Taylor', upper
    return box, None, aux


def vertex_quadratic_upper(intervals, slopes, weights):
    """Enumerate a coordinatewise convex majorant in at most 256 Gray-code vertices."""
    diagonal = [-4*E.S]*8
    interaction = [[0]*8 for _ in range(8)]
    for (i,j),w in zip(E.PAIRS,weights):
        if type(w) is not int or w < 0:
            raise ValueError('nonnegative exact upper edge weights')
        diagonal[i] += w
        diagonal[j] += w
        interaction[i][j] = interaction[j][i] = -2*w
    constant, clipped = 0, 0
    for i,d in enumerate(diagonal):
        if d < 0:
            lo,hi = intervals[i]
            minimum = 0 if lo <= 0 <= hi else min(lo*lo,hi*hi)
            constant += d*minimum
            diagonal[i] = 0
            clipped += 1
    values=[]
    for d,b,interval in zip(diagonal,slopes,intervals):
        values.append([d*x*x+E.S*(b[1] if x >= 0 else b[0])*x for x in interval])
    x=[v[0] for v in intervals]
    state=[0]*8
    rows=[sum(interaction[i][j]*x[j] for j in range(8)) for i in range(8)]
    value=constant+sum(v[0] for v in values)+sum(interaction[i][j]*x[i]*x[j] for i,j in E.PAIRS)
    best=value
    active=[i for i,v in enumerate(intervals) if v[0] != v[1]]
    for step in range(1,1 << len(active)):
        k=active[(step & -step).bit_length()-1]
        next_state=1-state[k]
        change=intervals[k][next_state]-x[k]
        value += values[k][next_state]-values[k][state[k]]+change*rows[k]
        x[k] += change
        state[k]=next_state
        for j in range(8):
            if j != k:
                rows[j] += interaction[j][k]*change
        best=max(best,value)
    return -(-best//(E.S*E.S)), 1 << len(active), clipped


def vertex_target_upper(box,qbox,radius,midpoint=None):
    if midpoint is None:
        midpoint,_=E.slab_center(box,radius,qbox[7])
    if len(midpoint) != 13 or not all(type(v) is int for v in midpoint):
        raise ValueError('exact rational reference center')
    center=E.positions([E.point(v) for v in midpoint])
    delta=[E.cs(z,c) for z,c in zip(qbox,center)]
    product,weights=1,[]
    for i,j in E.PAIRS:
        dx=center[i][0][0]-center[j][0][0]
        dy=center[i][1][0]-center[j][1][0]
        squared=dx*dx+dy*dy
        if not squared:
            return None
        product *= squared
        weights.append(E.rat(E.S**2,squared)[1])
    try:
        inverse,_=E.invert_edges(center)
    except ZeroDivisionError:
        return None
    gradient,_=E.force_gradient(center,inverse)
    mass=sum(z[0][0]**2+z[1][0]**2 for z in center)
    upper=log_ratio(product,7**7*E.S**56)[1]
    upper+=E.scale(E.sub(E.rat(7),E.rat(mass,E.S**2)),4)[1]
    vertices,clipped=0,0
    for component in [0,1]:
        bound,count,clip=vertex_quadratic_upper([v[component] for v in delta],
            [E.scale(E.neg(v[component]),2) for v in gradient],weights)
        upper+=bound
        vertices+=count
        clipped=max(clipped,clip)
    return {'log_Delta_over_target_upper':upper,'vertices_enumerated':vertices,
        'negative_diagonal_clips':clipped}


def centered_slice_max(intervals,weights,multipliers):
    """Maximize exchange-convex quadratic majorants on sum(x)=0 sliced boxes."""
    if len(intervals)!=8 or len(weights)!=28 or not multipliers:
        raise ValueError('eight intervals, 28 weights and at least one multiplier')
    if any(type(v) is not int or v<0 for v in weights):
        raise ValueError('nonnegative integer-grid edge weights')
    if any(type(v) is not int for v in multipliers) or any(lo>hi for lo,hi in intervals):
        raise ValueError('exact multipliers and ordered intervals')
    if sum(v[0] for v in intervals)>0 or sum(v[1] for v in intervals)<0:
        return None
    matrix=[[0]*8 for _ in range(8)]
    for (i,j),w in zip(E.PAIRS,weights):
        matrix[i][i]+=w
        matrix[j][j]+=w
        matrix[i][j]=matrix[j][i]=-w
    active=[i for i,(lo,hi) in enumerate(intervals) if lo!=hi]
    for eta in multipliers:
        if any(matrix[i][i]+matrix[j][j]-2*matrix[i][j]-2*eta<0
            for i in active for j in active if i<j):
            raise ValueError('exchange convexity is required for vertex acceptance')
    best=[None]*len(multipliers)
    tested,feasible=0,0

    def record(x,value,mass):
        nonlocal tested,feasible
        tested+=1
        if any(not lo<=v<=hi for v,(lo,hi) in zip(x,intervals)):
            return
        feasible+=1
        for k,eta in enumerate(multipliers):
            score=value-eta*mass
            best[k]=score if best[k] is None else max(best[k],score)

    if not active:
        x=[v[0] for v in intervals]
        value=sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
        record(x,value,sum(v*v for v in x))
    for free in active:
        other=[i for i in active if i!=free]
        x=[v[0] for v in intervals]
        x[free]=-sum(x[i] for i in range(8) if i!=free)
        state=[0]*8
        row=[sum(matrix[i][j]*x[j] for j in range(8)) for i in range(8)]
        value=sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
        mass=sum(v*v for v in x)
        record(x,value,mass)
        for step in range(1,1 << len(other)):
            k=other[(step & -step).bit_length()-1]
            next_state=1-state[k]
            change=intervals[k][next_state]-x[k]
            curvature=matrix[k][k]+matrix[free][free]-2*matrix[k][free]
            value+=2*change*(row[k]-row[free])+change*change*curvature
            mass+=2*change*(x[k]-x[free])+2*change*change
            x[k]+=change
            x[free]-=change
            state[k]=next_state
            for j in range(8):
                row[j]+=change*(matrix[j][k]-matrix[j][free])
            record(x,value,mass)
    if feasible==0:
        return None
    return {'upper':[-(-v//(E.S*E.S)) for v in best],
        'vertex_records_checked':tested,'feasible_vertex_records':feasible}


def centered_target_upper(box,qbox,radius,midpoint=None):
    """Global log bound retaining exact centroid slices and a mass multiplier."""
    if midpoint is None:
        midpoint,_=E.slab_center(box,radius,qbox[7])
    if len(midpoint)!=13 or not all(type(v) is int for v in midpoint):
        raise ValueError('exact rational reference center')
    center=E.positions([E.point(v) for v in midpoint])
    product,weights=1,[]
    for i,j in E.PAIRS:
        dx=center[i][0][0]-center[j][0][0]
        dy=center[i][1][0]-center[j][1][0]
        squared=dx*dx+dy*dy
        if not squared:
            return None
        product*=squared
        weights.append(E.rat(E.S**2,squared)[1])
    rows=[0]*8
    matrix=[[0]*8 for _ in range(8)]
    for (i,j),w in zip(E.PAIRS,weights):
        rows[i]+=w
        rows[j]+=w
        matrix[i][j]=matrix[j][i]=w
    pairs=set()
    for component in [0,1]:
        active=[i for i in range(8) if qbox[i][component][0]!=qbox[i][component][1]]
        pairs.update((i,j) for i in active for j in active if i<j)
    limit=min((rows[i]+rows[j]+2*matrix[i][j])//2 for i,j in pairs) if pairs else 4*E.S
    multipliers=sorted(set([min(4*E.S,limit),limit]))
    results=[]
    for component in [0,1]:
        result=centered_slice_max([v[component] for v in qbox],weights,multipliers)
        if result is None:
            return None
        results.append(result)
    constant=log_ratio(product,7**7*E.S**56)[1]-28*E.S
    candidates=[constant+7*eta+results[0]['upper'][k]+results[1]['upper'][k]
        for k,eta in enumerate(multipliers)]
    chosen=min(range(len(candidates)),key=lambda k:candidates[k])
    return {'log_Delta_over_target_upper':candidates[chosen],
        'mass_multiplier':multipliers[chosen], 'exchange_convexity_limit':limit,
        'mass_multiplier_choices':len(multipliers),
        'vertex_records_checked':sum(v['vertex_records_checked'] for v in results),
        'feasible_vertex_records':sum(v['feasible_vertex_records'] for v in results)}


def contract(box, settings, passes=5, algorithm=1):
    if type(algorithm) is not int or algorithm not in [1, 2, 3, 4]:
        raise ValueError('supported exact exclusion algorithm')
    box, reason, aux = _contract_v1(box, settings, passes)
    if box is None or algorithm == 1:
        return box, reason, aux
    upper = quadratic_target_upper(box, aux[0], settings['radius'])
    if upper is not None and upper['log_Delta_over_target_upper'] < 0:
        upper['center_source'] = 'midpoint'
        upper['legacy_midpoint_Taylor_defined'] = target_log_upper(box, aux[0], settings['radius']) is not None
        return None, 'mass_quadratic_log', upper
    midpoint = [(a+b)//2 for a, b in box]
    nearest = sorted(enumerate(REFERENCES['centers']),
        key=lambda item: sum((a-b)**2 for a, b in zip(midpoint, item[1])))[:2]
    for index, center in nearest:
        upper = quadratic_target_upper(box, aux[0], settings['radius'], center)
        if upper is not None and upper['log_Delta_over_target_upper'] < 0:
            upper['center_source'] = 'reference_'+str(index)
            upper['legacy_midpoint_Taylor_defined'] = target_log_upper(box, aux[0], settings['radius']) is not None
            return None, 'mass_quadratic_log', upper
    if algorithm == 2:
        return box,None,aux
    upper=vertex_target_upper(box,aux[0],settings['radius'])
    if upper is not None and upper['log_Delta_over_target_upper'] < 0:
        upper['center_source']='midpoint'
        upper['legacy_midpoint_Taylor_defined']=target_log_upper(box,aux[0],settings['radius']) is not None
        return None,'mass_quadratic_vertices',upper
    for index,center in nearest:
        upper=vertex_target_upper(box,aux[0],settings['radius'],center)
        if upper is not None and upper['log_Delta_over_target_upper'] < 0:
            upper['center_source']='reference_'+str(index)
            upper['legacy_midpoint_Taylor_defined']=target_log_upper(box,aux[0],settings['radius']) is not None
            return None,'mass_quadratic_vertices',upper
    if algorithm == 3:
        return box,None,aux
    upper=centered_target_upper(box,aux[0],settings['radius'])
    if upper is not None and upper['log_Delta_over_target_upper']<0:
        upper['center_source']='midpoint'
        return None,'mass_centered_quadratic',upper
    for index,center in nearest:
        upper=centered_target_upper(box,aux[0],settings['radius'],center)
        if upper is not None and upper['log_Delta_over_target_upper']<0:
            upper['center_source']='reference_'+str(index)
            return None,'mass_centered_quadratic',upper
    return box, None, aux
