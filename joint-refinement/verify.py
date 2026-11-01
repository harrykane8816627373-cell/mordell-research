"""Exact checks of rotation saturation, logarithm arithmetic and continuation guards."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import gzip
import json
import sys
import tempfile
import refinement as R
from check import check

E = R.E
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'power-moment'))
from complex_interval import C


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def force(q):
    return [sum((1/(z-w) for j, w in enumerate(q) if i != j), C()) for i, z in enumerate(q)]


def H_action(q, h):
    return [4*h[i]+sum(((h[i]-h[j])/(q[i]-q[j])**2
        for j in range(len(q)) if i != j), C()).conjugate() for i in range(len(q))]


def inner(a, b):
    return sum((x.conjugate()*y for x, y in zip(a, b)), C()).re


def rotation_checks():
    q0 = [C(1), C(-1), C(0, 1), C(0, -1), C(1, F(1, 2)),
          C(-1, -F(1, 2)), C(F(1, 2)), C(-F(1, 2))]
    q1 = list(q0)
    q1[0], q1[1] = F(3, 5)*q0[0]-F(4, 5)*q0[1], F(4, 5)*q0[0]+F(3, 5)*q0[1]
    rows = []
    for q in [q0, q1]:
        mass = sum(z.norm() for z in q)
        need(mass == 7 and len(set(q)) == 8, 'actual distinct mass-seven rational points')
        f = force(q)
        g = [4*z-v.conjugate() for z, v in zip(q, f)]
        u, v = [C(0, 1)*z for z in q], [C(0, 1)*z for z in g]
        Hu = H_action(q, u)
        need(Hu == v and inner(u, Hu) == 0, 'rotation row has zero quadratic form and equals rotated gradient')
        residual2 = sum(z.norm() for z in g)
        need(residual2 > 0 and sum(z.norm() for z in f) == 112+residual2,
             'inverse-square norm lower bound with exact force-residual square')
        row_budget = max(sum((F(1, (q[i]-q[j]).norm()) for j in range(8) if i != j), F())
                         for i in range(8))
        bound = 4+2*row_budget
        epsilon = 1/bound
        h = [a-epsilon*b for a, b in zip(u, v)]
        quadratic = inner(h, H_action(q, h))
        Hv = inner(v, H_action(q, v))
        need(quadratic == -2*epsilon*residual2+epsilon**2*Hv < 0,
             'explicit rational negative-curvature direction at a nonstationary normalized point')
        rows.append({'centered': sum(q, C()) == C(), 'residual_squared': str(residual2),
            'negative_quadratic': str(quadratic)})
    far = [100*z for z in q0]
    row_bound = 2*max(sum((F(1, (far[i]-far[j]).norm()) for j in range(8) if i != j), F())
                      for i in range(8))
    need(row_bound < 4 and sum(z.norm() for z in far) == 70000,
         'wrong-mass actual configuration has strictly positive full Hessian')
    return {'two_actual_normalized_witnesses': rows,
        'centering_not_required_for_rotation_theorem': True,
        'wrong_mass_strict_positive_Hessian_witness_mass': 70000,
        'wrong_mass_inverse_square_operator_bound': str(row_bound)}


def reference_log(a, b):
    ratio, exponent = F(a, b), 0
    while ratio < 1:
        ratio *= 2
        exponent -= 1
    while ratio >= 2:
        ratio /= 2
        exponent += 1

    def reduced(z):
        x = (z-1)/(z+1)
        m = 64
        total = 2*sum((x**(2*k+1)/F(2*k+1) for k in range(m)), F())
        tail = 2*x**(2*m+1)/((2*m+1)*(1-x*x))
        return total, total+tail
    lo, hi = reduced(ratio)
    l2, u2 = reduced(F(2))
    return (lo+exponent*l2, hi+exponent*u2) if exponent >= 0 else (lo+exponent*u2, hi+exponent*l2)


def logarithm_checks():
    E.setup(40)
    pairs = [(1, 1), (2, 1), (1, 2), (7, 10), (11, 7), (100, 1),
             (1, 10**100), (10**100, 1), (2**127-1, 2**100), (7**7, 2**30)]
    for a, b in pairs:
        lo, hi = R.log_ratio(a, b)
        reflo, refhi = reference_log(a, b)
        need(F(lo, E.S) <= reflo <= refhi <= F(hi, E.S),
             'integer logarithm encloses independent high-precision Fraction series')
    return {'independent_fraction_logarithm_cases': len(pairs), 'floating_point_acceptance': False}


def fixed(v):
    return E.rat(v.numerator, v.denominator)


def exact_log_product(q):
    product = F(1)
    for i, j in E.PAIRS:
        product *= (q[i]-q[j]).norm()
    return reference_log(product.numerator, product.denominator*7**7)


def Taylor_checks():
    E.setup(40)
    q = [C(1), C(-1), C(0, 1), C(0, -1), C(1, F(1, 2)),
         C(-1, -F(1, 2)), C(F(1, 2)), C(-F(1, 2))]
    rotated = list(q)
    n = 10001
    phase = C(F(n*n-1, n*n+1), F(2*n, n*n+1))
    rotated[4], rotated[5] = q[4]*phase, q[5]*phase
    need(sum(z.norm() for z in rotated) == 7 and sum(rotated, C()) == C(),
         'second actual point has exact center and mass')
    data = [[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]],
            [rotated[0].re]+[v for z in rotated[1:7] for v in [z.re,z.im]]]
    box = [(min(fixed(a)[0],fixed(b)[0]), max(fixed(a)[1],fixed(b)[1])) for a,b in zip(*data)]
    center = [int(a*E.S) for a in data[0]]
    # A second, deliberately outside center tests the reference-center form.
    outside = [a+(E.S//1024 if i%3 == 0 else -E.S//2048) for i,a in enumerate(center)]
    bounds = []
    for c in [center, outside]:
        result = R.target_log_upper(box, E.positions(box), E.model()['radius'], c)
        need(result is not None and result['log_Delta_over_target_upper'] < 0,
             'mass-adjusted Taylor bound excludes this entire normalized box')
        for actual in [q, rotated]:
            lower, _ = exact_log_product(actual)
            need(lower <= F(result['log_Delta_over_target_upper'], E.S),
                 'direct exact distance product is bounded by the Taylor certificate')
        bounds.append(str(F(result['log_Delta_over_target_upper'], E.S)))
    wide = [(a-10*E.S,a+10*E.S) for a in center]
    need(R.target_log_upper(wide, E.positions(wide), E.model()['radius'], center) is None,
         'a Taylor disk reaching a pole is never treated as an exclusion')
    return {'actual_mass_seven_points_checked': 2, 'arbitrary_reference_centers_checked': 2,
        'whole_box_log_upper_bounds': bounds, 'wide_pole_box_rejected_as_inconclusive': True}


def quadratic_checks():
    E.setup(40)
    cases = 0
    for d in [-4*E.S,-E.S//3,0,E.S//3,4*E.S]:
        for b in [(-2*E.S,E.S),(-E.S,-E.S),(0,0),(E.S,2*E.S),(-E.S//7,E.S//5)]:
            for interval in [(-3*E.S,-E.S),(-E.S,E.S),(E.S,2*E.S),(0,0),(0,2*E.S),(-2*E.S,0)]:
                actual = []
                for slope in map(lambda x:F(x,E.S),b):
                    coefficient = F(d,E.S)
                    lo, hi = map(lambda x:F(x,E.S),interval)
                    points = [lo,hi]
                    if coefficient < 0:
                        vertex = -slope/(2*coefficient)
                        if lo <= vertex <= hi:
                            points.append(vertex)
                    actual += [coefficient*x*x+slope*x for x in points]
                exact = max(actual)
                upper = F(R._quadratic_max(d,b,interval),E.S)
                need(exact <= upper < exact+F(1,E.S), 'independent Fraction quadratic endpoint/vertex maximum')
                cases += 1
    q = [C(1),C(-1),C(0,1),C(0,-1),C(1,F(1,2)),C(-1,-F(1,2)),C(F(1,2)),C(-F(1,2))]
    displacement = [C(F(i-3,11),F(2*i-5,13)) for i in range(8)]
    a = [(i+1)**2 for i in range(8)]
    diagonal = [F(-4)]*8
    lhs = -4*sum(z.norm() for z in displacement)
    square_remainder = F()
    for i,j in E.PAIRS:
        w = 1/(q[i]-q[j]).norm()
        diagonal[i] += w*(1+F(a[j],a[i]))
        diagonal[j] += w*(1+F(a[i],a[j]))
        lhs += w*(displacement[i]-displacement[j]).norm()
        square_remainder += w*(F(j+1,i+1)*displacement[i]+F(i+1,j+1)*displacement[j]).norm()
    rhs = sum(d*z.norm() for d,z in zip(diagonal,displacement))
    need(rhs-lhs == square_remainder >= 0, 'exact weighted-Laplacian diagonal majorizer identity')
    values = [q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]]
    center = [int(v*E.S) for v in values]
    point_box = [E.point(v) for v in center]
    lower,_ = exact_log_product(q)
    for scalar in [F(-1),F(1,2),F(2)]:
        arbitrary = [int(v*scalar) for v in center]
        upper = R.quadratic_target_upper(point_box,E.positions(point_box),E.model()['radius'],arbitrary)
        need(upper is not None and lower <= F(upper['log_Delta_over_target_upper'],E.S),
             'global majorant accepts arbitrary distant rational reference centers')
    line = [C(v) for v in [F(6,5),F(-6,5),F(9,10),F(-9,10),F(1),F(-1),F(1,2),F(-1,2)]]
    need(sum(line,C()) == C() and sum(z.norm() for z in line) == 7 and len(set(line)) == 8,
         'wide-box benchmark has an actual distinct centered mass-seven point')
    values = [line[0].re]+[v for z in line[1:7] for v in [z.re,z.im]]
    center = [round(v*E.S) for v in values]
    radius = E.S//16
    box = [(v-radius,v+radius) for v in center]
    c = E.positions([E.point(v) for v in center])
    qbox = [((z[0][0]-radius,z[0][0]+radius),(z[1][0]-radius,z[1][0]+radius)) for z in c]
    qbox[0] = box[0],E.ZERO
    for actual,interval in zip(line,qbox):
        need(all(F(interval[k][0],E.S) <= value <= F(interval[k][1],E.S)
            for k,value in enumerate([actual.re,actual.im])), 'actual benchmark is inside the whole physical box')
    need(R.target_log_upper(box,qbox,E.model()['radius'],center) is None,
         'legacy Taylor disk test is unavailable on this benchmark')
    upper = R.quadratic_target_upper(box,qbox,E.model()['radius'],center)
    need(upper is not None and upper['log_Delta_over_target_upper'] < 0,
         'global quadratic bound excludes entire normalized wide box without a Taylor disk')
    direct,_ = exact_log_product(line)
    need(direct <= F(upper['log_Delta_over_target_upper'],E.S), 'independent exact product of the physical benchmark')
    return {'independent_scalar_quadratic_maxima':cases,'weighted_Laplacian_square_identity':True,
        'arbitrary_distant_reference_centers':3,'nonempty_actual_wide_box_radius':str(F(radius,E.S)),
        'legacy_Taylor_unavailable_on_wide_box':True,
        'wide_box_log_upper':str(F(upper['log_Delta_over_target_upper'],E.S)),
        'wide_box_is_claimed_to_pass_all_old_global_cuts':False,'floating_point_acceptance':False}


def vertex_checks():
    E.setup(40)
    cases,vertices=0,0
    for index in range(6):
        weights=[(i+j+index+1)*E.S//(100 if index==0 else 7+index)
            for i,j in E.PAIRS]
        intervals=[((-2+(i%3))*E.S//7,(1+(i%4))*E.S//9) for i in range(8)]
        if index%2:
            intervals[0]=(E.S//3,E.S//3)
        slopes=[(-((i+1)*E.S)//13,(i+2)*E.S//17) for i in range(8)]
        upper,count,clipped=R.vertex_quadratic_upper(intervals,slopes,weights)
        rows=[-4*E.S]*8
        for (i,j),w in zip(E.PAIRS,weights):
            rows[i]+=w
            rows[j]+=w
        best=None
        for x in product(*[(a,) if a==b else (a,b) for a,b in intervals]):
            raw=sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
            raw-=4*E.S*sum(v*v for v in x)
            raw+=E.S*sum((b[1] if v>=0 else b[0])*v for b,v in zip(slopes,x))
            for d,v,(lo,hi) in zip(rows,x,intervals):
                if d<0:
                    minimum=0 if lo<=0<=hi else min(lo*lo,hi*hi)
                    raw+=(-d)*(v*v-minimum)
            best=raw if best is None else max(best,raw)
            vertices+=1
        need(upper==-(-best//(E.S*E.S)), 'Gray-code maximum equals independent full Cartesian enumeration')
        need(count==2**sum(lo!=hi for lo,hi in intervals)
             and clipped==sum(v<0 for v in rows), 'vertex and clipping counts')
        cases+=1
    line=[C(v) for v in [F(6,5),F(-6,5),F(9,10),F(-9,10),F(1),F(-1),F(1,2),F(-1,2)]]
    values=[line[0].re]+[v for z in line[1:7] for v in [z.re,z.im]]
    center=[round(v*E.S) for v in values]
    r=E.S//16
    box=[(v-r,v+r) for v in center]
    c=E.positions([E.point(v) for v in center])
    qbox=[((z[0][0]-r,z[0][0]+r),(z[1][0]-r,z[1][0]+r)) for z in c]
    qbox[0]=box[0],E.ZERO
    bound=R.vertex_target_upper(box,qbox,E.model()['radius'],center)
    lower,_=exact_log_product(line)
    need(bound is not None and lower <= F(bound['log_Delta_over_target_upper'],E.S) < 0,
         'physical wide-box target bound by full quadratic vertex majorant')
    return {'independent_Cartesian_enumeration_cases':cases,'independently_enumerated_vertices':vertices,
        'wide_box_log_upper':str(F(bound['log_Delta_over_target_upper'],E.S)),
        'maximum_vertices_per_physical_box':512,'floating_point_acceptance':False}


def centered_checks():
    E.setup(40)
    cases,feasible_records=0,0
    for index in range(8):
        weights=[(2*i+j+index+1)*E.S//(11+index) for i,j in E.PAIRS]
        intervals=[((-3+(i%3))*E.S//5,(1+(i%4))*E.S//7) for i in range(8)]
        if index%2:
            intervals[0]=(0,0)
        if index==6:
            intervals=[(v*E.S,v*E.S) for v in [-1,1,0,0,0,0,0,0]]
        if index==7:
            intervals=[(-E.S,E.S)]+[(0,0)]*7
        rows=[0]*8
        matrix=[[0]*8 for _ in range(8)]
        for (i,j),w in zip(E.PAIRS,weights):
            rows[i]+=w
            rows[j]+=w
            matrix[i][j]=matrix[j][i]=w
        active=[i for i,(lo,hi) in enumerate(intervals) if lo!=hi]
        limit=min((rows[i]+rows[j]+2*matrix[i][j])//2
            for i in active for j in active if i<j) if len(active)>1 else 4*E.S
        etas=sorted(set([-E.S,0,min(4*E.S,limit),limit]))
        actual=R.centered_slice_max(intervals,weights,etas)
        expected=[None]*len(etas)
        count=0
        free_choices=active or [None]
        for free in free_choices:
            indices=[i for i in range(8) if i!=free]
            choices=[(intervals[i][0],) if intervals[i][0]==intervals[i][1] else intervals[i] for i in indices]
            for values in product(*choices):
                x=[0]*8
                for i,v in zip(indices,values):
                    x[i]=v
                if free is not None:
                    x[free]=-sum(values)
                if any(not lo<=v<=hi for v,(lo,hi) in zip(x,intervals)) or sum(x)!=0:
                    continue
                count+=1
                graph=sum(w*(x[i]-x[j])**2 for (i,j),w in zip(E.PAIRS,weights))
                mass=sum(v*v for v in x)
                for k,eta in enumerate(etas):
                    value=graph-eta*mass
                    expected[k]=value if expected[k] is None else max(expected[k],value)
        need(actual is not None and actual['upper']==[-(-v//(E.S*E.S)) for v in expected],
             'centered Gray-code updates equal independent sliced Cartesian enumeration')
        need(actual['feasible_vertex_records']==count,'all feasible sliced vertex records accounted for')
        feasible_records+=count
        cases+=1
    need(R.centered_slice_max([(E.S,2*E.S)]*8,[E.S]*28,[0]) is None,
         'an empty centroid slice never produces a numeric target bound')
    try:
        R.centered_slice_max([(-E.S,E.S)]*2+[(0,0)]*6,[E.S]*28,[9*E.S])
    except ValueError:
        invalid_rejected=True
    else:
        raise ArithmeticError('unsafe exchange multiplier accepted')
    # At eta=9 for the unit complete graph, Q(t,-t,0,...,0)=-2*t^2.
    # The midpoint value 0 exceeds the two boundary values -2.
    def unsafe_value(t):
        x=[F(t),-F(t)]+[F(0)]*6
        return sum((x[i]-x[j])**2 for i,j in E.PAIRS)-9*sum(v*v for v in x)
    need(unsafe_value(0)==0 and unsafe_value(1)==unsafe_value(-1)==-2
        and unsafe_value(0)>max(unsafe_value(1),unsafe_value(-1)),
        'interior counterexample if exchange convexity were omitted')
    q=[C(1),C(-1),C(0,1),C(0,-1),C(1,F(1,2)),C(-1,-F(1,2)),C(F(1,2)),C(-F(1,2))]
    values=[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]]
    center=[int(v*E.S) for v in values]
    box=[E.point(v) for v in center]
    qbox=E.positions(box)
    bound=R.centered_target_upper(box,qbox,E.model()['radius'],center)
    lower,upper=exact_log_product(q)
    need(bound is not None and lower<=F(bound['log_Delta_over_target_upper'],E.S)<upper+F(1,10**6),
         'mass-seven physical point has a tight direct-log bound at its own reference center')
    c=[2*z for z in q]
    graph=sum((q[i]-q[j]).norm()/(c[i]-c[j]).norm() for i,j in E.PAIRS)
    for eta in [F(-3,2),F(0),F(4),F(17,3)]:
        need(7*eta+graph-eta*sum(z.norm() for z in q)==graph,
             'arbitrary exact mass multiplier leaves the physical global bound unchanged')
    line=[C(v) for v in [F(6,5),F(-6,5),F(9,10),F(-9,10),F(1),F(-1),F(1,2),F(-1,2)]]
    values=[line[0].re]+[v for z in line[1:7] for v in [z.re,z.im]]
    center=[round(v*E.S) for v in values]
    radius=E.S//12
    box=[(v-radius,v+radius) for v in center]
    cp=E.positions([E.point(v) for v in center])
    qbox=[((z[0][0]-radius,z[0][0]+radius),(z[1][0]-radius,z[1][0]+radius)) for z in cp]
    qbox[0]=box[0],E.ZERO
    centered=R.centered_target_upper(box,qbox,E.model()['radius'],center)
    old=R.vertex_target_upper(box,qbox,E.model()['radius'],center)
    need(centered is not None and old is not None
        and centered['log_Delta_over_target_upper']<old['log_Delta_over_target_upper']<0,
        'physical wide box improves after retaining centroid and mass multiplier')
    direct,_=exact_log_product(line)
    need(direct<=F(centered['log_Delta_over_target_upper'],E.S),'independent exact product in the centered wide box')
    return {'sliced_Cartesian_enumeration_cases':cases,'independent_feasible_vertex_records':feasible_records,
        'unsafe_exchange_multiplier_rejected':invalid_rejected,
        'interior_counterexample_without_exchange_condition':True,'mass_multiplier_identity_cases':4,
        'nonempty_actual_wide_box_radius':str(F(radius,E.S)),
        'centered_wide_box_log_upper':str(F(centered['log_Delta_over_target_upper'],E.S)),
        'previous_vertex_wide_box_log_upper':str(F(old['log_Delta_over_target_upper'],E.S)),
        'floating_point_acceptance':False}


def guard_checks():
    base = json.loads(gzip.decompress((ROOT/'joint-closure/certificates/partial_cover.json.gz').read_bytes()))
    nodes = [{'kind':'unknown'} for _ in base['pending']]
    toy = {'base_sha256':R.BASE_SHA256,'precision_bits':40,'settings':base['settings'],
        'roots':[[i,j] for j,(i,b) in enumerate(base['pending'])], 'nodes':nodes,
        'pending':[[j,b] for j,(i,b) in enumerate(base['pending'])],
        'complete_cover':False,'full_eight_point_exclusion_proved':False}
    rejected=0
    with tempfile.TemporaryDirectory(prefix='refinement-guard-',dir=ROOT) as temp:
        path=Path(temp)/'toy.gz'

        def save(data):
            path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        save(toy)
        need(check(path,False)['unresolved_regions']==21,'initial frontier remains unresolved')
        for kind in range(8):
            altered=json.loads(json.dumps(toy))
            if kind==0:
                altered['roots'].pop()
            elif kind==1:
                altered['full_eight_point_exclusion_proved']=True
            elif kind==2:
                altered['pending'].pop()
            elif kind==3:
                altered['pending'].pop(0)
                altered['nodes'][0]={'kind':'excluded','reason':'mass_log_Taylor'}
            elif kind in [4,5]:
                altered['pending'].pop(0)
                altered['nodes'][0]={'kind':'excluded','reason':'mass_quadratic_log',
                    'algorithm':99 if kind==4 else 2}
            elif kind==6:
                altered['pending'].pop(0)
                altered['nodes'][0]={'kind':'excluded','reason':'mass_quadratic_vertices','algorithm':3}
            else:
                altered['pending'].pop(0)
                altered['nodes'][0]={'kind':'excluded','reason':'mass_centered_quadratic','algorithm':4}
            save(altered)
            try:
                check(path,False)
            except ArithmeticError:
                rejected+=1
            else:
                raise ArithmeticError('invalid frontier accounting or unjustified Taylor exclusion accepted')
    need(rejected==8,'all eight false-proof variants rejected')
    return {'initial_21_regions_preserved':True,'false_proof_variants_rejected':rejected}


if __name__=='__main__':
    print(json.dumps({'passed':True,'rotation_saturation':rotation_checks(),
        'logarithm_arithmetic':logarithm_checks(),'Taylor_arithmetic':Taylor_checks(),
        'global_quadratic_arithmetic':quadratic_checks(),
        'quadratic_vertex_arithmetic':vertex_checks(),
        'centered_mass_slice_arithmetic':centered_checks(),
        'coverage_guards':guard_checks(),'full_eight_point_exclusion_proved':False},indent=2))
