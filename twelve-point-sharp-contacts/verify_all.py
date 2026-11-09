"""Exact contact map and reconstruction for the twelve-point Ward certificate."""
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


sym = load('contact_rank_symmetry', ROOT/'candidate-global-screen/verify_symmetry.py')
prior = sym.load_prior()
Ball, SCALE = prior.Ball, prior.SCALE


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def exact(x):
    if isinstance(x, Ball):
        return x
    q = F(x)
    return Ball.integer(q.numerator)/q.denominator


ZERO, ONE = exact(0), exact(1)


def abs_lower(x):
    if x.lo <= 0 <= x.hi:
        return 0
    return min(abs(x.lo), abs(x.hi))


def ball_power(x, k):
    result = ONE
    for _ in range(k):
        result = result*x
    return result


ZERO_MATRIX = ((ZERO, ZERO, ZERO),)*3


class Jet:
    def __init__(self, value, derivative=(ZERO, ZERO, ZERO), hessian=ZERO_MATRIX):
        self.v = exact(value)
        self.d = derivative
        self.h = hessian

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Jet) else Jet(value)

    def __add__(self, rhs):
        b = self.coerce(rhs)
        hessian = tuple(tuple(x+y for x,y in zip(row_a,row_b))
                        for row_a,row_b in zip(self.h,b.h))
        return Jet(self.v+b.v, tuple(x+y for x,y in zip(self.d,b.d)), hessian)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, tuple(-x for x in self.d),
                   tuple(tuple(-x for x in row) for row in self.h))

    def __sub__(self, rhs):
        return self+(-self.coerce(rhs))

    def __mul__(self, rhs):
        b = self.coerce(rhs)
        hessian = tuple(tuple(self.h[i][j]*b.v+self.v*b.h[i][j]
                              +self.d[i]*b.d[j]+self.d[j]*b.d[i]
                              for j in range(3)) for i in range(3))
        return Jet(self.v*b.v, tuple(x*b.v+y*self.v for x,y in zip(self.d,b.d)), hessian)

    __rmul__ = __mul__

    def __truediv__(self, rhs):
        return self*Jet(ONE/exact(rhs))

    def __pow__(self, k):
        result = Jet(1)
        for _ in range(k):
            result = result*self
        return result


def candidate():
    common = json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    row = next(r for r in common['rows'] if r['n'] == 12)
    inputs = json.loads((ROOT/'candidate-global-screen/symmetry_certificates.json').read_bytes())
    entry = next(r for r in inputs['rows'] if r['n'] == 12)
    certificate = sym.verify_row(entry, prior, common)
    error = F(certificate['root_refinement_error_upper'])
    radius = -((-error.numerator*SCALE)//error.denominator)
    center = [Ball.dyadic(v,common['center_exponent']) for v in row['center_numerators']]
    coordinates = [Ball(x.lo-radius,x.hi+radius) for x in center]
    coordinates.insert(1,ZERO)
    points = list(zip(coordinates[::2], coordinates[1::2]))
    mass_scale = exact(F(11,12))
    us = [mass_scale*(x.square()+y.square()) for x,y in points]
    layers = [v['point_indices'] for v in certificate['radial_layers']]
    need([len(v) for v in layers]==[3,3,6], 'fixed D3 populations')
    layer_index = {i:k for k,layer in enumerate(layers) for i in layer}
    nodes = [us[v[0]] for v in layers]
    need(exact(F(17,100)).hi < nodes[0].lo < nodes[0].hi < exact(F(18,100)).lo, 'inner radial interval')
    need(exact(F(107,100)).hi < nodes[1].lo < nodes[1].hi < exact(F(109,100)).lo, 'middle radial interval')
    anchor = certificate['D3_third_moment_nondegeneracy']['axis_point_index']
    reflected = certificate['D3_third_moment_nondegeneracy']['reflection_fixing_axis']
    w = points[anchor]
    aligned = [(x*w[0]+y*w[1],y*w[0]-x*w[1]) for x,y in points]
    middle_axis = [i for i in layers[1] if reflected[i]==i]
    need(len(middle_axis)==1 and aligned[middle_axis[0]][0].hi<0, 'middle triangle opposite inner axis')
    need(any(aligned[i][0].lo>0 and aligned[i][1].lo>0 for i in layers[2]), 'outer branch c,d positive')

    raw_moments = {0: [exact(12),exact(-1)]+[
        sum((ball_power(u-1,j) for u in us),ZERO) for j in range(2,5)]}
    for k in (3,6,9):
        powers=[]
        for a,b in aligned:
            p,q=ONE,ZERO
            for _ in range(k):
                p,q=p*a-q*b,p*b+q*a
            powers.append(p)
        raw_moments[k]=[sum((p*ball_power(u-1,j) for p,u in zip(powers,us)),ZERO) for j in range(5)]
    need(all(abs_lower(row[0])>0 for row in raw_moments.values()), 'every needed candidate pivot moment is nonzero')
    ratios={k:[v/row[0] for v in row] for k,row in raw_moments.items()}
    means=[mass_scale]+[sum((ball_power(u,j) for u in us),ZERO)/12 for j in range(2,5)]
    contacts=[]
    for obj in certificate['pair_contact_orbits']:
        i,j=obj['representative']
        if layer_index[i]>layer_index[j]:
            i,j=j,i
        u,v=nodes[layer_index[i]],nodes[layer_index[j]]
        a,b=points[i]
        c,d=points[j]
        t=mass_scale*(a*c+b*d)
        gram=u*v-t.square()
        physical_boundary=gram.lo<=0
        if physical_boundary:
            need({layer_index[i],layer_index[j]}=={0,1} and t.hi<0, 'unique antipodal inner-middle contact')
        need(u.lo>0 and v.hi<exact(F(121,72)).lo and (u+v-2*t).lo>exact(F(1,3)).hi, 'strict non-Gram pair-domain margins')
        contacts.append({'u':u,'v':v,'t':t,'same_layer':layer_index[i]==layer_index[j],
                         'boundary':physical_boundary,'pair':(i,j),'weight':obj['multiplicity']})
    need(len(contacts)==14 and sum(c['boundary'] for c in contacts)==1, '14 contacts, one Gram boundary')
    need(sum(c['same_layer'] for c in contacts)==6, 'six same-layer derivative duplications')
    need(sum(c['weight'] for c in contacts)==66, 'all unordered pairs accounted for')
    return contacts,nodes,means,ratios


def contact_columns(contact,nodes,means,ratios):
    u=Jet(contact['u'],(ONE,ZERO,ZERO));v=Jet(contact['v'],(ZERO,ONE,ZERO));t=Jet(contact['t'],(ZERO,ZERO,ONE))
    H=[Jet(1),t]
    for k in range(2,12):
        H.append(2*t*H[-1]-u*v*H[-2])
    up=[u**k for k in range(12)];vp=[v**k for k in range(12)]
    pu=[(u-1)**j for j in range(5)];pv=[(v-1)**j for j in range(5)]
    qu,qv=Jet(1),Jet(1)
    for node in nodes:
        qu=qu*(u-node);qv=qv*(v-node)
    radial_u=[u**j-means[j-1] for j in range(1,5)]+[10*qu**2]
    radial_v=[v**j-means[j-1] for j in range(1,5)]+[10*qv**2]
    names=[];columns=[]
    def add(name,k,cross,du,dv):
        names.append(name)
        columns.append(-H[k]*cross-(up[k]*du+vp[k]*dv)/22)
    for k in range(12):
        if k==0:
            fu,fv=radial_u,radial_v
        elif k in ratios:
            fu=[pu[j]-ratios[k][j] for j in range(1,5)]
            fv=[pv[j]-ratios[k][j] for j in range(1,5)]
        else:
            fu,fv=pu,pv
        for i in range(len(fu)):
            for j in range(i,len(fu)):
                if i==j:
                    add(f'Q{k}_{i}{j}',k,fu[i]*fv[j],fu[i]*fu[j],fv[i]*fv[j])
                else:
                    add(f'Q{k}_{i}{j}',k,fu[i]*fv[j]+fu[j]*fv[i],2*fu[i]*fu[j],2*fv[i]*fv[j])
    for i in (0,1):
        for j in range(i,5):
            if i==j:
                add(f'A1_{i}{j}',1,pu[i]*pv[j],pu[i]*pu[j],pv[i]*pv[j])
            else:
                add(f'A1_{i}{j}',1,pu[i]*pv[j]+pu[j]*pv[i],2*pu[i]*pu[j],2*pv[i]*pv[j])
    for k,constant in ((0,F(1,12)),(2,F(-3,4)),(3,F(-2,3))):
        ward_u,ward_v=pu[1]+constant,pv[1]+constant
        for j in range(0 if k==2 else 1,5):
            a,b=(pu[j],pv[j]) if k==2 else (pu[j]-ratios[k][j],pv[j]-ratios[k][j])
            add(f'A{k}_a{j}',k,ward_u*b+a*ward_v,2*ward_u*a,2*ward_v*b)
    need(len(columns)==187, '165 PSD and 22 first-order compatible Ward directions')
    return names,columns


def build():
    contacts,nodes,means,ratios=candidate()
    names=None;row_names=[];matrix=[]
    for index,contact in enumerate(contacts):
        current,col=contact_columns(contact,nodes,means,ratios)
        if names is None:
            names=current
        else:
            need(names==current,'consistent coefficient ordering')
        matrix.append([x.v for x in col]);row_names.append(f'{index}:value')
        if contact['boundary']:
            su=contact['t']/(2*contact['u']);sv=contact['t']/(2*contact['v'])
            matrix.extend([[x.d[0]+su*x.d[2] for x in col],[x.d[1]+sv*x.d[2] for x in col]])
            row_names.extend([f'{index}:tangent_u',f'{index}:tangent_v'])
        else:
            for j in ((0,2) if contact['same_layer'] else (0,1,2)):
                matrix.append([x.d[j] for x in col]);row_names.append(f'{index}:derivative_{j}')
    need(len(matrix)==49, '14 values and 35 independent-coordinate derivative rows')
    return names,row_names,matrix,contacts,ratios,nodes,means


def log_series(value):
    z = (value-1)/(value+1)
    need(0 <= z.lo <= z.hi < SCALE//2, 'valid logarithm series argument')
    square = z.square()
    term = z
    total = ZERO
    for j in range(100):
        total += 2*term/(2*j+1)
        term *= square
    remainder = 2*term/(201*(ONE-square))
    return Ball(total.lo, total.hi+remainder.hi)


def log_ball(value):
    need(value.lo > 0, 'positive contact distance')
    power = 0
    while value.lo < SCALE:
        value *= 2
        power -= 1
    while value.lo >= 2*SCALE:
        value /= 2
        power += 1
    need(value.hi <= 2*SCALE, 'contact distance stays within its dyadic log range')
    return log_series(value)+power*log_series(exact(2))


def right_hand_side(contacts):
    distances = [c['u']+c['v']-2*c['t'] for c in contacts]
    logs = [log_ball(d) for d in distances]
    B = sum((c['weight']*value for c,value in zip(contacts,logs)), ZERO)
    target = []
    for c,d,value in zip(contacts,distances,logs):
        u,v,t = c['u'],c['v'],c['t']
        h0 = (6*(u+v))/11+B/66-1
        target.append(value-h0)
        radial = ONE/d-exact(F(6,11))
        angular = -2/d
        if c['boundary']:
            target.extend([radial+t*angular/(2*u), radial+t*angular/(2*v)])
        else:
            target.extend([radial,angular] if c['same_layer'] else [radial,radial,angular])
    return target, B


def certify_minor(matrix, rows, columns, residual):
    a = [[matrix[i][j] for j in columns]+[residual[i]] for i in rows]
    size = len(rows)
    determinant = ONE
    lower_pivots = []
    for k in range(size):
        p = max(range(k,size), key=lambda i:abs_lower(a[i][k]))
        need(abs_lower(a[p][k])>0, f'exact pivot {k} excludes zero')
        if p != k:
            a[k],a[p] = a[p],a[k]
            determinant = -determinant
        pivot = a[k][k]
        determinant *= pivot
        lower_pivots.append(abs_lower(pivot))
        for i in range(k+1,size):
            factor = a[i][k]/pivot
            for j in range(k+1,size+1):
                a[i][j] -= factor*a[k][j]
    solution = [ZERO]*size
    for k in reversed(range(size)):
        solution[k] = (a[k][-1]-dot(a[k][k+1:size],solution[k+1:]))/a[k][k]
    return lower_pivots, determinant, solution


def dot(x,y):
    return sum((a*b for a,b in zip(x,y)), ZERO)


def dyadic(value):
    return exact(F.from_float(float.fromhex(value)))


def seed_parameters(seed, names):
    values = {}
    for k, block in enumerate(seed['Q_blocks']):
        Q = [list(map(dyadic,row)) for row in block]
        if k:
            E = [list(map(dyadic,row)) for row in seed['transforms'][k-1]]
            if k in (3,6,9):
                E = E[1:]
            Q = [[sum((E[i][a]*Q[a][b]*E[j][b] for a in range(len(Q))
                       for b in range(len(Q))), ZERO) for j in range(len(E))]
                 for i in range(len(E))]
        for i in range(len(Q)):
            for j in range(i,len(Q)):
                values[f'Q{k}_{i}{j}'] = Q[i][j]
    for name in names:
        if name.startswith('A'):
            prefix = {'0':'mass','2':'a2','3':'a3'}.get(name[1])
            source_name = name if prefix is None else f'{prefix}_{name[-1]}'
            values[name] = dyadic(seed['ward_parameters'][source_name])
    return [values[name] for name in names]


def ldl_lower(matrix, label):
    dim = len(matrix)
    L = [[ZERO]*dim for _ in range(dim)]
    D = []
    for i in range(dim):
        diag = matrix[i][i]-sum((L[i][k].square()*D[k] for k in range(i)),ZERO)
        need(diag.lo > 0, label+' is positive definite at pivot '+str(i))
        D.append(diag)
        for j in range(i+1,dim):
            L[j][i] = (matrix[j][i]-sum((L[j][k]*L[i][k]*D[k]
                                       for k in range(i)),ZERO))/diag
    return min(x.lo for x in D)


def check_grams(values):
    pivots = []
    for k in range(12):
        dim = 4 if k in (3,6,9) else 5
        block = [[values[f'Q{k}_{min(i,j)}{max(i,j)}'] for j in range(dim)]
                 for i in range(dim)]
        pivots.append(ldl_lower(block, f'corrected Gram block {k}'))
    return pivots


def check_contacts(contacts, nodes, means, ratios, names, coefficients):
    pivots = []
    normal_margin = None
    unit = [1,1,-2]
    for idx,c in enumerate(contacts):
        column_names, jets = contact_columns(c,nodes,means,ratios)
        need(names == column_names, 'second derivative uses identical basis')
        u,v,t = c['u'],c['v'],c['t']
        d = u+v-2*t
        inv_d = ONE/d
        H = [[dot(coefficients,[jet.h[i][j] for jet in jets])+unit[i]*unit[j]*inv_d.square()
              for j in range(3)] for i in range(3)]
        if c['boundary']:
            Gt = dot(coefficients,[jet.d[2] for jet in jets])+2*inv_d
            normal_margin = -(t*Gt).hi
            need(normal_margin > 0, 'strict positive inward derivative at the Gram boundary')
            tangents = [[ONE,ZERO,t/(2*u)],[ZERO,ONE,t/(2*v)]]
            mixed = t/(4*u*v)
            second_t = [[-t/(4*u.square()),mixed],[mixed,-t/(4*v.square())]]
            H = [[sum((tangents[i][a]*H[a][b]*tangents[j][b]
                       for a in range(3) for b in range(3)),ZERO)+Gt*second_t[i][j]
                  for j in range(2)] for i in range(2)]
        pivots.append(ldl_lower(H, f'contact Hessian {idx}'))
    need(normal_margin is not None, 'one boundary contact checked')
    return pivots, normal_margin


def check_negative_pair(B,nodes,means,ratios,names,coefficients):
    # A rational counterexample to the global majorant for this exact fit.
    u,v,t = F(941,1000),F(1419,1000),F(-979,1000)
    d = u+v-2*t
    need(F(47,50) < u < v < F(71,50) and t*t < u*v and d > F(1,3),
         'negative pair strictly inside the outer-outer physical domain')
    need(11*(u+v)+2*t < F(3427,100), 'negative pair also passes the deletion linear cut')
    c = dict(u=exact(u),v=exact(v),t=exact(t))
    column_names,jets = contact_columns(c,nodes,means,ratios)
    need(names == column_names, 'negative pair uses identical basis')
    gap = B/66-1+exact(F(6,11)*(u+v))+dot(coefficients,[jet.v for jet in jets])-log_ball(exact(d))
    need(gap.hi < exact(F(-1,4000)).lo, 'strict negative gap at the rational pair')
    return [str(x) for x in (u,v,t)]


def lower_power(number):
    return number.bit_length()-1-prior.PRECISION


def verify():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for path,sha in manifest.items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha, 'fixed input '+path)
    names,row_names,matrix,contacts,ratios,nodes,means = build()
    minor = json.loads((HERE/'minor.json').read_bytes())
    rows,columns = minor['selected_rows'],minor['selected_columns']
    need(len(set(rows)) == len(set(columns)) == 44, '44 distinct rows and columns')
    need([row_names[i] for i in rows] == minor['selected_row_names'] and
         [names[j] for j in columns] == minor['selected_column_names'],
         'fixed minor matches the contact ordering')
    seed = json.loads((HERE/'seed.json').read_bytes())
    parameters = seed_parameters(seed,names)
    target,B = right_hand_side(contacts)
    residual = [b-dot(row,parameters) for row,b in zip(matrix,target)]
    pivots,determinant,correction = certify_minor(matrix,rows,columns,residual)
    need(abs_lower(determinant)>0, 'interval determinant excludes zero')
    for index,increment in zip(columns,correction):
        parameters[index] += increment
    for row,b in zip(matrix,target):
        enclosure = dot(row,parameters)-b
        need(enclosure.lo <= 0 <= enclosure.hi, 'full contact row is interval-consistent')
    grams = check_grams(dict(zip(names,parameters)))
    local, normal = check_contacts(contacts,nodes,means,ratios,names,parameters)
    negative = check_negative_pair(B,nodes,means,ratios,names,parameters)
    return {
        'passed':True, 'n':12, 'contact_classes':len(contacts),
        'interior_contacts':13, 'gram_boundary_contacts':1,
        'linear_map_shape':[len(matrix),len(names)], 'certified_nonzero_minor_order':44,
        'contact_minor_min_pivot_strict_lower_power_of_two':lower_power(min(pivots))-1,
        'contact_minor_determinant_excludes_zero':True,
        'all_contact_values_and_first_derivatives_reconstructed':True,
        'exact_positive_definite_gram_blocks':len(grams),
        'minimum_gram_pivot_strict_lower_power_of_two':lower_power(min(grams))-1,
        'contacts_with_positive_definite_interior_or_tangent_hessian':len(local),
        'minimum_contact_pivot_strict_lower_power_of_two':lower_power(min(local))-1,
        'strict_inward_derivative_lower_power_of_two':lower_power(normal)-1,
        'rational_negative_pair':negative, 'negative_pair_gap_strict_upper':'-1/4000',
        'continuous_complement_covered':False, 'sharp_global_value_proved':False,
        'input_sha256':manifest,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args = parser.parse_args()
    report = verify()
    result = HERE/'verification_report.json'
    rendered = json.dumps(report,indent=2)+'\n'
    if args.record_report:
        result.write_bytes(rendered.encode('utf-8'))
    else:
        need(json.loads(result.read_bytes()) == report,'fresh replay agrees with saved report')
    print(rendered,end='')


if __name__=='__main__':
    main()


