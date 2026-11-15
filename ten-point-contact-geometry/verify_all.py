"""Replay the algebraic D2 representation and sharp contact rank for the n=10 candidate."""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path
from check_root_box import check as check_reduced_box

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
spec = importlib.util.spec_from_file_location('n12_reference',ROOT/'twelve-point-sharp-contacts/verify_all.py')
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
Ball, SCALE = ref.Ball, ref.SCALE
ZERO, ONE, exact = ref.ZERO, ref.ONE, ref.exact
abs_lower, power, need, Jet = ref.abs_lower, ref.ball_power, ref.need, ref.Jet


def candidate():
    common = json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    row = next(r for r in common['rows'] if r['n'] == 10)
    inputs = json.loads((ROOT/'candidate-global-screen/symmetry_certificates.json').read_bytes())
    entry = next(r for r in inputs['rows'] if r['n'] == 10)
    cert = ref.sym.verify_row(entry, ref.prior, common)
    error = F(cert['root_refinement_error_upper'])
    radius = -((-error.numerator*SCALE)//error.denominator)
    center = [Ball.dyadic(v,common['center_exponent']) for v in row['center_numerators']]
    coordinates = [Ball(x.lo-radius,x.hi+radius) for x in center]
    coordinates.insert(1,ZERO)
    points = list(zip(coordinates[::2], coordinates[1::2]))
    scale = exact(F(9,10))
    us = [scale*(x.square()+y.square()) for x,y in points]
    layers = [v['point_indices'] for v in cert['radial_layers']]
    need([len(v) for v in layers]==[2,4,4], 'D2 layer sizes')
    layer_of = {i:k for k,layer in enumerate(layers) for i in layer}
    nodes = [us[v[0]] for v in layers]
    anchor = cert['D2_rectangle_nondegeneracy']['inner_axis_point_index']
    w = points[anchor]
    aligned = [(x*w[0]+y*w[1],y*w[0]-x*w[1]) for x,y in points]
    rr = [scale*(x.square()-y.square())/(w[0].square()+w[1].square())
          for x,y in aligned]
    means = [scale]+[sum((power(u,j) for u in us),ZERO)/10 for j in range(2,5)]
    moments = {0:[exact(10),exact(-1)]+[sum((power(u-1,j) for u in us),ZERO) for j in range(2,5)]}
    for k in (2,4,6,8,10):
        powers = []
        for x,y in aligned:
            p,q = ONE,ZERO
            for _ in range(k):
                p,q=p*x-q*y,p*y+q*x
            powers.append(p)
        moments[k] = [sum((p*power(u-1,j) for p,u in zip(powers,us)),ZERO) for j in range(5)]
    need(all(abs_lower(m[0])>0 for m in moments.values()), 'nonzero raw harmonic pivots')
    ratios = {k:[e/m[0] for e in m] for k,m in moments.items()}
    contacts = []
    for obj in cert['pair_contact_orbits']:
        i,j = obj['representative']
        if layer_of[i]>layer_of[j]:
            i,j=j,i
        u,v = nodes[layer_of[i]],nodes[layer_of[j]]
        pi,pj = points[i],points[j]
        t = scale*(pi[0]*pj[0]+pi[1]*pj[1])
        boundary = j==entry['rotation_permutation'][i]
        same_layer = layer_of[i]==layer_of[j]
        if boundary:
            need(same_layer, 'antipodes on same layer')
            t = -u
        else:
            need((u*v-t.square()).lo>0, 'all nonantipodal contacts in Gram interior')
        d = u+v-2*t
        need(d.lo>exact(F(2,5)).hi and v.hi<exact(F(15403,10000)).lo, 'strict simple domain margins')
        q = 9-exact(F(9,8))*(u+v)-t/4
        need(q.lo>exact(6).hi,'q>6 so no nonlinear Schwarz cut is active at candidate')
        contacts.append(dict(u=u,v=v,t=t,d=d,q=q,boundary=boundary,same_layer=same_layer,
                             layers=(layer_of[i],layer_of[j]),pair=(i,j),weight=obj['multiplicity']))
    need(sum(c['weight'] for c in contacts)==45, 'all pairs counted')
    need(len(contacts)==15 and sum(c['boundary'] for c in contacts)==3 and
         sum(c['same_layer'] for c in contacts)==7, 'contact structure')
    return nodes, [rr[layer[0]] for layer in layers[1:]], means, ratios, contacts


def columns(contact,nodes,means,ratios,K,add_ward=True,degree=4):
    u=Jet(contact['u'],(ONE,ZERO,ZERO))
    v=Jet(contact['v'],(ZERO,ONE,ZERO))
    t=Jet(contact['t'],(ZERO,ZERO,ONE))
    H=[Jet(1),t]
    for k in range(2,K+1):
        H.append(2*t*H[-1]-u*v*H[-2])
    up,vp=([a**k for k in range(K+1)] for a in (u,v))
    pu,pv=([(a-1)**j for j in range(degree+1)] for a in (u,v))
    qu,qv=Jet(1),Jet(1)
    for node in nodes:
        qu,qv=qu*(u-node),qv*(v-node)
    radial_u=[u**j-means[j-1] for j in range(1,degree+1)]+[10*qu**2]
    radial_v=[v**j-means[j-1] for j in range(1,degree+1)]+[10*qv**2]
    names, result = [], []
    def add(name,k,cross,du,dv):
        names.append(name)
        result.append(-H[k]*cross-(up[k]*du+vp[k]*dv)/18)
    def sym_product(name,k,i,j,fu,fv):
        if i==j:
            add(name,k,fu[i]*fv[i],fu[i]*fu[i],fv[i]*fv[i])
        else:
            add(name,k,fu[i]*fv[j]+fu[j]*fv[i],2*fu[i]*fu[j],2*fv[i]*fv[j])
    for k in range(K+1):
        if k==0:
            fu,fv=radial_u,radial_v
        elif k in ratios:
            fu=[pu[j]-ratios[k][j] for j in range(1,degree+1)]
            fv=[pv[j]-ratios[k][j] for j in range(1,degree+1)]
        else:
            fu,fv=pu,pv
        for i in range(len(fu)):
            for j in range(i,len(fu)):
                sym_product(f'Q{k}_{i}{j}',k,i,j,fu,fv)
    if add_ward:
        for i in (0,1):
            for j in range(i,degree+1):
                sym_product(f'A1_{i}{j}',1,i,j,pu,pv)
        for k,c in ((0,F(1,10)),(2,F(-7,10)),(3,F(-3,5))):
            ward_u,ward_v=pu[1]+c,pv[1]+c
            for j in range(1 if k in ratios else 0,degree+1):
                a,b=pu[j],pv[j]
                if k in ratios:
                    a,b=a-ratios[k][j],b-ratios[k][j]
                add(f'A{k}_a{j}',k,ward_u*b+a*ward_v,2*ward_u*a,2*ward_v*b)
    return names,result


def build(nodes,means,ratios,contacts,K,add_ward=True,degree=4):
    names=None
    row_names,matrix=[],[]
    for i, c in enumerate(contacts):
        new, jets = columns(c,nodes,means,ratios,K,add_ward,degree)
        need(names is None or names==new, 'column order')
        names=new
        row_names.append(f'{i}:value')
        matrix.append([j.v for j in jets])
        if c['boundary']:
            row_names.append(f'{i}:tangent')
            matrix.append([j.d[0]-j.d[2]/2 for j in jets])
        else:
            for d in ((0,2) if c['same_layer'] else (0,1,2)):
                row_names.append(f'{i}:derivative_{d}')
                matrix.append([j.d[d] for j in jets])
    need(len(matrix)==50,'50 contact equations')
    return names,row_names,matrix


def lower_rank(matrix):
    m=[r.copy() for r in matrix]
    rows,columns=list(range(len(m))),list(range(len(m[0])))
    pivots=[]
    for k in range(min(len(m),len(m[0]))):
        best,i,j=max((abs_lower(m[i][j]),i,j)
                     for i in range(k,len(rows)) for j in range(k,len(columns)))
        if not best:
            break
        rows[k],rows[i]=rows[i],rows[k]
        m[k],m[i]=m[i],m[k]
        columns[k],columns[j]=columns[j],columns[k]
        for row in m:
            row[k],row[j]=row[j],row[k]
        pivots.append(best)
        for i in range(k+1,len(rows)):
            a=m[i][k]/m[k][k]
            for j in range(k+1,len(columns)):
                m[i][j]-=a*m[k][j]
    r=len(pivots)
    return r,rows[:r],columns[:r],F(min(pivots),SCALE)


def manifest_check():
    raw=(HERE/'source_manifest.json').read_bytes()
    manifest=json.loads(raw)
    need(set(manifest['files'])=={'.gitattributes','README.md','check_root_box.py',
                                'root_box_jacobian.json','verify_all.py'},'complete manifest own inputs')
    need(set(manifest['dependencies'])=={
        'validated-maxima/verify_all.py','validated-maxima/dyadic.py','validated-maxima/certificates.json',
        'candidate-global-screen/verify_symmetry.py','candidate-global-screen/symmetry_certificates.json',
        'twelve-point-sharp-contacts/verify_all.py'},'complete manifest proof dependencies')
    for base,key in ((HERE,'files'),(ROOT,'dependencies')):
        for name,sha in manifest[key].items():
            need(hashlib.sha256((base/name).read_bytes()).hexdigest()==sha,'fixed input '+name)
    return hashlib.sha256(raw).hexdigest()


def decimal_bounds(interval,places=15):
    d=10**places
    lo,hi=interval.lo*d//SCALE,-(-interval.hi*d//SCALE)
    def fmt(n):
        sign='-' if n<0 else ''
        whole,part=divmod(abs(n),d)
        return f'{sign}{whole}.{part:0{places}d}'
    return [fmt(lo),fmt(hi)]


def run():
    manifest_sha=manifest_check()
    nodes,xy,means,ratios,contacts=candidate()
    reduced_box,norm,population_gap=check_reduced_box()
    need(all(l<F(v.lo,SCALE)<=F(v.hi,SCALE)<r
             for (l,r),v in zip(reduced_box,nodes+xy)),'reduced and full roots coincide')
    a,b,c,x,y=nodes+xy
    need((b*y+c*x).hi<0,'the four cross-rectangle angular contacts are distinct')
    R=(c.square()-b.square()).square()+4*(x-y)*(x*c.square()-y*b.square())
    bterm=a.square()-2*a*x+b.square()
    cterm=a.square()-2*a*y+c.square()
    disc=16*(b.square()-x.square())*(c.square()-y.square())*bterm.square()*cterm.square()*R.square()
    delta=4**5*a*b.square()*c.square()*disc.square()
    direct=ONE
    for contact in contacts:
        direct*=power(contact['d'],contact['weight'])
    need(delta.lo<=direct.hi and direct.lo<=delta.hi,'objective product cross-check')
    need(delta.lo>exact(F(97,80)*9**9).hi,'strictly above the old feasible threshold')
    rank_report={}
    for ward in (False,True):
        names,rows,matrix=build(nodes,means,ratios,contacts,7,ward,4)
        rank,rr,cc,pivot=lower_rank(matrix)
        need(rank==44 and len(names)==(127 if ward else 105),'minimal harmonic interpolation rank')
        bound=F(735,10**8) if ward else F(827,10**8)
        need(pivot>bound,'strict pivot lower bound')
        rank_report['gram_ward' if ward else 'gram']={
            'matrix_shape':[len(rows),len(names)],'rank':rank,'analytic_row_relations':6,
            'minimum_pivot_strict_rational_lower_bound':str(bound),
            'minor_rows':[rows[i] for i in rr],'minor_columns':[names[i] for i in cc]}
    return {
        'passed':True,
        'scope':'certified ten-point candidate, stationary equality rigidity, and necessary contact geometry; no majorant or global maximum claimed',
        'source_manifest_sha256':manifest_sha,
        'normalization':'mass 9, except where stated',
        'reduced_root':{
            'parameter_names':list('abcxy'),'parameter_enclosures':[decimal_bounds(v) for v in nodes+xy],
            'dyadic_box_center_numerators':[118037,1084060,1216217,-797427,810046],
            'dyadic_box_center_and_radius_exponent':20,
            'contraction_norm_strict_upper_bound':'21/10000',
            'original_force_root_inside_box':True,
            'unique_nonsingular':True,
            'candidate_objective_mass9':decimal_bounds(delta,12),
            'candidate_objective_mass10':decimal_bounds(exact(F(10,9)**45)*delta,12)},
        'stationary_rigidity':{
            'population_counts':[2,4,4],'population_triples_checked':66,
            'other_population_mass_distance_lower_bound':str(population_gap),
            'penalty_features':3},
        'contacts':{
            'count':len(contacts),'gram_boundary_count':sum(c['boundary'] for c in contacts),
            'same_radius_count':sum(c['same_layer'] for c in contacts),
            'all_complement_mass_strictly_greater_than_6':True,
            'all_cross_rectangle_t_distinct':True,
            'multiplicities':[c['weight'] for c in contacts],
            'minimal_possible_polynomial_harmonic_degree':7,
            'minimal_possible_polynomial_degree_with_linear_over_d_row_term':7,
            'contact_rank':rank_report},
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    output=(json.dumps(run(),indent=2)+'\n').encode('utf-8')
    if args.record_report:
        (HERE/'verification_report.json').write_bytes(output)
    sys.stdout.buffer.write(output)


if __name__=='__main__':
    main()
