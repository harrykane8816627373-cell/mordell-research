"""Certify a ten-point contact fit with cyclic row and second-variation terms."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CONTACT_SOURCE = ROOT/'ten-point-contact-geometry'
sys.path.insert(0, str(CONTACT_SOURCE))
spec = importlib.util.spec_from_file_location('ten_point_contact_geometry',CONTACT_SOURCE/'verify_all.py')
geometry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geometry)

sharp = geometry.ref
Ball, Jet, exact, ZERO, ONE = geometry.Ball, geometry.Jet, geometry.exact, geometry.ZERO, geometry.ONE
need, dot = sharp.need, sharp.dot
BANDS = ((F(21,250), F(159,1000)), (F(89,100), F(1271,1000)))
HARMONICS = 11
RIGIDITY_WEIGHT = F(1,2**40)


def read(path):
    return json.loads((ROOT/path).read_bytes())


def angle_classes():
    bands = read('ten-point-angular-radii/verification_report.json')
    gap = read('ten-point-outer-gap/verification_report.json')
    pairing = read('ten-point-tight-pairing/verification_report.json')
    for row in (bands,gap,pairing):
        need(row['passed'] and row['normalization_mass']==9 and
             row['applies_to_every_global_maximizer'], 'global geometric hypotheses')
    need((bands['exact_inner_count'],bands['exact_outer_count'])==(2,8), 'two inner and eight outer')
    for label,band in zip(('inner','outer'),BANDS):
        need(tuple(map(F,bands[f'{label}_squared_radius_open_interval']))==band, 'radial bands')
    lo,hi = map(F,gap['outer_successive_angle_cosine_open_interval'])
    near_opposite = F(pairing['pairing']['cyclic_rank_pairing_cosine_strict_upper'])
    need(F(1,2)<lo<hi<1, 'three successive angles are less than pi')
    cosine3 = lambda c: 4*c**3-3*c
    spans = [(lo,hi),(2*lo**2-1,2*hi**2-1),(cosine3(lo),cosine3(hi)),(F(-1),near_opposite)]
    need(all(spans[i+1][1]<spans[i][0] for i in range(3)), 'disjoint cyclic classes')
    return dict(enumerate(spans,1))


def classify(contacts, spans):
    counts = [0]*6
    for c in contacts:
        u,v,t = (c[s] for s in 'uvt')
        for label,z in zip(c['layers'],(u,v)):
            lower,upper = BANDS[int(label>0)]
            need(exact(lower).hi<z.lo<=z.hi<exact(upper).lo, 'candidate is in the radial band')
        if c['layers'][1]==0:
            rank=0
        elif c['layers'][0]==0:
            rank=1
        elif c['boundary']:
            rank=5
        else:
            cosine=t/sharp.sym.sqrt_ball(u*v,Ball,geometry.SCALE)
            ranks=[i+1 for i,(a,b) in spans.items() if exact(a).hi<cosine.lo<=cosine.hi<exact(b).lo]
            need(len(ranks)==1, 'candidate contact lies strictly in an angular class')
            rank=ranks[0]
        counts[rank]+=c['weight']
        c['class']=rank
    need(counts==[1,16,8,8,8,4], 'exact contact populations by class')
    return counts


def radial_features(u, band, unit):
    lo,hi = BANDS[band]
    s = (2*u-unit*exact(lo+hi))/exact(hi-lo)
    raw = [unit,s,2*s*s-unit]
    return raw+[0*unit]*3 if band==0 else [0*unit]*3+raw


def real_powers(radius, component):
    values = [ONE,component]
    for _ in range(2,HARMONICS//2+1):
        values.append(2*component*values[-1]-radius.square()*values[-2])
    return values


def projections(nodes, xy, pivots):
    rows = [radial_features(u,int(i>0),ONE) for i,u in enumerate(nodes)]
    powers = [real_powers(u,v) for u,v in zip(nodes,[nodes[0],*xy])]
    result = []
    for k,pivot in enumerate(pivots):
        if k%2:
            need(pivot is None, 'odd moments vanish before projection')
            result.append(None)
            continue
        need(isinstance(pivot,int) and 0<=pivot<6, 'valid even projection pivot')
        moments = [sum((rows[i][j]*powers[i][k//2]*n for i,n in enumerate((2,4,4))),ZERO)
                   for j in range(6)]
        need(sharp.abs_lower(moments[pivot])>0, 'candidate moment pivot excludes zero')
        result.append((pivot,[v/moments[pivot] for v in moments]))
    return result


def project(row, projection):
    if projection is None:
        return row
    pivot,ratios = projection
    return [v-row[pivot]*ratios[i] for i,v in enumerate(row) if i!=pivot]


def reciprocal(jet):
    inverse = ONE/jet.v
    first,second = -inverse.square(),2*inverse.square()*inverse
    return Jet(inverse,tuple(x*first for x in jet.d),
               tuple(tuple(jet.h[i][j]*first+jet.d[i]*jet.d[j]*second for j in range(3))
                     for i in range(3)))


def squared_features(row):
    return [(1+int(i!=j))*row[i]*row[j] for i,j in combinations_with_replacement(range(len(row)),2)]


def feature_columns(u,v,t,unit,bands,cls,nodes,transforms,cleared):
    d = u+v-2*t
    bi,bj = bands
    ru,rv = radial_features(u,bi,unit),radial_features(v,bj,unit)
    fu,fv = ([project(row,p) for p in transforms] for row in (ru,rv))
    harmonics = [unit,t]
    for _ in range(2,HARMONICS+1):
        harmonics.append(2*t*harmonics[-1]-u*v*harmonics[-2])
    values = []
    for k,(a,b) in enumerate(zip(fu,fv)):
        for i,j in combinations_with_replacement(range(len(a)),2):
            mixed = a[i]*b[j]+(a[j]*b[i] if i!=j else 0*unit)
            diagonal = (u**k*a[i]*a[j]+v**k*b[i]*b[j])*(1+int(i!=j))
            values.append(-harmonics[k]*mixed-diagonal/18)
    for k,c0 in ((0,F(9,10)),(1,0),(1,1),(2,F(17,10)),(3,F(8,5))):
        wu,wv = (u-c0,v-c0) if c0 else (unit,unit)
        for x,y in zip(fu[k],fv[k]):
            values.append(-harmonics[k]*(wu*y+wv*x)-(u**k*wu*x+v**k*wv*y)/9)
    values.extend(F(8-bi,9)*x-bj*x+F(8-bj,9)*y-bi*y for x,y in zip(ru,rv))
    values.extend((v-1+u/9)*x+(u-1+v/9)*y for x,y in zip(fu[0],fv[0]))
    for k in range(1,4):
        factor = int(cls==k)-2*int(cls==4)
        values.extend(factor*(x+y) for x,y in zip(ru[3:],rv[3:]))

    a,b,c_out = nodes
    def zero_at_nodes(s,band):
        return (s-a)/F(49,1000) if band==0 else (s-b)*(s-c_out)/F(467,2000)**2
    gu = [zero_at_nodes(u,bi)*x for x in ru]
    gv = [zero_at_nodes(v,bj)*x for x in rv]
    if cleared:
        values = [value*d**2 for value in values]
        values.extend(d*((u-t-5*u*d/9)*x+(v-t-5*v*d/9)*y) for x,y in zip(gu,gv))
    else:
        inv = reciprocal(d)
        values.extend(((u-t)*inv-5*u/9)*x+((v-t)*inv-5*v/9)*y for x,y in zip(gu,gv))
    norm = [(u*x+v*y)*F(5,9) for x,y in zip(squared_features(gu),squared_features(gv))]
    radial = squared_features([(u-t)*x+(v-t)*y for x,y in zip(gu,gv)])
    angular = squared_features([x-y for x,y in zip(gu,gv)])
    second = [x-(u*v-t*t)*y for x,y in zip(radial,angular)]
    if cleared:
        norm = [value*d**2 for value in norm]
    else:
        second = [value*inv**2 for value in second]
    values.extend(-n-s for n,s in zip(norm,second))
    values.extend(-n+s for n,s in zip(norm,second))
    return values


def columns(c,nodes,transforms):
    coords = [Jet(c[key],tuple(ONE if i==j else ZERO for i in range(3))) for j,key in enumerate('uvt')]
    return feature_columns(*coords,Jet(1),tuple(int(label>0) for label in c['layers']),
                           c['class']-1,nodes,transforms,False)


def rigidity_term(u,v,nodes):
    def feature(s):
        a,b,c = nodes
        return ((s-a)*(s-b)*(s-c))**2
    fu,fv = feature(u),feature(v)
    return -(fu*fv+(fu**2+fv**2)/18)*RIGIDITY_WEIGHT


def build(contacts,nodes,transforms):
    matrix,target = [],[]
    logs = [sharp.log_ball(c['d']) for c in contacts]
    budget = sum((c['weight']*value for c,value in zip(contacts,logs)),ZERO)
    for c,log_d in zip(contacts,logs):
        cols = columns(c,nodes,transforms)
        matrix.append([j.v for j in cols])
        target.append(log_d-budget/45+1-(c['u']+c['v'])*exact(F(5,9)))
        inv = ONE/c['d']
        rad = inv-exact(F(5,9))
        if c['boundary']:
            matrix.append([j.d[0]-j.d[2]/2 for j in cols])
            target.append(rad+inv)
        else:
            for i in ((0,2) if c['same_layer'] else (0,1,2)):
                matrix.append([j.d[i] for j in cols])
                target.append(rad if i<2 else -2*inv)
    need(len(matrix)==len(target)==50, 'full value and first-derivative system')
    return matrix,target,budget


def correct(matrix,target,seed):
    parameters = [exact(F(v)) for v in seed['parameters']]
    rows,cols = seed['minor_rows'],seed['minor_columns']
    need(len(set(rows))==len(set(cols))==44 and all(0<=i<50 for i in rows) and
         all(0<=i<len(parameters) for i in cols), '44 by 44 minor indices')
    need(all(len(row)==len(parameters) for row in matrix), 'complete coefficients')
    residual = [v-dot(row,parameters) for row,v in zip(matrix,target)]
    pivots,determinant,adjust = sharp.certify_minor(matrix,rows,cols,residual)
    need(sharp.abs_lower(determinant)>0, 'sharp contact determinant excludes zero')
    for i,d in zip(cols,adjust):
        parameters[i]+=d
    differences = [dot(row,parameters)-v for row,v in zip(matrix,target)]
    need(all(e.lo<=0<=e.hi for e in differences), 'interval check of the dependent rows')
    return parameters,min(pivots)


def check_grams(parameters):
    size = len(parameters)
    offset,base,second = 0,[],[]
    for k,dim in enumerate([5,6]*6+[6,6]):
        if k==12:
            offset=size-42
        matrix=[[ZERO]*dim for _ in range(dim)]
        for i,j in combinations_with_replacement(range(dim),2):
            matrix[i][j]=matrix[j][i]=parameters[offset]
            offset+=1
        pivot=sharp.ldl_lower(matrix,f'Gram block {k}')
        (base if k<=HARMONICS else second).append(pivot)
    need(offset==size, 'all second-variation coefficient blocks inspected')
    return base,second


def gap_hessian(c,nodes,transforms,weights):
    u,v,t = (c[x] for x in 'uvt')
    inv = ONE/(u+v-2*t)
    jets = columns(c,nodes,transforms)
    axis = [1,1,-2]
    hessian = [[dot(weights,[jet.h[i][j] for jet in jets])+axis[i]*axis[j]*inv.square()
                for j in range(3)] for i in range(3)]
    inward = dot(weights,[jet.d[2] for jet in jets])+2*inv
    extra = rigidity_term(Jet(u,(ONE,ZERO,ZERO)),Jet(v,(ZERO,ONE,ZERO)),nodes)
    hessian = [[hessian[i][j]+extra.h[i][j] for j in range(3)] for i in range(3)]
    return hessian,inward


def tangent_hessian(hessian,inward,u,v,t):
    tangents = [[ONE,ZERO,t/(2*u)],[ZERO,ONE,t/(2*v)]]
    mixed = t/(4*u*v)
    second = [[-t/(4*u.square()),mixed],[mixed,-t/(4*v.square())]]
    return [[sum((tangents[i][a]*hessian[a][b]*tangents[j][b]
                  for a in range(3) for b in range(3)),ZERO)+inward*second[i][j]
             for j in range(2)] for i in range(2)]


def local_boxes(contacts,nodes,transforms,weights,radii):
    need(len(radii)==len(contacts), 'one dyadic radius per contact')
    result=[]
    for i,(c,radius) in enumerate(zip(contacts,radii)):
        need(0<radius<=F(1,256) and radius.denominator & (radius.denominator-1)==0,'dyadic local radius')
        u,v,t = (c[x] for x in 'uvt')
        n = exact(radius).hi
        ur,vr = u+Ball(-n,n),v+Ball(-n,n)
        width = 2*n if c['boundary'] else n
        tr = t+Ball(-width,width)
        need(min(ur.lo,vr.lo,(ur+vr-2*tr).lo)>0, 'positive local radii and distance')
        hessian,inward = gap_hessian(dict(c,u=ur,v=vr,t=tr),nodes,transforms,weights)
        if c['boundary']:
            face = -sharp.sym.sqrt_ball(ur*vr,Ball,geometry.SCALE)
            need(inward.lo>0 and tr.lo<=face.lo<=face.hi<=tr.hi, 'inward sign and Gram projection')
            hessian=tangent_hessian(hessian,inward,ur,vr,face)
        pivot=sharp.ldl_lower(hessian,f'local contact {i}')
        box = [[F(z.hi,geometry.SCALE)-radius,F(z.lo,geometry.SCALE)+radius] for z in (u,v,t)]
        need(all(lo<hi for lo,hi in box), 'nonempty covered box')
        result.append({'index':i, 'class':c['class'], 'radius':str(radius),
                       'cover_box':[list(map(str,pair)) for pair in box],
                       'hessian_pivot_lower_power':sharp.lower_power(pivot)-1})
    return result


def negative_pair(nodes,transforms,weights,budget,spans):
    u0,v0,t0 = F(119,100),F(119,100),F(3511,5000)
    lo,hi = spans[1]
    need(BANDS[1][0]<u0==v0<BANDS[1][1] and lo<t0/u0<hi and u0+v0-2*t0>F(2,5),
         'rational pair strictly inside the radius, angle and distance bands')
    u,v,t = map(exact,(u0,v0,t0))
    c = {'u':u,'v':v,'t':t,'layers':(1,1),'class':2}
    jets = columns(c,nodes,transforms)
    h = (budget/45-1+(u+v)*exact(F(5,9))+dot(weights,[j.v for j in jets])+
         rigidity_term(Jet(u),Jet(v),nodes).v)
    gap = h-sharp.log_ball(u+v-2*t)
    ceiling = F(-1,40000)
    need(gap.hi<exact(ceiling).lo,'strict negative gap for the current seed')
    return {'class':2,'point':list(map(str,(u0,v0,t0))),'strict_gap_upper':str(ceiling)}


def exact_kernel():
    manifest = read('ten-point-sharp-contacts/source_manifest.json')
    for path,digest in manifest.items():
        need(sha256((ROOT/path).read_bytes()).hexdigest()==digest,'fixed input '+path)
    spans=angle_classes()
    nodes,xy,_,_,contacts=geometry.candidate()
    counts=classify(contacts,spans)
    seed=read('ten-point-sharp-contacts/seed.json')
    need(len(seed['harmonic_pivots'])==HARMONICS+1,'one transform per harmonic')
    transforms=projections(nodes,xy,seed['harmonic_pivots'])
    matrix,target,budget=build(contacts,nodes,transforms)
    weights,minor=correct(matrix,target,seed)
    grams,variations=check_grams(weights)
    boxes=local_boxes(contacts,nodes,transforms,weights,list(map(F,seed['local_radii'])))
    obstruction=negative_pair(nodes,transforms,weights,budget,spans)
    report = {'passed':True,'mass':9,'angular_class_cosines':{str(k):list(map(str,r)) for k,r in spans.items()},
            'pair_counts_by_class':counts,'contact_matrix_shape':[len(matrix),len(weights)],
            'contact_rank':44,'contact_minor_pivot_lower_power':sharp.lower_power(minor)-1,
            'positive_gram_blocks':len(grams),'gram_pivot_lower_power':sharp.lower_power(min(grams))-1,
            'positive_second_variation_blocks':len(variations),
            'second_variation_pivot_lower_power':sharp.lower_power(min(variations))-1,
            'radial_rigidity_weight':str(RIGIDITY_WEIGHT),
            'local_boxes':boxes,'rational_negative_pair':obstruction,
            'continuous_complement_covered':False,'sharp_global_value_proved':False}
    return report,(weights,nodes,transforms,budget,spans)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    report,_=exact_kernel()
    path=HERE/'verification_report.json'
    encoded=json.dumps(report,indent=2)+'\n'
    if args.record_report:
        path.write_bytes(encoded.encode())
    else:
        need(json.loads(path.read_bytes())==report,'matches frozen exact replay')
    print(encoded,end='')


if __name__=='__main__':
    main()
