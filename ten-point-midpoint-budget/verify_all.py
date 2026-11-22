"""Certify a joint midpoint bound from the four opposite-rank outer edges and the inner edge."""
import argparse
import importlib.util
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from outer_endpoint import prove_endpoint, prove_monotonicity

HERE = Path(__file__).resolve().parent
LAMBDA = F(69,500)
MU = -F(14,25000)
CENTERS = F(6,25), F(219,100)
COMPENSATIONS = F(93,200000), F(491,200000)
MIDPOINT_NORM_CAP = F(149,500)
BANDS = ((F(37,500),F(43,250)),(F(169,200),F(164,125)))
D = {(1,0,0): F(1), (0,1,0): F(1), (0,0,1): F(-2)}


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def require(condition,message):
    if not condition:
        raise ArithmeticError(message)


def load(path,name):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_inputs(repo):
    manifest=read(HERE/'source_manifest.json')
    for key,base in (('files',HERE),('dependencies',repo)):
        for relative,expected in manifest[key].items():
            require(digest(base/relative)==expected,'fixed input '+relative)
    radii=read(repo/'ten-point-refined-radii/verification_report.json')
    require(radii['passed'] and radii['applies_to_every_global_maximizer'] and radii['unresolved']==0,
            'accepted global radius and inner midpoint constraints')
    require((radii['exact_inner_count'],radii['exact_outer_count'])==(2,8) and
            tuple(tuple(map(F,radii[k])) for k in ('inner_squared_radius_open_interval','outer_squared_radius_open_interval'))==BANDS
            and F(radii['midpoint']['inner_midpoint_squared_strict_upper'])==F(57121,4000000),
            'fixed inner and outer bands and inner midpoint bound')
    condition=read(repo/'ten-point-conditioned-inner/verification_report.json')
    require(condition['passed'] and condition['applies_to_every_global_maximizer'] and condition['unresolved']==0
            and tuple(map(F,condition['inner_squared_distance_open_interval']))==(F(41,100),F(11,20))
            and F(condition['inner_angle_cosine_strict_upper'])==-F(77,100),'accepted inner-pair distance and angle')
    paired=read(repo/'ten-point-tight-pairing/verification_report.json')
    require(paired['passed'] and paired['applies_to_every_global_maximizer'] and paired['normalization_mass']==9
            and F(paired['pairing']['cyclic_rank_pairing_cosine_strict_upper'])==-F(31,32)
            and F(paired['pairing']['paired_squared_distance_strict_lower'])==F(10647,3200)
            and paired['pairing']['closed_arc_occupancy_upper']==4,'accepted four opposite-rank outer edges')
    return read(repo/'ten-point-curved-kernel/verification_report.json')


def nonmatching_compensation(repo):
    report=read(repo/'ten-point-angular-upper/verification_report.json')
    require(report['passed'] is True and report['normalization_mass']==9 and report['unresolved']==0,
            'accepted exact angular cover')
    expected_classes=[([0,1,0],16,-F(97,2000000)),([1,1,1],8,-F(103,2000000)),
                      ([1,1,2],8,-F(477,10000000)),([1,1,3],8,-F(479,10000000))]
    selected=report['classes'][1:5]
    require([(row['class'],row['weight'],F(row['deficit_floor'])) for row in selected]==expected_classes,
            'verified overlapping neither the inner pair nor the four outer matching edges')
    return -sum(weight*floor for kind,weight,floor in expected_classes)


def necessary_inner(curved,box):
    (ul,uh),(vl,vh),(tl,th)=box
    if curved.prune(box,True):
        return False
    return not (uh+vh-2*tl<F(41,100) or ul+vl-2*th>F(11,20)
                or tl**2<F(77,100)**2*ul*vl or ul+vl+2*tl>F(239,1000)**2)


def shifted(curved,polynomial,band):
    affine={(0,0,0):-MU*CENTERS[band],(1,0,0):LAMBDA/4+MU,
            (0,1,0):LAMBDA/4+MU,(0,0,1):LAMBDA/2}
    correction=curved.multiply(D,affine)
    result={key:polynomial.get(key,F(0))-correction.get(key,F(0))
            for key in polynomial.keys() | correction.keys()}
    u,v,t=(F(1,8),F(1,7),F(-1,10)) if band==0 else (F(1),F(11,10),F(-1))
    du=u+v-2*t
    expected=du*(LAMBDA*(u+v+2*t)/4+MU*(u+v-CENTERS[band]))
    require(expected==sum(coeff*u**i*v**j*t**k for (i,j,k),coeff in correction.items()),
            'exact comparison with midpoint and radius correction')
    return result


def cover_inner(curved,pack,poly,shape,stored):
    root=(BANDS[0],BANDS[0],(F(-27,200),-F(1,10)))
    widths=[b-a for a,b in root]
    stack=[(root,curved.upper.initial(poly,shape,root),0)]
    nodes=[]
    max_depth=0
    c=COMPENSATIONS[0]
    while stack:
        box,packed,depth=stack.pop()
        max_depth=max(depth,max_depth)
        if stored is not None:
            require(len(nodes)<len(stored),'all children of each split are present')
            tag=stored[len(nodes)]
        elif not necessary_inner(curved,box):
            tag='P'
        elif pack.dominates(packed,box,c):
            tag='A'
        elif pack.curved_dominates(packed,box,c):
            tag='Q'
        else:
            require(depth<66,f'locate failed box at depth 66: {box}')
            scale=[(b-a)/width for (a,b),width in zip(box,widths)]
            tag=str(max(range(3),key=lambda axis: scale[axis]))
        if tag in ('0','1','2'):
            stack.extend(curved.split(pack,packed,box,depth,int(tag)))
        else:
            require((tag=='P' and not necessary_inner(curved,box)) or
                    (tag=='A' and pack.dominates(packed,box,c)) or
                    (tag=='Q' and pack.curved_dominates(packed,box,c)),
                    f'valid leaf {tag} in the inner band')
        nodes.append(tag)
    require(stored is None or len(nodes)==len(stored),'cover tree consumed exactly')
    counts=Counter(nodes)
    require(sum(counts[tag] for tag in '012')+1==sum(counts[tag] for tag in 'PAQ'),
            'binary tree conservation')
    return ''.join(nodes), {'visited':len(nodes),'max_depth':max_depth,'counts':dict(sorted(counts.items()))}


def run(repo,build):
    curved_path=repo/'ten-point-curved-kernel'
    curved=load(curved_path/'verify_all.py','midpoint_budget_curved')
    upper=curved.upper
    report=read(curved_path/'verification_report.json') if build else check_inputs(repo)
    kernel_path=curved_path/'kernel.json'
    data=read(kernel_path)
    require((data['n'],data['mass'])==(10,9) and report['passed'] and report['unresolved']==0
            and report['kernel_sha256']==digest(kernel_path) and report['PSD_blocks']==10,
            'accepted curved kernel of mass nine and positive definite Gram blocks')
    require(CENTERS[0]+4*CENTERS[1]==9,'the affine radius contribution sums to zero')
    coefficients,force=upper.corrections(data)
    old,shape=upper.polynomial(data,coefficients,force)
    require(bool(force),'the curved kernel clears denominator d')
    params=curved.parsed(data)
    pack=curved.Packed(shape)
    datafile=HERE/'cover.json'
    trees=[None,None,None] if build else read(datafile)
    require(len(trees)==3,'inner cover and outer monotonicity and endpoint covers')
    inner=shifted(curved,curved.root_polynomial(old,shape,params,0,0),0)
    outer=shifted(curved,curved.root_polynomial(old,shape,params,1,1),1)
    trees[0],inner_proof=cover_inner(curved,pack,inner,shape,trees[0])
    if build:
        print('inner',len(trees[0]),flush=True)
    trees[1]=prove_monotonicity(curved,outer,BANDS[1],trees[1])
    if build:
        print('outer monotonicity',len(trees[1]),flush=True)
    trees[2]=prove_endpoint(curved,outer,BANDS[1],COMPENSATIONS[1],trees[2])
    proof={'inner':inner_proof,
           'outer_monotonicity':{'visited':len(trees[1]),'counts':dict(sorted(Counter(trees[1]).items()))},
           'outer_endpoint':{'visited':len(trees[2]),'counts':dict(sorted(Counter(trees[2]).items()))}}
    if build:
        datafile.write_text(json.dumps(trees,separators=(',',':'))+'\n',encoding='ascii',newline='\n')
    b=upper.hx(data['B'])
    witness=upper.witness()
    candidate_lower=upper.actual_root_lower()
    require(candidate_lower>witness,'certified candidate lower bound exceeds the integer witness')
    nonmatching=nonmatching_compensation(repo)
    budget=b-upper.log.log_lower(candidate_lower)+nonmatching+COMPENSATIONS[0]+4*COMPENSATIONS[1]
    require(budget<LAMBDA*MIDPOINT_NORM_CAP**2,'strict midpoint bound by summing all 45 edges')
    result={'passed':True,'applies_to_every_global_maximizer':True,'normalization_mass':9,
            'strict_joint_midpoint_norm_upper':str(MIDPOINT_NORM_CAP),
            'lambda':str(LAMBDA),'mu':str(MU),'pair_gap_compensations':list(map(str,COMPENSATIONS)),
            'nonmatching_compensation':str(nonmatching),
            'log_budget_margin_decimal_display':float(LAMBDA*MIDPOINT_NORM_CAP**2-budget),
            'cover':proof,'unresolved':0}
    saved=HERE/'verification_report.json'
    content=json.dumps(result,indent=2)+'\n'
    if build:
        saved.write_text(content,encoding='ascii',newline='\n')
    else:
        require(result==read(saved),'saved proof report')
    print(content,end='')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=HERE.parent)
    parser.add_argument('--build-cover',action='store_true')
    args=parser.parse_args()
    run(args.repo.resolve(),args.build_cover)
