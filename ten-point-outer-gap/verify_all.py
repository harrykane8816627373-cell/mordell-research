"""Exclude wide successive outer gaps and strengthen the angular Hessian bound."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import importlib.util
import json
from math import acos, factorial, floor, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
spec = importlib.util.spec_from_file_location('outer_gap_angular', ROOT/'ten-point-angular-upper/verify_all.py')
angular = importlib.util.module_from_spec(spec)
spec.loader.exec_module(angular)
curved, upper, need = angular.curved, angular.upper, angular.need

WEIGHT = F(21,1000)
SHORT = F(381,25000)
LONG = F(553,20000)
PAIR_FLOOR = F(11,5000)
COSINE = (F(59,100), F(39,50))
LO, HI = angular.COSINES[0]
CASES = {
    'short': (F(1,2), HI, SHORT, WEIGHT),
    'long': (LO, F(1,2), LONG, WEIGHT),
    'too-wide': (F(1,2), COSINE[0], PAIR_FLOOR, F(0)),
    'too-narrow': (COSINE[1], HI, PAIR_FLOOR, F(0)),
}


@lru_cache(maxsize=None)
def angle_lower(index):
    """Proposal by inverse cosine, accepted by a rational lower bound on cosine."""
    q = F(index,512)
    need(-1 < q < 1, 'arccos sample in physical range')
    candidate = F(floor(acos(float(q))*10**7)-1,10**7)
    need(0 < candidate < 2, 'cosine decreasing on the range used')
    # Last term is negative. The alternating tail decreases for 0 < candidate < 2.
    lower = sum((-1)**k*candidate**(2*k)/factorial(2*k) for k in range(12))
    need(lower > q, 'certified lower endpoint for arccos')
    return candidate


def outside(a, b, box):
    (ul,uh),(vl,vh),(tl,th)=box
    if angular.outside(1,1,1,box):
        return True
    if a>0 and (th<=0 or th**2 < a**2*ul*vl):
        return True
    if a<0 and th<0 and th**2 > a**2*uh*vh:
        return True
    return b>0 and tl>0 and tl**2 > b**2*uh*vh


def threshold(box, high, target, weight):
    if not weight:
        return -target
    (ul,uh),(vl,vh),(_,th)=box
    square = ul*vl if th>=0 else uh*vh
    units = 10**9
    root = isqrt((square.numerator*units**2)//square.denominator)
    root_bound = F(root if th>=0 else root+1, units)
    q = 512*min(high,th/root_bound)
    index = -(-q.numerator//q.denominator)
    # A slightly outlying negative box may still overlap a feasible point with cos >= LO.
    lo = LO*512
    index = max(index, -(-lo.numerator//lo.denominator))
    return weight*angle_lower(index)-target


def root_boxes(lo,hi):
    r = angular.BANDS[1]
    left = lo*r[int(lo<0)]
    right = hi*r[1]
    t = ((left,F(0)),(F(0),right)) if left<0 else ((left,right),)
    return [(r,r,ends) for ends in t]


def cover(name, poly, shape, stored):
    lo,hi,target,weight = CASES[name]
    packer = curved.Packed(shape)
    boxes = root_boxes(lo,hi)
    need(stored is None or len(stored)==len(boxes), 'all signed root boxes')
    trees, maximum = [], 0
    for idx,box in enumerate(boxes):
        original=box
        stack=[(box,upper.initial(poly,shape,box),0)]
        proof=None if stored is None else stored[idx]
        nodes=[]
        for item in proof if proof is not None else iter(lambda:bool(stack),False):
            need(bool(stack), 'node must correspond to a pending box')
            box,packed,depth=stack.pop()
            maximum=max(depth,maximum)
            if proof is None:
                if outside(lo,hi,box):
                    tag='P'
                else:
                    c=threshold(box,hi,target,weight)
                    if packer.dominates(packed,box,c):
                        tag='A'
                    elif packer.curved_dominates(packed,box,c):
                        tag='Q'
                    else:
                        need(depth<56,f'failed to cover {name} at box {box}')
                        widths=(F(1),F(1),F(2))
                        tag=str(max(range(3),key=lambda a:(box[a][1]-box[a][0])/widths[a]))
            else:
                tag=item
            if tag in '012':
                stack.extend(curved.split(packer,packed,box,depth,int(tag)))
            elif proof is not None:
                if tag=='P':
                    accepted=outside(lo,hi,box)
                else:
                    c=threshold(box,hi,target,weight)
                    accepted=(tag=='A' and packer.dominates(packed,box,c)
                              or tag=='Q' and packer.curved_dominates(packed,box,c))
                need(accepted,f'accepted {tag} leaf in {name}')
            nodes.append(tag)
        need(not stack and len(nodes)==1+2*sum(n in '012' for n in nodes), f'complete cover of {original}')
        trees.append(''.join(nodes))
        if proof is None:
            print(name,idx,len(nodes),flush=True)
    return trees, {'nodes':sum(map(len,trees)),'depth':maximum}


def setup():
    data=angular.read(curved.HERE/'kernel.json')
    A,force=upper.corrections(data)
    need(bool(force), 'denominator-cleared curved kernel')
    base,shape=upper.polynomial(data,A,force)
    params=curved.parsed(data)
    poly=curved.root_polynomial(base,shape,params,1,1)
    curved.check_root(base,poly,params,1,1)
    return data,poly,shape


def combine_budgets(data):
    report=angular.read(angular.HERE/'verification_report.json')
    need(report['passed'] and report['unresolved']==0 and report['normalization_mass']==9,
         'accepted angular floors')
    rows=report['classes']
    need(tuple(tuple(x['class']) for x in rows)==angular.CLASSES and
         tuple(x['weight'] for x in rows)==angular.WEIGHTS, 'complete angular floor counts')
    floors=[-F(row['deficit_floor']) for row in rows]
    all_edges=sum(c*n for c,n in zip(floors,angular.WEIGHTS))
    gap=upper.hx(data['B'])-upper.log.log_lower(upper.actual_root_lower())
    need(LONG>SHORT and 7*SHORT+LONG-WEIGHT*F(44,7)>gap+all_edges-8*floors[2],
         'one long edge exceeds the remaining budget after summing all eight successive gaps')
    need(PAIR_FLOOR>gap+all_edges-floors[2], 'one forbidden edge exceeds the other 44 edge budget')


def curvature():
    ol,oh=angular.BANDS[1]
    q=F(41,40)
    need((ol+oh)**2 < 4*q*q*ol*oh, 'q=(u+v)/(2 sqrt(uv)) bounded on the full outer band')
    lo,hi=COSINE
    need(2-q*q-q*hi>0 and hi*q+hi**2<2, 'curvature increases in cosine and decreases in q')
    base=1/(1+q)
    one_step,two_steps=F(52,25),F(37,50)
    doubled=2*lo*lo-1
    need(lo>F(1,2) and (1-q*lo)/(q-lo)**2>one_step
         and (1-q*doubled)/(q-doubled)**2>two_steps,
         'strict lower bounds for consecutive and two-step edges')
    old=angular.read(ROOT/'ten-point-outer-angles/verification_report.json')
    need(old['passed'] and old['normalization_mass']==9 and old['outer_count']==8
         and F(old['refined_bands']['fixed_mean_outer_angular_negative_hessian_lower'])==F(173,100),
         'accepted angular curvature inputs')
    need(F(17,12)**2>2, 'frequency-one cycle gap is greater than 7/12')
    additional=min((one_step-base)*F(7,12)+2*(two_steps-base), 2*(one_step-base))
    outer=8*base+additional
    near=F(113,250)
    far=F(23,100)
    need(angular.BANDS[0][1]/ol < near**2 and far**2 < angular.BANDS[0][0]/oh
         and 2*near/(1-near)**2 < F(151,50), 'every near inner-outer curvature is below 151/50')
    need((1-F(77,100))/2 < F(17,50)**2 and 2*far/(1+far**2) > F(17,50),
         'for each outer point the farther inner direction has nonpositive curvature')
    need(outer-F(219,100)>F(317,100) and outer-F(151,50)>F(117,50),
         'both strengthened mean-zero angular Hessian bounds')


def verify(build, selected):
    if not build:
        manifest=angular.read(HERE/'source_manifest.json')
        for base,kind in ((HERE,'files'),(ROOT,'dependencies')):
            for name,h in manifest[kind].items():
                need(sha256((base/name).read_bytes()).hexdigest()==h,'fixed input '+name)
        angular.check_domains()
        curved.curvature_selftest()
        upper.old.selftest()
    data,poly,shape=setup()
    result={}
    for name in CASES if selected is None else (selected,):
        path=HERE/(name+'.json')
        frozen=None if build else angular.read(path)
        trees,result[name]=cover(name,poly,shape,frozen)
        if build:
            path.write_text(json.dumps(trees)+'\n',encoding='ascii',newline='\n')
    combine_budgets(data)
    curvature()
    if selected is not None:
        return {'passed':True,'covered_case':selected,**result[selected]}
    return {'passed':True,'applies_to_every_global_maximizer':True,'normalization_mass':9,
            'outer_successive_angle_cosine_open_interval':list(map(str,COSINE)),
            'outer_mean_zero_negative_angle_hessian_lower_with_old_min_distance':'317/100',
            'outer_mean_zero_negative_angle_hessian_lower_on_successive_gap_polytope':'117/50',
            'half_turn_from_inner_antipodes_and_matched_radii_without_branch_condition':True,
            'nodes':sum(r['nodes'] for r in result.values()),'cases':result,'unresolved':0}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',choices=CASES)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    result=verify(args.build is not None,args.build)
    path=HERE/'verification_report.json'
    output=json.dumps(result,indent=2)+'\n'
    if args.record_report:
        need(args.build is None,'replay all cases to record a verified report')
        path.write_text(output,encoding='ascii',newline='\n')
    elif args.build is None:
        need(result==angular.read(path),'saved proof report')
    print(output,end='')


if __name__=='__main__':
    main()
