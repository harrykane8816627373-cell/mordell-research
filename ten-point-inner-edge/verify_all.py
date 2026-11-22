"""Certify the inner edge directly from the 44 complementary pair budgets."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MIDPOINT = F(87, 1000)
LIMIT = F(1073, 500000)
DISTANCE = F(43, 100), F(53, 100)
ANGLE = -F(87, 100)
PARALLEL = F(7, 125)
STIFFNESS = F(53, 2)


def need(test, detail):
    if not test:
        raise ArithmeticError(detail)


def read(path):
    return json.loads(path.read_bytes())


def load(path):
    spec = importlib.util.spec_from_file_location('inner_edge_angular',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def setup(repo):
    manifest = read(HERE/'source_manifest.json')
    for base,key in ((HERE,'files'),(repo,'dependencies')):
        for name,h in manifest[key].items():
            need(sha256((base/name).read_bytes()).hexdigest()==h,'fixed input '+name)
    angular = load(repo/'ten-point-angular-upper/verify_all.py')
    angular.check_domains()
    curved,upper = angular.curved,angular.upper
    data = read(curved.HERE/'kernel.json')
    rad = read(repo/'ten-point-angular-radii/verification_report.json')
    mid = read(repo/'ten-point-inner-midpoint/verification_report.json')
    ang = read(repo/'ten-point-angular-upper/verification_report.json')
    need(all(x['passed'] and x['normalization_mass']==9 and x['unresolved']==0
             for x in (rad,mid,ang)) and rad['applies_to_every_global_maximizer'] and
         mid['applies_to_every_global_maximizer'] and F(mid['inner_midpoint_norm_strict_upper'])==MIDPOINT,
         'new radial and midpoint domain')
    need((rad['exact_inner_count'],rad['exact_outer_count'])==(2,8), 'populations')
    inner,outer = (tuple(map(F,rad[k])) for k in ('inner_squared_radius_open_interval','outer_squared_radius_open_interval'))
    counts = [(row['class'],row['weight']) for row in ang['classes']]
    need(counts==[([0,0,0],1),([0,1,0],16),([1,1,1],8),([1,1,2],8),([1,1,3],8),([1,1,4],4)],
         'all 44 other edges present with multiplicity')
    loss = -sum(row['weight']*F(row['deficit_floor']) for row in ang['classes'][1:])
    need(upper.hx(data['B'])-upper.log.log_lower(upper.actual_root_lower())+loss<LIMIT,
         'strict upper bound on the deficit of the unique inner edge')
    A,force = upper.corrections(data)
    need(bool(force),'denominator-cleared curved kernel')
    base,shape = upper.polynomial(data,A,force)
    params = curved.parsed(data)
    poly = curved.root_polynomial(base,shape,params,0,0)
    curved.check_root(base,poly,params,0,0)
    return angular,inner,outer,shape,poly


def geometry(box):
    (ul,uh),(vl,vh),(tl,th) = box
    d0,d1 = ul+vl-2*th,uh+vh-2*tl
    imbalance = max(abs(uh-vl),abs(vh-ul))
    return (DISTANCE[0]<d0 and d1<DISTANCE[1] and th<0 and th**2>ANGLE**2*uh*vh
            and imbalance**2<4*d0*PARALLEL**2)


def cover(context, proof):
    angular,inner,_,shape,poly = context
    curved,upper = angular.curved,angular.upper
    box = (inner,inner,angular.root_boxes(0,0,0)[0][2])
    pack = curved.Packed(shape)
    widths = [b-a for a,b in box]
    todo = [(box,upper.initial(poly,shape,box),0)]
    tags = []
    max_depth = 0
    while todo:
        box,packed,depth = todo.pop()
        max_depth = max(max_depth,depth)
        if proof is not None:
            need(len(tags)<len(proof),'complete inner-edge tree')
            action = proof[len(tags)]
        elif (angular.outside(0,0,0,box) or box[0][0]+box[1][0]+2*box[2][0]>4*MIDPOINT**2):
            action = 'P'
        elif geometry(box):
            action = 'G'
        elif pack.dominates(packed,box,-LIMIT):
            action = 'A'
        elif pack.curved_dominates(packed,box,-LIMIT):
            action = 'Q'
        else:
            need(depth<54,'unresolved inner-edge box '+str([[float(x) for x in b] for b in box]))
            action = str(max(range(3),key=lambda i:(box[i][1]-box[i][0])/widths[i]))
        if action in ('0','1','2'):
            todo.extend(curved.split(pack,packed,box,depth,int(action)))
        elif proof is not None:
            valid = (action=='P' and (angular.outside(0,0,0,box) or box[0][0]+box[1][0]+2*box[2][0]>4*MIDPOINT**2)
                     or action=='G' and geometry(box) or action=='A' and pack.dominates(packed,box,-LIMIT)
                     or action=='Q' and pack.curved_dominates(packed,box,-LIMIT))
            need(valid,'valid recorded geometry or deficit leaf')
        tags.append(action)
    need(proof is None or len(proof)==len(tags),'no extra recorded leaves')
    return ''.join(tags),max_depth


def stiffness(outer):
    a,b = outer
    W = DISTANCE[1]/4
    need(F(9,10)**2<a and MIDPOINT<F(9,100) and b+MIDPOINT**2<F(13,10), 'g increasing and bounded')
    need(MIDPOINT/(2*F(9,10)**3)<F(1,10) and 3*F(9,10)**2>F(13,100), 'strict convexity of inverse square')

    def bound(u):
        q = 10**7
        r = F(isqrt(u.numerator*q*q//u.denominator),q)
        g = (r-MIDPOINT)**2-W
        need(r>MIDPOINT and g>0,'positive lower denominator')
        return 1/g**2

    fa,fb = map(bound,outer)
    mass = 9-2*(W+MIDPOINT**2)
    need(fa>fb and 8*a<mass<8*b,'decreasing endpoint chord with necessary outer mass')
    need(1/W**2-2*(8*fa+(fb-fa)*(mass-8*a)/(b-a))>STIFFNESS,'strict inner-edge curvature')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=ROOT)
    parser.add_argument('--record',action='store_true')
    args = parser.parse_args()
    ctx = setup(args.repo.resolve())
    stiffness(ctx[2])
    tree_path = HERE/'cover.txt'
    proof = None if args.record else tree_path.read_text(encoding='ascii').strip()
    tags,max_depth = cover(ctx,proof)
    report = {'passed':True,'normalization_mass':9,'applies_to_every_global_maximizer':True,
              'inner_squared_distance_open_interval':list(map(str,DISTANCE)),
              'inner_angle_cosine_strict_upper':str(ANGLE),'inner_midpoint_axis_projection_strict_upper':str(PARALLEL),
              'inner_squared_half_difference_negative_hessian_lower':str(STIFFNESS),
              'inner_edge_deficit_strict_upper':str(LIMIT), 'nodes':len(tags),'maximum_depth':max_depth,'unresolved':0}
    save = HERE/'verification_report.json'
    out = json.dumps(report,indent=2)+'\n'
    if args.record:
        tree_path.write_text(tags+'\n',encoding='ascii',newline='\n')
        save.write_text(out,encoding='ascii',newline='\n')
    else:
        need(report==read(save),'frozen inner-edge geometry report')
    print(out,end='')


if __name__=='__main__':
    main()
