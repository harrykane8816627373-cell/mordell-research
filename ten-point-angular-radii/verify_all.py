"""Verify marked-radius exclusions from angular-class row bounds."""

import argparse
import importlib.util
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
spec = importlib.util.spec_from_file_location('n10_row_angular', ROOT/'ten-point-angular-upper/verify_all.py')
angular = importlib.util.module_from_spec(spec)
spec.loader.exec_module(angular)
curved, upper, need = angular.curved, angular.upper, angular.need

CLASSES = ((0, 1), (1, 2, 3, 4, 5))
COUNTS = ((1, 8), (2, 2, 2, 2, 1))
BANDS = angular.BANDS
NEW_BANDS = ((F(84, 1000), F(159, 1000)), (F(89, 100), F(1271, 1000)))
REGIONS = {
    'inner-low': (0, BANDS[0][0], NEW_BANDS[0][0]),
    'inner-high': (0, NEW_BANDS[0][1], BANDS[0][1]),
    'outer-low': (1, BANDS[1][0], NEW_BANDS[1][0]),
    'outer-high': (1, NEW_BANDS[1][1], BANDS[1][1]),
}
DISTANCE = {(1,0,0): F(1), (0,1,0): F(1), (0,0,1): F(-2)}
ONE, U, V, T = ({(0,0,0): F(1)}, {(1,0,0): F(1)}, {(0,1,0): F(1)}, {(0,0,1): F(1)})


def combine(*terms):
    output = {}
    for weight, poly in terms:
        for key, coeff in poly.items():
            output[key] = output.get(key, F(0)) + weight*coeff
    return {key: value for key, value in output.items() if value}


def multiply(a, b):
    return curved.multiply(a, b)


def ward_numerators():
    u2, u3 = multiply(U,U), multiply(multiply(U,U),U)
    h2 = combine((2,multiply(T,T)), (-1,multiply(U,V)))
    h3 = combine((4,multiply(T,multiply(T,T))), (-3,multiply(multiply(U,V),T)))
    numerators = [multiply(DISTANCE, combine((1,V),(-1,ONE),(F(1,9),U))),
                  multiply(DISTANCE, combine((1,T),(F(1,9),U))),
                  combine((1,U),(-1,T),(-F(5,9),multiply(U,DISTANCE))),
                  multiply(DISTANCE, combine((1,multiply(V,T)),(F(1,9),u2)))]
    for power, harmonic, ward in ((u2,h2,F(17,10)), (u3,h3,F(8,5))):
        full = combine((1,multiply(combine((1,V),(-ward,ONE)),harmonic)),
                       (F(1,9),multiply(combine((1,U),(-ward,ONE)),power)))
        numerators.append(multiply(DISTANCE, full))
    return numerators


def at(poly, u, v, t):
    return sum(a*u**i*v**j*t**k for (i,j,k),a in poly.items())


def crosscheck(numerators):
    for u,v,t in ((F(1,9),F(1),F(-1,4)), (F(11,10),F(9,8),F(1,4)),
                  (F(1,8),F(2,17),F(-11,100))):
        d=u+v-2*t
        expected=(v-1+u/9,t+u/9,(u-t)/d-5*u/9,v*t+u**2/9,
                  (v-F(17,10))*(2*t**2-u*v)+(u-F(17,10))*u**2/9,
                  (v-F(8,5))*(4*t**3-3*u*v*t)+(u-F(8,5))*u**3/9)
        need(tuple(at(poly,u,v,t)/d for poly in numerators)==expected, 'marked Ward row expansion')


def setup():
    data=angular.read(curved.HERE/'kernel.json')
    report=angular.read(angular.HERE/'verification_report.json')
    need(report['passed'] and report['unresolved']==0 and report['normalization_mass']==9,
         'accepted angular-class deficit bounds')
    need(tuple(tuple(c['class']) for c in report['classes'])==angular.CLASSES
         and tuple(c['weight'] for c in report['classes'])==angular.WEIGHTS, 'fixed classes and counts')
    floors=tuple(-F(c['deficit_floor']) for c in report['classes'])
    gap=upper.hx(data['B'])-upper.log.log_lower(upper.actual_root_lower())
    total=sum(w*c for w,c in zip(angular.WEIGHTS,floors))
    caps=tuple(gap+total-sum(w*floors[category] for w,category in zip(counts,classes))
               for counts,classes in zip(COUNTS,CLASSES))
    ward,force=upper.corrections(data)
    need(bool(force), 'cleared distance denominator')
    base,shape=upper.polynomial(data,ward,force)
    params=curved.parsed(data)
    band_polys={(i,j):curved.root_polynomial(base,shape,params,i,j) for i in range(2) for j in range(2)}
    for i,j in band_polys:
        curved.check_root(base,band_polys[i,j],params,i,j)
    numerators=ward_numerators()
    crosscheck(numerators)
    return band_polys,shape,numerators,caps


def outside(category, box):
    i,j,k=angular.CLASSES[category]
    a,b,t=box
    sorted_pair=((min(a[0],b[0]),min(a[1],b[1])), (max(a[0],b[0]),max(a[1],b[1])), t)
    return angular.outside(i,j,k,sorted_pair)


def cover(poly,shape,packer,box,category,tree):
    stack=[(box,upper.initial(poly,shape,box),0)]
    counts=Counter()
    sequence=[]
    depth_max=0
    i,j,_=angular.CLASSES[category]
    widths=sorted([BANDS[i][1]-BANDS[i][0],BANDS[j][1]-BANDS[j][0]])
    if category==1 and box[0][0]>BANDS[0][1]:
        widths.reverse()
    widths.append(F(1) if category==1 else F(2) if category else F(1,6))
    for item in tree if tree is not None else iter(lambda: bool(stack),False):
        need(bool(stack), 'remaining box for stored node')
        box,packed,depth=stack.pop()
        depth_max=max(depth_max,depth)
        if tree is not None:
            tag=item
        elif outside(category,box):
            tag='P'
        elif packer.dominates(packed,box,F(0)):
            tag='A'
        elif packer.curved_dominates(packed,box,F(0)):
            tag='Q'
        else:
            if depth>=58:
                import math
                point=[sum(x)/2 for x in box]
                d=point[0]+point[1]-2*point[2]
                raise ArithmeticError(f'unresolved class {category}, {[float(x) for x in point]}, margin {float(at(poly,*point)/d)-math.log(d)}')
            axis=max(range(3),key=lambda a: (box[a][1]-box[a][0])/widths[a])
            tag=str(axis)
        if tag in '012':
            stack.extend(curved.split(packer,packed,box,depth,int(tag)))
        elif tree is not None:
            accepted=(tag=='P' and outside(category,box) or tag=='A' and packer.dominates(packed,box,F(0))
                      or tag=='Q' and packer.curved_dominates(packed,box,F(0)))
            need(accepted, f'valid {tag} leaf for class {category} in {box}')
        counts[tag]+=1
        sequence.append(tag)
    need(not stack and len(sequence)==2*sum(counts[c] for c in '012')+1,'complete subtree')
    return ''.join(sequence),depth_max


def region(name,built,build=None,progress=False):
    mark,start,end = REGIONS[name]
    polynomials,shape,numerators,caps=built
    packer=curved.Packed(shape)
    filename=HERE/f'{name}.json'
    rows=angular.read(Path(build)/f'{name}.json' if build is not None else filename)
    cursor=start
    nodes, maximum=0,0
    min_margin=None
    for entry in rows:
        lo,hi=map(F,entry['range'])
        need(lo==cursor and lo<hi<=end,'gapless cover of the marked-radius interval')
        cursor=hi
        weights=tuple(map(F,entry['multipliers']))
        floors=tuple(map(F,entry['floor']))
        classes=CLASSES[mark]
        need(len(weights)==len(numerators) and len(floors)==len(classes),'complete row multipliers')
        margin=sum(n*c for n,c in zip(COUNTS[mark],floors))-caps[mark]
        need(margin>F(1,1000000),'strict regional exclusion margin')
        min_margin=margin if min_margin is None else min(min_margin,margin)
        trees=[]
        for k,category in enumerate(classes):
            i,j,step=angular.CLASSES[category]
            neighbor=i+j-mark
            base=polynomials[mark,neighbor]
            poly=combine((1,base),(-floors[k],DISTANCE),*[(-v,p) for v,p in zip(weights,numerators)])
            root_intervals=[box[2] for box in angular.root_boxes(i,j,step)]
            stored=[None]*len(root_intervals) if build is not None else entry['trees'][k]
            need(len(stored)==len(root_intervals),'all signed root boxes')
            class_trees=[]
            for t_interval,tree in zip(root_intervals,stored):
                out,depth=cover(poly,shape,packer,((lo,hi),BANDS[neighbor],t_interval),category,tree)
                class_trees.append(out)
                nodes+=len(out)
                maximum=max(maximum,depth)
            trees.append(class_trees)
        entry['trees']=trees
        if progress:
            print(name, str(hi), nodes, flush=True)
    need(cursor==end,'last marked-radius interval endpoint')
    if build is not None:
        filename.write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8',newline='\n')
    return {'interval':[str(start),str(end)],'cells':len(rows),'nodes':nodes,'maximum_depth':maximum,
            'strict_margin_decimal_display':float(min_margin)}


def check_inputs():
    manifest=angular.read(HERE/'source_manifest.json')
    for base,group in ((HERE,'files'),(ROOT,'dependencies')):
        for name,digest in manifest[group].items():
            need(sha256((base/name).read_bytes()).hexdigest()==digest,'fixed input '+name)
    angular.check_domains()
    upper.old.selftest()
    curved.curvature_selftest()


def inner_curvature():
    """Rational check of strict concavity in the inner squared half-difference."""
    a,b=NEW_BANDS[1]
    midpoint=F(239,2000)
    wmax=F(11,80)
    scale=10**6
    def upper_inverse(u):
        r=F(isqrt(u.numerator*scale**2//u.denominator),scale)
        g=(r-midpoint)**2-wmax
        need(r>midpoint and g>0, 'positive lower bound for the complex quadratic denominator')
        return 1/g**2
    fa,fb=map(upper_inverse,(a,b))
    need(F(9,10)**2<a and midpoint<F(3,25) and b+midpoint**2<F(13,10),
         'g prime greater than 13/15 and g below 13/10')
    need(midpoint/(2*F(9,10)**3)<F(1,10) and 3*F(13,15)**2>F(13,100),
         'inverse-square denominator is convex over the outer band')
    mass_floor=9-2*(wmax+midpoint**2)
    need(fa>fb and 8*a<mass_floor<8*b, 'decreasing endpoint chord and outer mass lower bound')
    total=8*fa+(mass_floor-8*a)*(fb-fa)/(b-a)
    need(1/wmax**2-2*total>15, 'strict inner squared half-difference concavity')
    return '15'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path)
    parser.add_argument('--region',choices=REGIONS)
    parser.add_argument('--progress',action='store_true')
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    if args.build is None:
        check_inputs()
    built=setup()
    if args.region:
        print(json.dumps(region(args.region,built,args.build,args.progress),indent=2))
        return
    need(args.build is None,'build one region at a time')
    regions={name:region(name,built,progress=args.progress) for name in REGIONS}
    result={'passed':True,'applies_to_every_global_maximizer':True,'normalization_mass':9,
            'exact_inner_count':2,'exact_outer_count':8,
            'inner_squared_radius_open_interval':list(map(str,NEW_BANDS[0])),
            'outer_squared_radius_open_interval':list(map(str,NEW_BANDS[1])),
            'inner_squared_half_difference_negative_hessian_lower':inner_curvature(),
            'nodes':sum(r['nodes'] for r in regions.values()),'unresolved':0,'regions':regions}
    path=HERE/'verification_report.json'
    output=json.dumps(result,indent=2)+'\n'
    if args.record_report:
        path.write_text(output,encoding='utf-8',newline='\n')
    else:
        need(result==angular.read(path),'saved proof report')
    print(output,end='')


if __name__=='__main__':
    main()
