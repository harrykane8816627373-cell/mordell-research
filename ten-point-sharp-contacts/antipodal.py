"""Replay the exact zero-compensation cover of opposite outer pairs."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
import json
from math import comb, prod

import cover

model = cover.model
SCALE = model.geometry.SCALE


@lru_cache(None)
def radial_power(a,b):
    """Coefficients of ((x-y)/2)^a ((x+y)/2)^b."""
    terms=[]
    for j in range(a+b+1):
        value=sum((-1)**i*comb(a,i)*comb(b,j-i) for i in range(max(0,j-b),min(a,j)+1))
        if value:
            terms.append((a+b-j,j,F(value,2**(a+b))))
    return terms


def regular_coordinates(poly):
    """Substitute u=p², v=q², t=-pq+e, where p=(x-y)/2, q=(x+y)/2."""
    result={}
    for (i,j,k),coefficient in poly.items():
        for m in range(k+1):
            term=coefficient*((-1)**(k-m)*comb(k,m))
            for a,b,c in radial_power(2*i+k-m,2*j+k-m):
                key=(a,b,m)
                result[key]=result.get(key,model.ZERO)+term*model.exact(c)
    for power,value in result.items():
        if power[1]%2:
            model.need(value.lo<=0<=value.hi,'outer kernel is symmetric under swapping u and v')
    return {power:value for power,value in result.items() if power[1]%2==0}


def falling(n,k):
    return prod(range(n-k+1,n+1)) if k<=n else 0


class Packed(cover.Packed):
    def __init__(self,shape):
        super().__init__(shape)
        self.log_masks={}
        for i in range(5):
            for j in range(3):
                mask=sum((x&e)*falling(k,i)*falling(l,j) for k,x in enumerate(self.axisones[0])
                         for l,e in enumerate(self.axisones[2]) if k>=i and l>=j)
                self.log_masks[i,j]=mask,falling(self.degrees[0],i)*falling(self.degrees[2],j)

    def dominates(self,control,box):
        (xl,xh),_,(el,eh)=box
        a=xh*xh-2*el
        origin=xl*xl-2*el
        w=xh-xl
        movement={(1,0):2*xl*w,(2,0):w*w,(0,1):-2*(eh-el)}
        q=F(cover.engine.log_upper_lanes(a),cover.S)+F(3,2)
        threshold={key:2*(q*origin-a)*s for key,s in movement.items()}
        threshold[0,0]=(q*origin-2*a)*origin+a*a/2
        for i,j in combinations_with_replacement(movement,2):
            key=(i[0]+j[0],i[1]+j[1])
            threshold[key]=threshold.get(key,F(0))+q*movement[i]*movement[j]*(1+int(i!=j))
        limit=cover.BIAS*self.ones
        lower=upper=cover.BIAS
        for power,coefficient in threshold.items():
            mask,denominator=self.log_masks[power]
            lane=cover.packed_engine.ceilq(cover.S*coefficient/denominator)
            limit+=mask*lane
            lower+=min(0,lane*denominator)
            upper+=max(0,lane*denominator)
        model.need(0<=lower<=upper<1<<(cover.W-2),'quadratic logarithm lane bounds')
        return ((control+self.high-limit)&self.high)==self.high


def physical_box(box):
    (xl,xh),(yl,yh),(el,eh)=box
    pl,ph=(xl-yh)/2,(xh-yl)/2
    ql,qh=(xl+yl)/2,(xh+yh)/2
    model.need(pl>0 and ql>0,'positive outer moduli')
    return ((pl**2,ph**2),(ql**2,qh**2),(-ph*qh+el,-pl*ql+eh)),ph*qh


def evaluate(poly,point):
    values=[model.exact(x) for x in point]
    return sum((v*prod((model.sharp.ball_power(x,i) for x,i in zip(values,key)),start=model.ONE)
                for key,v in poly.items()),model.ZERO)


def derivative(poly,axis):
    return {p[:axis]+(p[axis]-1,)+p[axis+1:]:v*p[axis] for p,v in poly.items() if p[axis]}


def ball(span):
    return model.Ball(model.exact(span[0]).lo,model.exact(span[1]).hi)


def centered_value(poly,span):
    middle=sum(span)/2
    x=model.exact(middle)
    delta=ball((span[0]-middle,span[1]-middle))
    coefficients=[model.ZERO]*(1+max(p[0] for p in poly))
    for (k,),a in poly.items():
        for j in range(k+1):
            coefficients[j]+=a*comb(k,j)*model.sharp.ball_power(x,k-j)
    return sum((a*model.sharp.ball_power(delta,i) for i,a in enumerate(coefficients)),model.ZERO)


def nonnegative(poly,box):
    packed,control=cover.initial_data(poly,box)
    return ((control+packed.high-cover.BIAS*packed.ones)&packed.high)==packed.high


def contact_tubes(poly,nodes):
    d=cover.Polynomial({(2,0,0):model.ONE,(0,0,1):model.exact(-2)})
    normal=(d*cover.Polynomial(derivative(poly,2))+4*cover.Polynomial(poly)+2*d**2).terms
    axial={(i,):v for (i,j,k),v in poly.items() if j==k==0}
    transverse={(i,j-2):v for (i,j,k),v in poly.items() if j>=2 and k==0}
    axial2=derivative(derivative(axial,0),0)
    axial3=derivative(axial2,0)
    tubes,radii=[],[]
    for node in nodes[1:]:
        anchor=2*model.sharp.sym.sqrt_ball(node,model.Ball,SCALE)
        middle=F((anchor.lo+anchor.hi)//2,SCALE)
        for power in range(8,19):
            radius=F(1,2**power)
            span=(middle-radius,middle+radius)
            window=ball(span)
            x0=model.exact(middle)
            center=centered_value(axial2,(middle,middle))-(24*model.sharp.log_ball(x0)+14)*x0.square()
            slope=centered_value(axial3,span)-(48*model.sharp.log_ball(window)+52)*window
            curvature=center+slope*ball((-radius,radius))
            near=(span,(F(0),radius),(F(0),radius))
            if curvature.lo<=0 or not nonnegative(transverse,near[:2]) or not nonnegative(normal,near):
                continue
            model.need(model.exact(span[0]).hi<anchor.lo<=anchor.hi<model.exact(span[1]).lo and
                       span[0]**2>2*radius,'tube contains the exact anchor and has positive distance')
            tubes.append(near)
            radii.append(str(radius))
            break
        else:
            raise ArithmeticError('strict Gram-boundary contact tube')
    return tubes,radii


def verify(kernel,report,node_limit=100000):
    cls=5
    _,nodes,_,_,spans=kernel
    band=model.BANDS[1]
    angle_budget=1+spans[4][1]
    raw=cover.numerator(cls,kernel)
    polynomial=regular_coordinates(raw)
    x,y,e=F(21,10),F(1,20),F(1,100)
    p,q=(x-y)/2,(x+y)/2
    old,new=evaluate(raw,(p*p,q*q,-p*q+e)),evaluate(polynomial,(x,y,e))
    model.need(old.lo<=new.hi and new.lo<=old.hi,'exact check of the coordinate expansion')
    lower,upper=F(943,1000),F(1128,1000)
    model.need(lower**2<band[0]<band[1]<upper**2 and 0<angle_budget<1,'outer modulus and angle bounds')
    root=((2*lower,2*upper),(F(0),upper-lower),(F(0),angle_budget*band[1]))
    tubes,radii=contact_tubes(polynomial,nodes)
    boxes=[tuple(tuple(map(F,p)) for p in row['cover_box']) for row in report['local_boxes'] if row['class']==cls]
    shape=tuple(max(p[i] for p in polynomial)+1 for i in range(3))
    packed=Packed(shape)
    control=packed.initial({p:F(v.lo,SCALE) for p,v in polynomial.items()},root)

    def exclude(box):
        (u,v,_),product=physical_box(box)
        return (u[1]<band[0] or v[1]<band[0] or u[0]>band[1] or v[0]>band[1] or
                box[2][0]>angle_budget*product)

    def inside(box):
        return any(cover.nested(box,q) for q in tubes) or any(cover.nested(physical_box(box)[0],q) for q in boxes)

    counts=cover.split_cover(packed,root,control,exclude,inside,packed.dominates,node_limit)
    model.need(counts['unresolved']==0,'complete opposite outer-pair cover')
    return {'class':cls,'tensor_shape':list(shape),'contact_tube_radii':radii,**counts}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    contact_report,kernel=model.exact_kernel()
    result=verify(kernel,contact_report)
    path=model.HERE/'outer_antipodal_cover_report.json'
    encoded=json.dumps(result,indent=2)+'\n'
    if args.record_report:
        path.write_bytes(encoded.encode())
    else:
        model.need(json.loads(path.read_bytes())==result,'matches frozen opposite-pair replay')
    print(encoded,end='')


if __name__=='__main__':
    main()
