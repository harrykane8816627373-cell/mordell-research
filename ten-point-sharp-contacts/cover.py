"""Test zero-compensation coverage for each cyclic class of the sharp contact fit."""

import argparse
from fractions import Fraction as F
from itertools import combinations
from math import lcm,prod
import json

import verify_all as model

engine = model.sharp.load('n12_packed_cover',model.ROOT/'twelve-point-sharp-global/verify_all.py')
packed_engine=engine.cover
S=packed_engine.S
BIAS=packed_engine.BIAS
W=packed_engine.W
ZERO=model.ZERO
CLASSES=[(0,0),(0,1)]+[(1,1)]*4
DMIN=F(2,5)


class Polynomial:
    def __init__(self,terms):
        self.terms={power:value for power,value in terms.items() if value.lo or value.hi}

    @staticmethod
    def constant(value):
        return Polynomial({(0,0,0):model.exact(value)})

    @staticmethod
    def coerce(value):
        return value if isinstance(value,Polynomial) else Polynomial.constant(value)

    def __add__(self,rhs):
        terms=self.terms.copy()
        for power,value in self.coerce(rhs).terms.items():
            terms[power]=terms.get(power,ZERO)+value
        return Polynomial(terms)

    __radd__=__add__

    def __neg__(self):
        return Polynomial({power:-value for power,value in self.terms.items()})

    def __sub__(self,rhs):
        return self+-self.coerce(rhs)

    def __rsub__(self,lhs):
        return self.coerce(lhs)+-self

    def __mul__(self,rhs):
        terms={}
        for a,x in self.terms.items():
            for b,y in self.coerce(rhs).terms.items():
                power=tuple(i+j for i,j in zip(a,b))
                terms[power]=terms.get(power,ZERO)+x*y
        return Polynomial(terms)

    __rmul__=__mul__

    def __truediv__(self,rhs):
        return self*(model.ONE/model.exact(rhs))

    def __pow__(self,k):
        result=self.constant(1)
        for _ in range(k):
            result*=self
        return result


def numerator(cls,kernel):
    weights,nodes,transforms,budget,_=kernel
    axes=[Polynomial({tuple(int(i==j) for i in range(3)):model.ONE}) for j in range(3)]
    u,v,t=axes
    unity=Polynomial.constant(1)
    d=u+v-2*t
    base=unity*(budget/45-1)+(u+v)*F(5,9)+model.rigidity_term(u,v,nodes)
    columns=model.feature_columns(u,v,t,unity,CLASSES[cls],cls-1,nodes,transforms,True)
    value=base*d**2+sum((column*weight for column,weight in zip(columns,weights)),Polynomial.constant(0))
    # Check the denominator-cleared evaluation against the Jet-based one.
    points=[(F(1,9),F(2,15),F(-1,10)),(F(7,50),F(11,10),F(-2,5)),(F(1),F(6,5),F(2,3))]
    point=points[min(cls,2)]
    u,v,t=map(model.exact,point)
    physical=dict(u=u,v=v,t=t,layers=CLASSES[cls],**{'class':cls})
    expected=(budget/45-1+(u+v)*model.exact(F(5,9))+
              model.rigidity_term(model.Jet(u),model.Jet(v),nodes).v+
              model.dot(weights,[jet.v for jet in model.columns(physical,nodes,transforms)]))*(u+v-2*t).square()
    actual=sum((a*model.sharp.ball_power(u,i)*model.sharp.ball_power(v,j)*model.sharp.ball_power(t,k)
                for (i,j,k),a in value.terms.items()),ZERO)
    model.need(actual.lo<=expected.hi and expected.lo<=actual.hi,'expanded numerator crosscheck')
    return value.terms


class Packed(packed_engine.Packed):
    def __init__(self,shape):
        super().__init__(shape)
        self.degrees=tuple(n-1 for n in shape)
        self.linear=[sum(i*m for i,m in enumerate(masks)) for masks in self.axisones]
        self.square=[sum(i*(i-1)*m for i,m in enumerate(masks)) for masks in self.axisones]
        self.cross={(i,j):sum(x*y*(a&b) for x,a in enumerate(self.axisones[i])
                             for y,b in enumerate(self.axisones[j])) for i,j in combinations(range(len(shape)),2)}

    def initial(self,coefficients,box):
        denominator=lcm(*(a.denominator for a in coefficients.values()))
        tensor=[0]*prod(self.shape)
        for power,a in coefficients.items():
            index=0
            for dimension,exponent in zip(self.shape,power):
                index=index*dimension+exponent
            tensor[index]=a.numerator*(denominator//a.denominator)
        for axis,dimension in enumerate(self.shape):
            matrix,divisor=packed_engine.transform(dimension-1,*box[axis])
            stride=prod(self.shape[axis+1:])
            group=dimension*stride
            transformed=[0]*len(tensor)
            for base in range(0,len(tensor),group):
                for offset in range(stride):
                    values=tensor[base+offset:base+offset+group:stride]
                    for i,row in enumerate(matrix):
                        transformed[base+i*stride+offset]=sum(x*y for x,y in zip(row,values))
            tensor=transformed
            denominator*=divisor
        packed=0
        for value in reversed(tensor):
            lane=S*value//denominator+BIAS
            model.need(0<=lane<1<<(W-2),'initial coefficient lane bounds')
            packed=(packed<<W)|lane
        return packed

    def dominates(self,control,box):
        (ul,uh),(vl,vh),(tl,th)=box
        a=uh+vh-2*tl
        q=F(engine.log_upper_lanes(a),S)+F(3,2)
        d0=ul+vl-2*tl
        widths=(uh-ul,vh-vl,-2*(th-tl))
        return self.upper_quadratic(control,q,-2*a,a*a/2,d0,widths)

    def upper_quadratic(self,control,q,linear,constant,origin,widths):
        bias=BIAS+packed_engine.ceilq(S*((q*origin+linear)*origin+constant))
        bound=bias*self.ones
        lo=hi=bias

        def add(coefficient,mask,denominator):
            nonlocal bound,lo,hi
            value=packed_engine.ceilq(S*coefficient/denominator)
            bound+=value*mask
            lo+=min(0,value*denominator)
            hi+=max(0,value*denominator)

        for i,w in enumerate(widths):
            n=self.degrees[i]
            add((2*q*origin+linear)*w,self.linear[i],n)
            add(q*w*w,self.square[i],n*(n-1))
        for (i,j),mask in self.cross.items():
            add(2*q*widths[i]*widths[j],mask,self.degrees[i]*self.degrees[j])
        model.need(0<=lo<=hi<1<<(W-2),'quadratic threshold lane bounds')
        return ((control+self.high-bound)&self.high)==self.high


def roots(cls,spans):
    i,j=CLASSES[cls]
    u,v=model.BANDS[i],model.BANDS[j]
    if cls==0:
        return [(u,v,(-u[1],F(0)))]
    if cls==1:
        r=F(9,20)
        model.need(r*r>u[1]*v[1],'mixed dot support')
        return [(u,v,(-r,F(0))),(u,v,(F(0),r))]
    lo,hi=spans[cls-1]
    a=u[1 if lo<0 else 0]*lo
    b=u[1 if hi>0 else 0]*hi
    return [(u,v,(a,F(0))),(u,v,(F(0),b))] if a<0<b else [(u,v,(a,b))]


def outside(box,cls,spans,inner):
    (ul,uh),(vl,vh),(tl,th)=box
    if ul>vh or max(tl,-th,0)**2>uh*vh or uh+vh-2*tl<DMIN:
        return True
    if cls==0:
        dlo,dhi,cos=inner
        return (uh+vh-2*tl<dlo or ul+vl-2*th>dhi or tl>=0 or tl*tl<cos*cos*ul*vl)
    if cls==1:
        return False
    lo,hi=spans[cls-1]
    if lo>=0 and (th<0 or th*th<lo*lo*ul*vl):
        return True
    if lo<0 and th<0 and th*th>lo*lo*uh*vh:
        return True
    if hi<=0 and (tl>0 or tl*tl<hi*hi*ul*vl):
        return True
    return hi>0 and tl>0 and tl*tl>hi*hi*uh*vh


def inner_bounds():
    report=model.read('ten-point-inner-edge/verification_report.json')
    model.need(report['passed'] and report['unresolved']==0 and report['normalization_mass']==9 and
               report['applies_to_every_global_maximizer'],'certified inner pair bounds')
    return (*map(F,report['inner_squared_distance_open_interval']),F(report['inner_angle_cosine_strict_upper']))


def split_cover(packed,root,control,exclude,local,accept,node_limit):
    counts={key:0 for key in ('visited','split','geometry','local','bernstein','unresolved')}
    deepest=0
    scales=[b-a for a,b in root]
    stack=[(root,control,0)]
    while stack:
        box,patch,depth=stack.pop()
        deepest=max(deepest,depth)
        counts['visited']+=1
        if exclude(box):
            counts['geometry']+=1
        elif local(box):
            counts['local']+=1
        elif accept(patch,box):
            counts['bernstein']+=1
        elif depth>=90 or counts['visited']>=node_limit:
            counts['unresolved']+=1
            if counts['visited']>=node_limit:
                counts['unresolved']+=len(stack)
                break
        else:
            axis=max(range(len(box)),key=lambda i:(box[i][1]-box[i][0])/scales[i])
            a,b=box[axis]
            m=(a+b)/2
            left,right=packed.split(patch,axis)
            stack.append((box[:axis]+((m,b),)+box[axis+1:],right,depth+1))
            stack.append((box[:axis]+((a,m),)+box[axis+1:],left,depth+1))
            counts['split']+=1
    if counts['unresolved']==0:
        model.need(counts['visited']==1+2*counts['split'] and
                   counts['visited']==sum(v for k,v in counts.items() if k!='visited'),'complete cover accounting')
    return {'depth':deepest,**counts}


def initial_data(poly,box):
    shape=tuple(1+max(p[i] for p in poly) for i in range(len(box)))
    packed=Packed(shape)
    model.need(len(box)!=3 or box[2][1]<=0 or box[2][0]>=0, 'fixed dot-product sign in a root box')
    negative=len(box)==3 and box[2][1]<=0
    lower={p:F(a.hi if negative and p[2]%2 else a.lo,model.geometry.SCALE) for p,a in poly.items()}
    return packed,packed.initial(lower,box)


def nested(box,near):
    return all(c<=a<=b<=d for (a,b),(c,d) in zip(box,near))


def inner_cover(kernel,report,node_limit):
    """Use monotonicity toward the Gram face to handle the only inner edge."""
    raw=numerator(0,kernel)
    source=Polynomial(raw)
    dt=Polynomial({(a,b,c-1):v*c for (a,b,c),v in raw.items() if c})
    d=Polynomial({(1,0,0):model.ONE,(0,1,0):model.ONE,(0,0,1):model.exact(-2)})
    normal=(d*dt+4*source+2*d**2).terms
    dlo,dhi,cos=inner_bounds()
    band=model.BANDS[0]
    projected_dhi=dhi+2*(1+cos)*band[1]
    root=(band,band,(-band[1],cos*band[0]))
    packed,control=initial_data(normal,root)
    normal_count=split_cover(packed,root,control,lambda b:outside(b,0,{},(dlo,projected_dhi,cos)),
                             lambda b:False,lambda p,b:((p+packed.high-BIAS*packed.ones)&packed.high)==packed.high,
                             node_limit)

    face={}
    for (a,b,c),v in raw.items():
        key=(2*a+c,2*b+c)
        face[key]=face.get(key,ZERO)+(-v if c%2 else v)
    edge=(F(289,1000),F(399,1000))
    model.need(edge[0]**2<band[0]<band[1]<edge[1]**2 and edge[0]>F(1,4),'square-root enclosure')
    root=(edge,edge)
    packed,control=initial_data(face,root)
    near=next(row for row in report['local_boxes'] if row['class']==0)
    local=tuple(tuple(map(F,b)) for b in near['cover_box'])

    def exclude(box):
        (pl,ph),(ql,qh)=box
        return (pl>qh or ph**2<band[0] or qh**2<band[0] or pl**2>band[1] or ql**2>band[1] or
                (ph+qh)**2<dlo or (pl+ql)**2>projected_dhi)

    def inside(box):
        (pl,ph),(ql,qh)=box
        return nested(((pl**2,ph**2),(ql**2,qh**2),(-ph*qh,-pl*ql)),local)

    def accept(control,box):
        (pl,ph),(ql,qh)=box
        a=ph+qh
        la=F(engine.log_upper_lanes(a),S)
        q=(12*la+7)*a*a
        first=(8*la+2)*a**3
        value=2*a**4*la
        return packed.upper_quadratic(control,q,first-2*q*a,value-first*a+q*a*a,
                                      pl+ql,(ph-pl,qh-ql))

    face_count=split_cover(packed,root,control,exclude,inside,accept,node_limit)
    return {'class':0,'normal':normal_count,'face':face_count,
            'unresolved':normal_count['unresolved']+face_count['unresolved']}


def cover_class(cls,kernel,report,node_limit):
    _,_,_,_,spans=kernel
    poly=numerator(cls,kernel)
    inner=inner_bounds()
    local=[tuple(tuple(map(F,b)) for b in row['cover_box']) for row in report['local_boxes'] if row['class']==cls]
    boxes=roots(cls,spans)
    counts=[]
    for root in boxes:
        packed,control=initial_data(poly,root)
        counts.append(split_cover(packed,root,control,lambda b:outside(b,cls,spans,inner),
                                  lambda b:any(nested(b,near) for near in local),packed.dominates,node_limit))
    return {'class':cls,'roots':counts,'tensor_shape':list(packed.shape),
            'unresolved':sum(row['unresolved'] for row in counts)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--explore',type=int,nargs='+',choices=range(1,6),
                      help='run an exploratory cover of the selected pair classes')
    mode.add_argument('--outer-step-two',action='store_true',
                      help='replay the complete cover for outer pairs two positions apart')
    parser.add_argument('--node-limit',type=int,default=100000)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    if args.node_limit<1 or (args.explore and args.record_report):
        parser.error('positive node limit required and exploration cannot record the proof report')
    report,kernel=model.exact_kernel()
    if args.explore:
        for cls in args.explore:
            print(json.dumps(cover_class(cls,kernel,report,args.node_limit)),flush=True)
        return
    if args.outer_step_two:
        result=cover_class(3,kernel,report,args.node_limit)
        path=model.HERE/'outer_step_two_cover_report.json'
    else:
        result=inner_cover(kernel,report,args.node_limit)
        path=model.HERE/'inner_cover_report.json'
    model.need(result['unresolved']==0,'complete zero-compensation pair cover')
    encoded=json.dumps(result,indent=2)+'\n'
    if args.record_report:
        path.write_bytes(encoded.encode())
    else:
        model.need(json.loads(path.read_bytes())==result,'matches frozen inner edge replay')
    print(encoded,end='')


if __name__=='__main__':
    main()
