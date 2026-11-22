"""Exact monotonicity and angular endpoint cover for the outer midpoint estimate."""
from fractions import Fraction as F
from math import isqrt, prod

KAPPA = F(31,32)
EXACT_D = {(1,0,0): F(1), (0,1,0): F(1), (0,0,1): F(-2)}
RATIONAL_BOX_SCALE=10**6


def merge(left,right,scale=1):
    return {key: left.get(key,F(0))+scale*right.get(key,F(0)) for key in left.keys()|right.keys()}


def choose(shape,counts,axes):
    return max(axes,key=lambda axis: F(1,1<<counts[axis])*shape[axis])


def traverse(root, initial, axes, test, halves, demand, stored):
    zero=(0,0,0)
    stack=[(root, initial, zero)]
    tags=[]
    while stack:
        box,pack,splits=stack.pop()
        tag=test(box,pack) if stored is None else stored[len(tags)]
        if tag is None:
            demand(sum(splits)<64,'outer endpoint split depth exceeded')
            tag=str(choose([b-a for a,b in root],splits,axes))
        if tag in ('0','1','2'):
            axis=int(tag)
            demand(axis in axes,'allowed split axis')
            l,r=box[axis]
            mid=(l+r)/2
            left,right=list(box),list(box)
            left[axis]=(l,mid)
            right[axis]=(mid,r)
            child=tuple(v+(i==axis) for i,v in enumerate(splits))
            low,high=halves(pack,axis)
            stack.extend(((tuple(right),high,child),(tuple(left),low,child)))
        else:
            demand(tag in ('P','A') and tag==test(box,pack),'outer endpoint leaf check')
        tags.append(tag)
    demand(stored is None or len(stored)==len(tags),'complete outer endpoint tree')
    return ''.join(tags)


def prove_monotonicity(curved,poly,band,stored):
    """The derivative of C/d-log(d) in t has the sign of d C_t+2 C+2d."""
    derivative={(i,j,k-1): k*v for (i,j,k),v in poly.items() if k}
    sign=merge(curved.multiply(EXACT_D,derivative),poly,2)
    sign=merge(sign,EXACT_D,2)
    sign={key:-v for key,v in sign.items() if v}
    shape=tuple(max(key[i] for key in sign)+1 for i in range(3))
    box=(band,band,(-band[1],-KAPPA*band[0]))
    upper=curved.upper
    pack=upper.old.Packed(shape)
    initial=upper.initial(sign,shape,box)
    def test(box,lane):
        (ul,uh),(vl,vh),(tl,th)=box
        if ul>vh or tl**2<KAPPA**2*ul*vl or th**2>uh*vh:
            return 'P'
        if ((lane+pack.high-upper.BIAS*pack.ones)&pack.high)==pack.high:
            return 'A'
        return None
    return traverse(box,initial,(0,1,2),test,pack.split,upper.need,stored)


def endpoint(poly):
    out={}
    for (i,j,k),a in poly.items():
        key=(2*i+k,2*j+k,0)
        out[key]=out.get(key,F(0))+a*(-KAPPA)**k
    return {key:v for key,v in out.items() if v}


class EndpointPacked:
    """Bounds for degree-elevated Bernstein coefficients in base 2^96 lanes."""
    W=96
    WSCALE=1<<32

    def __init__(self, shape, upper):
        self.upper=upper
        self.shape=shape
        count=prod(shape)
        self.ones=((1<<(self.W*count))-1)//((1<<self.W)-1)
        self.high=self.ones<<(self.W-1)
        self.full=[]
        self.valid=[]
        self.active=[]
        for axis,dim in enumerate(shape):
            stride=prod(shape[axis+1:])
            group=dim*stride
            base=sum(1<<(self.W*i) for i in range(0,count,group))*((1<<(self.W*stride))-1)//((1<<self.W)-1)
            one=[base<<(self.W*j*stride) for j in range(dim)]
            self.full.append([x*((1<<self.W)-1) for x in one])
            self.valid.append([sum(one[:dim-j])*((1<<(self.W-1))-1) for j in range(dim)])
            self.active.append([sum(one[:dim-j]) for j in range(dim)])

    def init(self,poly,box,up=False):
        if up:
            poly={key:-a for key,a in poly.items()}
        short=self.upper.initial(poly,self.shape,box)
        b=short.to_bytes(8*prod(self.shape),'little')
        expanded=int.from_bytes(b''.join(b[i:i+8]+b'\0'*4 for i in range(0,len(b),8)),'little')
        if up:
            bias=self.upper.BIAS*self.ones
            self.upper.need(((self.high+bias-expanded)&self.high)==self.high,
                            'positive upper coefficient input for d or d squared')
            return bias-expanded
        return expanded

    def split_one(self,B,axis,up):
        deg=self.shape[axis]-1
        shift=self.W*prod(self.shape[axis+1:])
        A=B
        left=B&self.full[axis][0]
        right=B&self.full[axis][deg]
        for j in range(1,deg+1):
            add=self.active[axis][j] if up else 0
            A=((A+(A>>shift)+add)>>1)&self.valid[axis][j]
            left|=(A&self.full[axis][0])<<(j*shift)
            right|=A&self.full[axis][deg-j]
        return left,right

    def halves(self,B,axis):
        return tuple(zip(*(self.split_one(B[i],axis,up=(i!=0)) for i in range(3))))

    def compare(self,B,eta,c):
        upper=self.upper
        a=F(upper.old.log_upper_scaled(eta),upper.S)-1-c
        b=1/eta
        upper.need(0<a<2 and 0<b<1,'nonnegative outward tangent coefficients')
        a_up=upper.old.ceilq(a*self.WSCALE)
        b_up=upper.old.ceilq(b*self.WSCALE)
        # Input upper coefficients for D and D² lie below 2^60.
        # C's biased lower coefficients are below 2^62. Multiplication
        # by 2^32 and the positive sum below both stay below bit 95.
        threshold=a_up*B[1]+b_up*B[2]+upper.BIAS*self.WSCALE*self.ones
        return ((self.WSCALE*B[0]+self.high-threshold)&self.high)==self.high


def prove_endpoint(curved,poly,band,c,stored):
    upper=curved.upper
    end_poly=endpoint(poly)
    end_d=endpoint(EXACT_D)
    end_d2=curved.multiply(end_d,end_d)
    shape=tuple(max(key[i] for key in end_poly.keys()|end_d2.keys())+1 for i in range(3))
    denom=RATIONAL_BOX_SCALE
    lo=F(isqrt(int(band[0]*denom**2)),denom)
    hi=F(isqrt(int(band[1]*denom**2))+1,denom)
    upper.need(0<lo*lo<=band[0]<band[1]<=hi*hi,'outward radius enclosures')
    box=((lo,hi),(lo,hi),(F(0),F(1)))
    pack=EndpointPacked(shape,upper)
    initial=(pack.init(end_poly,box),pack.init(end_d,box,True),pack.init(end_d2,box,True))
    def test(box,B):
        (ul,uh),(vl,vh),_=box
        if ul>vh:
            return 'P'
        dlo=ul**2+vl**2+2*KAPPA*ul*vl
        dhi=uh**2+vh**2+2*KAPPA*uh*vh
        if pack.compare(B,(dlo+dhi)/2,c):
            return 'A'
        return None
    return traverse(box,initial,(0,1),test,pack.halves,upper.need,stored)
