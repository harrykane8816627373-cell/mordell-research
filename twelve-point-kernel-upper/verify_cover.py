"""Standard-library exact Bernstein replay using packed unbounded integers.
Discovery supplies only a subdivision tree; every input, split, discard and
log-tangent domination is recomputed with integers/Fraction. No float bounds.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,lcm,prod
from functools import lru_cache
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parent.parent;S=1<<40;BIAS=1<<60;W=64

def need(ok,label):
 if not ok:raise ArithmeticError(label)
def ceilq(q):return -((-q.numerator)//q.denominator)
def hx(a):return F.from_float(float.fromhex(a))
def polymul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def features(data,k):
 s=data['feature_order']
 if k==0:
  nodes=list(map(hx,data['radial_nodes']));means=list(map(hx,data['radial_moment_centers']));q=[F(1)]
  for a in nodes:q=polymul(q,[-a,F(1)])
  M=[[F(0)]*s for _ in range(2*len(nodes)+1)]
  for j in range(1,s):M[j][j-1]=1;M[0][j-1]=-means[j-1]
  for i,a in enumerate(polymul(q,q)):M[i][-1]=10*a
  return M
 E=[[hx(a) for a in row] for row in data['transforms'][k-1]]
 return [[sum((F(comb(j,i)*(-1)**(j-i))*E[j][r] for j in range(i,s)),F(0)) for r in range(len(E[0]))] for i in range(s)]

# For each k, G_k(u,v)=F_k(u)^T Q_k F_k(v).  The stored polynomial is
# h = n(u+v)/(2(n-1)) + 2b/(n-1) - sum_k H_k G_k(u,v)
#     - [sum_k u^k G_k(u,u) + sum_k v^k G_k(v,v)]/(2(n-1)).
# The diagonal identity H_k(u,u,u)=u^k gives every correction exactly.
def polynomial(data):
 modes=data['harmonics_through'];degree=modes+2*(data['feature_order']-1);shape=(degree+1,degree+1,modes+1);C={};H=[{(0,0,0):1},{(0,0,1):1}]
 def add(index,value):C[index]=C.get(index,F(0))+value
 for k in range(2,modes+1):
  q={}
  for (a,b,c),v in H[-1].items():q[(a,b,c+1)]=q.get((a,b,c+1),0)+2*v
  for (a,b,c),v in H[-2].items():q[(a+1,b+1,c)]=q.get((a+1,b+1,c),0)-v
  H.append(q)
 n=data['n'];m=n-1;B=hx(data['B']);b0=(B-F(n*m,2))/n
 add((0,0,0),2*b0/m);add((1,0,0),F(n,2*m));add((0,1,0),F(n,2*m))
 for k,raw in enumerate(data['Q_blocks']):
  Q=[[hx(a) for a in row] for row in raw];d=len(Q);need(all(Q[i][j]==Q[j][i] for i in range(d) for j in range(d)),'exact symmetric Q')
  P=features(data,k);G=[[sum((P[i][a]*Q[a][b]*P[j][b] for a in range(d) for b in range(d)),F(0)) for j in range(len(P))] for i in range(len(P))]
  for i,row in enumerate(G):
   for j,g in enumerate(row):
    for (a,b,c),v in H[k].items():add((a+i,b+j,c),-v*g)
    add((i+j+k,0,0),-g/(2*m));add((0,i+j+k,0),-g/(2*m))
 C={i:v for i,v in C.items() if v}
 need(all(C.get((j,i,k),F(0))==v for (i,j,k),v in C.items()),'exact u/v symmetry of h')
 return C,shape

def transform(deg,lo,hi):
 M=[]
 for i in range(deg+1):
  row=[]
  for j in range(deg+1):
   row.append(sum((F(comb(i,k),comb(deg,k))*comb(j,k)*lo**(j-k)*(hi-lo)**k for k in range(min(i,j)+1)),F(0)))
  M.append(row)
 den=lcm(*(x.denominator for row in M for x in row));return [[x.numerator*(den//x.denominator) for x in row] for row in M],den

def initial_coefficients(data):
 C,shape=polynomial(data);den=lcm(*(v.denominator for v in C.values()));N=prod(shape);A=[0]*N
 for (i,j,k),v in C.items():A[(i*shape[1]+j)*shape[2]+k]=v.numerator*(den//v.denominator)
 U=F(121,72);box=((F(0),U),(F(0),U),(-U,U))
 for axis in range(3):
  M,d=transform(shape[axis]-1,*box[axis]);stride=prod(shape[axis+1:]);group=shape[axis]*stride;out=[0]*N
  for base in range(0,N,group):
   for offset in range(stride):
    vals=[A[base+j*stride+offset] for j in range(shape[axis])]
    for i,row in enumerate(M):out[base+i*stride+offset]=sum(x*y for x,y in zip(row,vals))
  A=out;den*=d
 values=[S*a//den+BIAS for a in A]
 need(0<=min(values)<=max(values)<1<<62,'packed coefficient lanes nonnegative and below 2^62')
 packed=sum(v<<(W*i) for i,v in enumerate(values))
 return packed,shape,box,{'initial_minimum_lane':min(values),'initial_maximum_lane':max(values),'common_denominator_bits':den.bit_length(),'coefficient_count':N,'exact_h_symmetric':True}

# Each integer stores N coefficients as base-2^64 digits.  Every input digit
# is between 0 and 2^62, and subdivision takes only convex midpoint averages;
# therefore pairwise sums are below 2^63 and cannot carry into another digit.
# A global right shift by one would transfer each upper digit's odd bit into
# the next digit's bit 63.  The low-63-bit-per-digit mask removes precisely
# these transferred bits; the retained digits are floor((a+b)/2), hence safe
# lower Bernstein coefficients.  All arithmetic uses unbounded Python ints.
#
# To compare all digits A_i >= T_i, put bit 63 in every digit before subtracting
# T.  Since 0<=A_i,T_i<2^62, subtraction has no inter-digit borrow, and the
# sentinel bit survives exactly when A_i>=T_i.  Lane bounds are checked before
# every such comparison; no floating-point bound or machine-int overflow is
# used anywhere in the accepting path.
class Packed:
 def __init__(self,shape):
  self.shape=shape;N=prod(shape);self.ones=((1<<(W*N))-1)//((1<<W)-1);self.high=self.ones<<(W-1);self.axisones=[];self.full=[];self.valid=[]
  for axis,dim in enumerate(shape):
   stride=prod(shape[axis+1:]);group=dim*stride;repeat=sum(1<<(W*i) for i in range(0,N,group));inside=((1<<(W*stride))-1)//((1<<W)-1)
   ones=[(repeat*inside)<<(W*i*stride) for i in range(dim)];self.axisones.append(ones);self.full.append([o*((1<<W)-1) for o in ones]);self.valid.append([sum(ones[:dim-j])*((1<<(W-1))-1) for j in range(dim)])
 def split(self,B,axis):
  deg=self.shape[axis]-1;shift=W*prod(self.shape[axis+1:]);A=B;left=B&self.full[axis][0];right=B&self.full[axis][deg]
  for j in range(1,deg+1):
   A=(((A+(A>>shift))>>1)&self.valid[axis][j]);left|=(A&self.full[axis][0])<<(j*shift);right|=A&self.full[axis][deg-j]
  return left,right
 def dominates(self,B,box,c):
  (ul,uh),(vl,vh),(tl,th)=box;dlo=max(F(1,3),ul+vl-2*th);dhi=uh+vh-2*tl;a=(dlo+dhi)/2
  need(a>0,'positive log tangent point')
  base=BIAS+ceilq(S*(F(ul+vl-2*tl,a)-1-c))+log_upper_scaled(a)
  T=base*self.ones;lo=hi=base
  for axis,slope in enumerate((uh-ul,vh-vl,-2*(th-tl))):
   vals=[ceilq(F(S*slope*i,a*(self.shape[axis]-1))) for i in range(self.shape[axis])];lo+=min(vals);hi+=max(vals);T+=sum(v*o for v,o in zip(vals,self.axisones[axis]))
  need(0<=lo<=hi<1<<62,'tangent threshold lane bounds prevent inter-lane borrow/carry')
  return ((B+self.high-T)&self.high)==self.high

def log_series(a,terms=18):
 z=(a-1)/(a+1);need(0<=z<=F(1,3),'log range reduction to [1,2]');z2=z*z;p=z;s=F(0)
 for j in range(terms):s+=2*p/(2*j+1);p*=z2
 return s,s+2*p/((2*terms+1)*(1-z2))
LOG2=log_series(F(2))
@lru_cache(maxsize=100000)
def log_upper_scaled(a):
 k=0;y=a
 while y<1:y*=2;k-=1
 while y>=2:y/=2;k+=1
 low,upper=log_series(y);upper+=k*LOG2[1 if k>=0 else 0]
 return ceilq(S*upper)

def prune(box):
 (ul,uh),(vl,vh),(tl,th)=box
 minabs=tl if tl>0 else -th if th<0 else F(0)
 return ul>vh or minabs**2>uh*vh or uh+vh-2*tl<F(1,3)

def selftest():
 # Validate packed midpoint subdivision against direct scalar integer arithmetic.
 import random
 shape=(4,3,2);p=Packed(shape);rng=random.Random(20261108);A=[BIAS+rng.randrange(-100000,100000) for _ in range(prod(shape))];B=sum(a<<(W*i) for i,a in enumerate(A))
 for axis in range(3):
  L,R=p.split(B,axis);stride=prod(shape[axis+1:]);group=shape[axis]*stride;expectedL=[0]*len(A);expectedR=[0]*len(A)
  for base in range(0,len(A),group):
   for offset in range(stride):
    q=[A[base+j*stride+offset] for j in range(shape[axis])];left=[q[0]];right=[q[-1]]
    while len(q)>1:q=[(a+b)//2 for a,b in zip(q,q[1:])];left.append(q[0]);right.append(q[-1])
    for j,(l,r) in enumerate(zip(left,reversed(right))):expectedL[base+j*stride+offset]=l;expectedR[base+j*stride+offset]=r
  need(L==sum(a<<(W*i) for i,a in enumerate(expectedL)) and R==sum(a<<(W*i) for i,a in enumerate(expectedR)),'packed/scalar split identity')

def polynomial_crosscheck(data):
 # Independent evaluation of the prior encoded-kernel expression at exact
 # rational points is a regression check, not a substitute for the formulas
 # above or for the complete Bernstein replay.
 import importlib.util
 spec=importlib.util.spec_from_file_location('upper_prior_kernel_acceptor',ROOT/'candidate-kernel/verify_all.py');prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
 C,shape=polynomial(data);n=data['n'];B=hx(data['B']);b=(B-F(n*(n-1),2))/n
 points=[(F(1,4),F(2,3),-F(1,5)),(F(7,6),F(3,2),-F(3,4)),tuple(map(F,data['test_pair'])),(F(0),F(121,72),F(0))]
 for u,v,t in points:
  def p(a):return (F(n,2)*a+b-prior.encoded_kernel(a,a,a,data)/2)/(n-1)
  direct=p(u)+p(v)-prior.encoded_kernel(u,v,t,data)
  expanded=sum((q*u**i*v**j*t**k for (i,j,k),q in C.items()),F(0))
  need(expanded==direct,'independent exact h polynomial evaluation')
 return len(points)

def verify():
    """Replay the complete fixed tree; return deterministic exact proof facts."""
    selftest()
    c=F(21,10000)
    path=ROOT/'candidate-kernel/data/n12_exact_interpretation.json'
    raw=path.read_bytes();data=json.loads(raw)
    crosschecks=polynomial_crosscheck(data)
    need((data['n'],data['mass'],data['feature_order'],data['harmonics_through'])==(12,11,5,11),
         'fixed mass eleven n12 kernel')
    B,shape,box,initial=initial_coefficients(data)
    need(shape==(20,20,12),'fixed Bernstein tensor degrees')
    p=Packed(shape)
    actions_path=Path(__file__).with_name('cover.txt')
    tree_raw=actions_path.read_bytes();actions=tree_raw.decode('ascii').strip()
    stack=[(box,B,0)]
    splits=accepted=pruned=visited=0;maxdepth=0
    for action in actions:
        need(bool(stack),'tree action has a remaining box')
        box,B,depth=stack.pop();visited+=1;maxdepth=max(depth,maxdepth)
        if action=='P':
            need(prune(box),'exact domain discard');pruned+=1
        elif action=='A':
            need(p.dominates(B,box,c),'exact continuous Bernstein/log-tangent leaf');accepted+=1
        else:
            need(action in '012','fixed split axis')
            axis=int(action);left,right=p.split(B,axis);mid=sum(box[axis])/2
            bl=list(box);br=list(box)
            bl[axis]=(box[axis][0],mid);br[axis]=(mid,box[axis][1])
            stack.append((tuple(br),right,depth+1));stack.append((tuple(bl),left,depth+1));splits+=1
    need(not stack and visited==2*splits+1 and accepted+pruned==splits+1,
         'complete binary partition conservation')
    return {'passed':True,'scope':'exact continuous h-log(d)>=-21/10000 on the full stated necessary domain',
            'n':12,'mass':11,'U':'121/72','dmin':'1/3','c':str(c),
            'B':str(hx(data['B'])),'visited':visited,'split':splits,'accepted':accepted,
            'domain_discarded':pruned,'pending':0,'maximum_depth':maxdepth,
            'input_sha256':hashlib.sha256(raw).hexdigest(),
            'tree_sha256':hashlib.sha256(tree_raw).hexdigest(),
            'arithmetic':'stdlib Fraction and packed unbounded Python integers; 64-bit digit positions, no machine integers',
            'initial_coefficients':initial,'independent_exact_polynomial_evaluations':crosschecks,
            'sharp_candidate_bound_proved':False}


def main():
    print(json.dumps(verify(),indent=2))


if __name__=='__main__':main()
