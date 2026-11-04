"""Outward-rounded Bernstein covering for the exact auxiliary polynomial.

Every discarded box is either outside the domain, inside a rigorously convex
contact neighborhood, or has a positive Bernstein lower bound after subtracting
a tangent upper bound to log. Unfinished frontiers are explicitly retained.
"""
from pathlib import Path
import sys,json,time,math,hashlib,itertools,collections
import numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data'
sys.path.insert(0,str(ROOT.parent.parent/'friedman_full_audit_20261023/deps'))
from flint import arb,ctx,fmpq
ctx.prec=160;np.seterr(over='raise',invalid='raise')
source=OUT/'global_pair_bound_exact_contacts.json';data=source.read_bytes();sha=hashlib.sha256(data).hexdigest();raw=json.loads(data)
a=(2*arb.pi()/7).cos();cc=[a,2*a*a-1,-2*a*a-a+arb(1)/2];logs=[(2-2*c).log()for c in cc]
def fld(coefs):
    val=arb(0)
    for q in coefs:val=val*a+arb(fmpq(q))
    return val
def scalar(value):return fld(value[0])+sum((fld(c)*l for c,l in zip(value[1:],logs)),arb(0))
encoded={tuple(e['powers']):e['value']for e in raw['polynomial_terms']}
assert all(encoded[(j,i,k)]==value for (i,j,k),value in encoded.items())
P={e:scalar(v)for e,v in encoded.items()}
local=json.loads((OUT/'global_pair_bound_local_certificate.json').read_text())
assert local['coefficient_sha256']==sha and local['local_positivity_proved']
contacts=[[arb(0),arb(1),arb(0)]]+[[arb(1),arb(1),c]for c in cc]
neighborhoods=[(p,rec['radius'])for p,rec in zip(contacts,local['records'])]
neighborhoods.append(([arb(1),arb(0),arb(0)],neighborhoods[0][1]))
U=math.ceil((18411/10000)*2**20)/2**20;UA=arb(U);MIN=arb(2)/7
smax=UA*8/7
assert 882*smax**3-2352*smax**2+1755*smax-1496>0 and smax>arb('21/10')
DEGREES=tuple(max(e[j]for e in P)for j in range(3));SHAPE=tuple(n+1 for n in DEGREES)
DOWN=-np.inf;UP=np.inf
def lower(x):return np.nextafter(float(x.lower()),DOWN)
def upper(x):return np.nextafter(float(x.upper()),UP)

def initial_coefficients():
    # Substitute u=U*x, v=U*y, t=-U+2U*z into the exact polynomial.
    power={}
    for (i,j,k),c in P.items():
        for ell in range(k+1):
            e=(i,j,ell);val=c*UA**(i+j)*math.comb(k,ell)*(-UA)**(k-ell)*(2*UA)**ell
            power[e]=power.get(e,arb(0))+val
    lo=np.empty(SHAPE);hi=np.empty(SHAPE)
    for ix in np.ndindex(SHAPE):
        val=arb(0)
        for e,c in power.items():
            if all(j<=i for i,j in zip(ix,e)):
                q=fmpq(1)
                for i,j,n in zip(ix,e,DEGREES):q*=fmpq(math.comb(i,j),math.comb(n,j))
                val+=c*arb(q)
        lo[ix]=lower(val);hi[ix]=upper(val)
    return lo,hi

def split_coefficients(lo,hi,axis):
    ll=np.moveaxis(lo,axis,0).copy();hh=np.moveaxis(hi,axis,0).copy();n=ll.shape[0]-1
    leftl=np.empty_like(ll);lefth=np.empty_like(hh);rightl=np.empty_like(ll);righth=np.empty_like(hh)
    leftl[0]=ll[0];lefth[0]=hh[0];rightl[n]=ll[n];righth[n]=hh[n]
    for r in range(1,n+1):
        # Each elementary addition and multiplication rounds outward.
        ll=np.nextafter(np.nextafter(ll[:-1]+ll[1:],DOWN)*.5,DOWN)
        hh=np.nextafter(np.nextafter(hh[:-1]+hh[1:],UP)*.5,UP)
        leftl[r]=ll[0];lefth[r]=hh[0];rightl[n-r]=ll[-1];righth[n-r]=hh[-1]
    return tuple(np.moveaxis(q,0,axis)for q in [leftl,lefth,rightl,righth])

def tangent_upper_coefficients(box):
    lo_u,hi_u,lo_v,hi_v,lo_t,hi_t=box
    center_d=(lo_u+hi_u+lo_v+hi_v)/2-lo_t-hi_t
    d0=max(center_d,.25) # Any strictly positive exact binary number is valid.
    d=arb(d0);inv=1/d;cl=d.log()-1
    invlo,invhi=lower(inv),upper(inv);ch=upper(cl)
    grids=[]
    for axis,(aa,bb)in enumerate([(lo_u,hi_u),(lo_v,hi_v),(lo_t,hi_t)]):
        n=DEGREES[axis];vals=[arb(aa)+(arb(bb)-arb(aa))*arb(fmpq(i,n))for i in range(n+1)]
        grids.append((np.array([lower(x)for x in vals]),np.array([upper(x)for x in vals])))
    ul,uh=grids[0];vl,vh=grids[1];tl,th=grids[2]
    # D=u+v-2t; only its upper endpoint is needed since 1/d0 is positive.
    Dh=np.nextafter(np.nextafter(uh[:,None,None]+vh[None,:,None],UP)-np.nextafter(2*tl[None,None,:],DOWN),UP)
    product=np.maximum(np.nextafter(Dh*invlo,UP),np.nextafter(Dh*invhi,UP))
    return np.nextafter(ch+product,UP)

state_path=OUT/'global_pair_bound_bernstein_frontier.json'
seconds=180
for arg in sys.argv[1:]:
    if arg.isdigit():seconds=int(arg)
resume='--resume'in sys.argv
if resume and state_path.exists():
    state=json.loads(state_path.read_text());assert state['coefficient_sha256']==sha and state['U']==U
    queue=collections.deque([(tuple(x['box']),np.array(x['lower']),np.array(x['upper']))for x in state['frontier']]);counts=state['counts']
else:
    l,h=initial_coefficients();queue=collections.deque([((0.,U,0.,U,-U,U),l,h)])
    counts={'visited':0,'outside':0,'local':0,'positive':0,'split':0}
print('BERNSTEIN START',SHAPE,'frontier',len(queue),flush=True)
start=time.monotonic()
while queue and time.monotonic()-start<seconds:
    box,bl,bh=queue.pop();lu,hu,lv,hv,lt,ht=box;counts['visited']+=1
    dupper=arb(hu)+arb(hv)-2*arb(lt);tabs=0 if lt<=0<=ht else min(abs(lt),abs(ht))
    if lu>hv or dupper<MIN or arb(tabs)**2>arb(hu)*arb(hv):counts['outside']+=1;continue
    bounds=[(lu,hu),(lv,hv),(lt,ht)]
    if any(all(arb(lo)>=c-eps and arb(hi)<=c+eps for (lo,hi),c in zip(bounds,p))for p,eps in neighborhoods):
        counts['local']+=1;continue
    tangent=tangent_upper_coefficients(box)
    if np.min(np.nextafter(bl-tangent,DOWN))>0:counts['positive']+=1;continue
    widths=[hi-lo for lo,hi in bounds];axis=int(np.argmax(widths));lo,hi=bounds[axis];mid=(lo+hi)/2
    if mid==lo or mid==hi:raise RuntimeError('Unresolved box at floating endpoint resolution')
    l1,h1,l2,h2=split_coefficients(bl,bh,axis)
    b1=list(box);b2=list(box);b1[2*axis+1]=mid;b2[2*axis]=mid
    # Exact midpoints are necessary for the half-interval subdivision above.
    assert arb(mid)==(arb(lo)+arb(hi))/2,(counts,box,axis,mid)
    queue.append((tuple(b1),l1,h1));queue.append((tuple(b2),l2,h2));counts['split']+=1
    if counts['visited']%2000==0:print('BERNSTEIN',counts,'pending',len(queue),'seconds',round(time.monotonic()-start,1),flush=True)
state={'coefficient_sha256':sha,'U':U,'minimum_distance_squared':'2/7','degrees':DEGREES,'counts':counts,
       'frontier':[{'box':box,'lower':l.tolist(),'upper':h.tolist()}for box,l,h in queue],
       'global_majorant_proved':not queue,'local_certificate':'global_pair_bound_local_certificate.json','seconds_this_run':time.monotonic()-start}
state_path.write_text(json.dumps(state))
print('BERNSTEIN STOP',counts,'UNRESOLVED',len(queue),'GLOBAL PROVED',not queue,flush=True)
