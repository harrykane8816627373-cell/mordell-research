"""Outward-rounded Bernstein covering for the exact auxiliary polynomial.

Every discarded box is either outside the domain, inside a rigorously convex
contact neighborhood, or has a positive Bernstein lower bound after subtracting
a tangent upper bound to log. Unfinished frontiers are explicitly retained.
"""
from pathlib import Path
import sys,json,time,math,hashlib,itertools,collections
import numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT
if not __debug__:raise RuntimeError('Run without -O.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
from flint import arb,ctx,fmpq
from arb_bernstein import coefficients as arb_coefficients, strictly_positive as arb_positive
from arb_bernstein import coefficients_and_cubic_residual
from cubic_log_bernstein import CubicLogUpper
from physical_boundary_bernstein import correction_upper,exact_q_coefficients
ctx.prec=160;np.seterr(over='raise',invalid='raise')
source=OUT/'global_pair_bound_exact_contacts.json';data=source.read_bytes();sha=hashlib.sha256(data).hexdigest();raw=json.loads(data)
a=arb(2).sqrt();cc=[a/2,arb(0),-a/2,arb(-1)];logs=[(2-a).log(),arb(2).log()]
def fld(coefs):
    val=arb(0)
    for q in coefs:val=val*a+arb(fmpq(q))
    return val
def scalar(value):return fld(value[0])+sum((fld(c)*l for c,l in zip(value[1:],logs)),arb(0))
encoded={tuple(e['powers']):e['value']for e in raw['polynomial_terms']}
assert all(encoded[(j,i,k)]==value for (i,j,k),value in encoded.items())
P={e:scalar(v)for e,v in encoded.items()}
local_name='global_pair_bound_local_preconditioned_certificate.json'
local=json.loads((OUT/local_name).read_text())
definition=json.loads((OUT/'global_pair_bound_exact_definition_certificate.json').read_text())
assert definition['coefficient_sha256']==sha and definition['kernel_definition_exact']and definition['contact_values_and_gradients_exact']
assert local['coefficient_sha256']==sha and local['local_positivity_proved']and local['all_matrices_positive_definite']
assert [r['contact']for r in local['records']]==list(range(5))
contacts=[[arb(0),arb(1),arb(0)]]+[[arb(1),arb(1),c]for c in cc]
neighborhoods=[(p,rec['radius'])for p,rec in zip(contacts,local['records'])]
neighborhoods.append(([arb(1),arb(0),arb(0)],neighborhoods[0][1]))
U=1860281/1048576;UA=arb(U);MIN=arb(1)/4
radial_rho=(225-arb(18465).sqrt())/90
radial_L=arb(14372)/405-arb(20)/81*arb(18465).sqrt()
beta=(radial_rho+(radial_rho*radial_rho+4*radial_L.sqrt()).sqrt())/2
assert UA>beta and radial_L>0 and radial_rho>0
DEGREES=tuple(max(e[j]for e in P)for j in range(3));SHAPE=tuple(n+1 for n in DEGREES)
cubic=CubicLogUpper(DEGREES)
DOWN=-np.inf;UP=np.inf
def lower(x):return np.nextafter(float(x.lower()),DOWN)
def upper(x):return np.nextafter(float(x.upper()),UP)

P_numeric=[(e,float(c.mid()))for e,c in P.items()]
Pt_numeric=[((i,j,k-1),c*k)for (i,j,k),c in P_numeric if k]
def normal_multiplier(box):
    u,v,t=[(box[2*j]+box[2*j+1])/2 for j in range(3)]
    d=u+v-2*t
    if d<=0 or abs(t)<1e-8:return 0.
    derivative=sum(c*u**i*v**j*t**k for (i,j,k),c in Pt_numeric)+2/d
    lam=derivative/(-2*t)
    return lam if math.isfinite(lam)and lam>0 else 0.
def remaining_radial_lower(s):
    candidates=[]
    for k in range(8):
        r=8-s-k*radial_rho
        candidates.append(arb((r*r).lower())if r>0 or r<0 else arb(0))
    return radial_rho*radial_rho*min(candidates)/28

def strict_failure_at_center(box):
    xyz=[(box[2*j]+box[2*j+1])/2 for j in range(3)];u,v,t=xyz;d=u+v-2*t
    if d<.25 or u*v<t*t:return None
    hint=sum(c*u**e[0]*v**e[1]*t**e[2]for e,c in P_numeric)-math.log(d)
    if hint>=-1e-9:return None
    x,y,z=map(arb,xyz);dist=x+y-2*z
    if not (dist>=MIN and x*y-z*z>=0):return None
    xr=x-radial_rho;yr=y-radial_rho
    remainder_upper=radial_rho*radial_rho*min(arb(((8-x-y-k*radial_rho)*(8-x-y-k*radial_rho)).upper())for k in range(8))/28
    if not x*x*xr*xr+y*y*yr*yr+remainder_upper<radial_L:return None
    aa=x*x+y*y;ss=x+y;V=arb(256)/27;wstar=aa+arb(5)/4*(8-ss)
    if wstar>=V:w=V
    elif wstar<=V:w=wstar
    else:return None
    if not 2*(x*x*x+y*y*y+(w-aa)*(w-aa)/(8-ss))-5*w+arb(200)/9<0:return None
    # Evaluate at exact binary coordinates with directed ball arithmetic.
    value=sum((c*x**e[0]*y**e[1]*z**e[2]for e,c in P.items()),arb(0))-dist.log()
    if not value<0:return None
    return {'point_binary_floats':xyz,'strict_residual_ball':str(value),'meaning':'Counterexample to this auxiliary majorant on its necessary pair domain, not a nine-point counterexample.'}

def initial_coefficients():
    balls=arb_coefficients(P,DEGREES,(0.,U,0.,U,-U,U))
    lo=np.empty(SHAPE);hi=np.empty(SHAPE)
    for ix in np.ndindex(SHAPE):
        val=balls[ix]
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
    counts={'visited':0,'outside':0,'moment_outside':0,'quartic_outside':0,'local':0,'positive':0,'cubic_positive':0,'arb_recomputed':0,'arb_positive':0,'split':0}
counts.setdefault('cubic_positive',0)
counts.setdefault('arb_refreshed',0)
counts.setdefault('discrete_radial_outside',0)
counts.setdefault('physical_correction_positive',0)
counts.setdefault('quadratic_correction_positive',0)
counts.setdefault('arb_correction_positive',0)
print('BERNSTEIN START',SHAPE,'frontier',len(queue),flush=True)
start=time.monotonic()
stop_reason='time_limit';failure=None
while queue and time.monotonic()-start<seconds:
    box,bl,bh=queue.pop();lu,hu,lv,hv,lt,ht=box;counts['visited']+=1
    dupper=arb(hu)+arb(hv)-2*arb(lt);tabs=0 if lt<=0<=ht else min(abs(lt),abs(ht))
    if lu>hv or dupper<MIN or arb(tabs)**2>arb(hu)*arb(hv):counts['outside']+=1;continue
    # Necessary joint-radius condition for EVERY local maximum, including degenerate ones.
    ub=arb(lu).union(arb(hu));vb=arb(lv).union(arb(hv));ss=ub+vb;aa=ub**2+vb**2;bb=ub**3+vb**3
    ur=ub-radial_rho;vr=vb-radial_rho
    quartic_pair=ub*ub*ur*ur+vb*vb*vr*vr
    if quartic_pair>radial_L:
        counts['quartic_outside']+=1;continue
    # This optional necessary test is useful only near the quartic boundary.
    # Skipping it elsewhere keeps a larger domain, never removes a valid box.
    if quartic_pair>radial_L-arb(1)/10 and quartic_pair+remaining_radial_lower(ss)>radial_L:
        counts['discrete_radial_outside']+=1;continue
    V=arb(256)/27;wstar=aa+arb(5)/4*(8-ss)
    if wstar>=V and 2*(bb+(V-aa)**2/(8-ss))-5*V+arb(200)/9>0:
        counts['moment_outside']+=1;continue
    bounds=[(lu,hu),(lv,hv),(lt,ht)]
    if any(all(arb(lo)>=c-eps and arb(hi)<=c+eps for (lo,hi),c in zip(bounds,p))for p,eps in neighborhoods):
        counts['local']+=1;continue
    tangent=tangent_upper_coefficients(box)
    if np.min(np.nextafter(bl-tangent,DOWN))>0:counts['positive']+=1;continue
    higher=cubic.upper(box)
    if np.min(np.nextafter(bl-higher,DOWN))>0:counts['cubic_positive']+=1;continue
    tb=arb(lt).union(arb(ht))
    qbox=ub*vb-tb*tb
    if not qbox>0:
        lam=normal_multiplier(box)
        if lam:
            correction=correction_upper(box,DEGREES,lam)
            corrected=np.nextafter(np.nextafter(bl-higher,DOWN)-correction,DOWN)
            if np.min(corrected)>0:counts['physical_correction_positive']+=1;continue
            proved=False
            for kappa in [1.,4.]:
                if not arb(kappa)*qbox<=arb(lam):continue
                correction=correction_upper(box,DEGREES,lam,kappa)
                corrected=np.nextafter(np.nextafter(bl-higher,DOWN)-correction,DOWN)
                if np.min(corrected)>0:proved=True;break
            if proved:counts['quadratic_correction_positive']+=1;continue
    failure=strict_failure_at_center(box)
    if failure:
        queue.append((box,bl,bh));stop_reason='strict_auxiliary_counterexample';break
    widths=[hi-lo for lo,hi in bounds];axis=int(np.argmax(widths));lo,hi=bounds[axis];mid=(lo+hi)/2
    if max(widths)<.001 and (np.max(bh-bl)>1e-12 or max(widths)<1e-6):
        counts['arb_recomputed']+=1
        exact_coefficients,exact_residual=coefficients_and_cubic_residual(P,DEGREES,box)
        if all(x>0 for x in exact_residual.flat):counts['arb_positive']+=1;continue
        if not qbox>0:
            lam=normal_multiplier(box)
            if lam:
                qc,q2c=exact_q_coefficients(box,DEGREES);proved=False
                for kappa in [0.,1.,4.]:
                    if kappa and not arb(kappa)*qbox<=arb(lam):continue
                    if all(r-arb(lam)*q+arb(kappa)*q2>0 for r,q,q2 in zip(exact_residual.flat,qc.flat,q2c.flat)):
                        proved=True;break
                if proved:counts['arb_correction_positive']+=1;continue
        # Refresh both endpoints even when subdivision is still needed. This
        # prevents inherited rounding widths from dominating every descendant.
        refreshed_lo=np.array([lower(x)for x in exact_coefficients.flat]).reshape(SHAPE)
        refreshed_hi=np.array([upper(x)for x in exact_coefficients.flat]).reshape(SHAPE)
        if np.any(refreshed_lo>bh)or np.any(refreshed_hi<bl):
            queue.append((box,bl,bh));stop_reason='inconsistent_coefficient_enclosures';break
        bl,bh=refreshed_lo,refreshed_hi
        counts['arb_refreshed']+=1
    if mid==lo or mid==hi:
        queue.append((box,bl,bh));stop_reason='unresolved_at_float_endpoint_resolution';break
    if not arb(mid)==(arb(lo)+arb(hi))/2:
        queue.append((box,bl,bh));stop_reason='unresolved_at_inexact_binary_midpoint';break
    l1,h1,l2,h2=split_coefficients(bl,bh,axis)
    b1=list(box);b2=list(box);b1[2*axis+1]=mid;b2[2*axis]=mid
    # Exact midpoints are necessary for the half-interval subdivision above.
    assert arb(mid)==(arb(lo)+arb(hi))/2,(counts,box,axis,mid)
    queue.append((tuple(b1),l1,h1));queue.append((tuple(b2),l2,h2));counts['split']+=1
    if counts['visited']%2000==0:print('BERNSTEIN',counts,'pending',len(queue),'seconds',round(time.monotonic()-start,1),flush=True)
state={'coefficient_sha256':sha,'U':U,'minimum_distance_squared':'1/4','additional_domain':['TWO_RADIUS_DOMAIN.txt','QUARTIC_RADIAL_MOMENT_PROOF.txt','DISCRETE_RADIAL_REMAINDER_PROOF.txt'],'degrees':DEGREES,'counts':counts,
       'frontier':[{'box':box,'lower':l.tolist(),'upper':h.tolist()}for box,l,h in queue],
       'global_majorant_proved':not queue,'local_certificate':local_name,'seconds_this_run':time.monotonic()-start,
       'stop_reason':'complete'if not queue else stop_reason,'strict_auxiliary_counterexample':failure}
state_path.write_text(json.dumps(state))
print('BERNSTEIN STOP',counts,'UNRESOLVED',len(queue),'GLOBAL PROVED',not queue,flush=True)
