"""Exploratory SDP with boundary grids and bounded Bernstein violation feedback.

This is numerical discovery only. It does not accept a global certificate.
Frozen input features are retained exactly; only Gram matrices and B vary.
"""
from pathlib import Path
import argparse,importlib.util,json,math,time
import numpy as np
import cvxpy as cp
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[2]
SCRIPT=Path(__file__).resolve().parent
HERE=ROOT/'work/jacobian-search'

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

prior=load(ROOT/'candidate-kernel/discover_sdp.py','cut_prior')
bern=load(SCRIPT/'bernstein_discovery.py','cut_bern')
U=121/72;DMIN=1/3;N=66

def frozen_kernel():
    data=json.loads((ROOT/'candidate-kernel/data/n12_exact_interpretation.json').read_bytes())
    x,nodes,_=prior.target(12);K=prior.Kernel(x,nodes,5,11)
    K.nodes=np.array([float.fromhex(v) for v in data['radial_nodes']])
    K.moments=np.array([float.fromhex(v) for v in data['radial_moment_centers']])
    K.transforms=[np.array([[float.fromhex(v) for v in row] for row in E]) for E in data['transforms']]
    K.dims=[5]+[E.shape[1] for E in K.transforms]
    K.logtarget=float.fromhex(data['B'])
    return K,data

def canonical(points):
    p=np.asarray(points,dtype=float).reshape(-1,3)
    p[:,:2]=np.sort(p[:,:2],axis=1)
    p=np.unique(np.round(p,13),axis=0)
    return p[(p[:,0]>=0)&(p[:,1]<=U+1e-12)&(p[:,0]*p[:,1]>=p[:,2]**2-1e-12)&(p[:,0]+p[:,1]-2*p[:,2]>=DMIN-1e-12)]

def initial_points():
    random=prior.sample_domain(2200,20261108,n=12)
    pts=[*random]
    # The former random/SLSQP search missed two distinct troughs on v=U.
    for u in np.linspace(0,U,61):
        for rho in np.linspace(-1,1,121):
            t=rho*math.sqrt(u*U)
            if u+U-2*t>=DMIN:pts.append([u,U,t])
    for u in np.linspace(0,U,31):
        for v in np.linspace(u,U,31):
            for t in (-math.sqrt(u*v),math.sqrt(u*v),(u+v-DMIN)/2):
                if abs(t)<=math.sqrt(u*v)+1e-13 and u+v-2*t>=DMIN-1e-13:pts.append([u,v,t])
    pts.extend([[.5264,U,.9033],[1.1945,U,-.7353],[1.2388,U,-1.1359]])
    for name in ('float-bernstein-c0.002.json','float-bernstein-c0.00065.json'):
        for entry in json.loads((SCRIPT/'seeds'/name).read_bytes()).get('unresolved',[]):
            pts.append([sum(pair)/2 for pair in entry['box']])
    return canonical(pts)

def export(data,Qs,eps):
    out=dict(data);out['B']=float(float.fromhex(data['B'])+eps).hex()
    out['Q_blocks']=[[[float(v).hex() for v in row] for row in ((Q+Q.T)/2)] for Q in Qs]
    return out

def check_bounds(data,K,Qs,eps,c=2e-5,maxnodes=200000,seconds=12):
    start=time.monotonic();C=bern.coefficients(data);box=((0.,U),(0.,U),(-U,U));B=bern.bernstein(C,box)
    stack=[(box,B,0)];visited=accepted=pruned=0;maxdepth=0;leaves=[]
    while stack and visited<maxnodes and time.monotonic()-start<seconds:
        box,B,depth=stack.pop();visited+=1;maxdepth=max(depth,maxdepth)
        (ul,uh),(vl,vh),(tl,th)=box
        amin=tl if tl>0 else -th if th<0 else 0
        if ul>vh+1e-12 or amin*amin>uh*vh+1e-12 or uh+vh-2*tl<DMIN-1e-12:
            pruned+=1;continue
        dlo=max(DMIN,ul+vl-2*th);dhi=uh+vh-2*tl;a=(dlo+dhi)/2
        R=B-bern.dcoeff(box,B.shape)/a
        lower=float(R.min())+1-math.log(a)
        if lower>=-c:accepted+=1;continue
        if depth>=18:
            # Save actual feasible points, not just an interval lower bound.
            for u in (ul,(ul+uh)/2,uh):
                for v in (vl,(vl+vh)/2,vh):
                    upper=min(math.sqrt(u*v),(u+v-DMIN)/2)
                    lower_t=-math.sqrt(u*v)
                    if upper>=lower_t:
                        t=min(upper,max(lower_t,(tl+th)/2))
                        leaves.append([u,v,t])
            continue
        axis=int(np.argmax([float(np.ptp(R,axis=j).max()) for j in range(3)]));mid=sum(box[axis])/2
        L,RB=bern.split(B,axis);bl=list(box);br=list(box)
        bl[axis]=(box[axis][0],mid);br[axis]=(mid,box[axis][1])
        stack.append((tuple(br),RB,depth+1));stack.append((tuple(bl),L,depth+1))
    random=prior.sample_domain(18000,20261109,n=12)
    samples=canonical([*leaves,*random])
    if len(samples):
        vals=K.residual(samples,Qs,eps);order=np.argsort(vals);chosen=[]
        for i in order:
            if len(chosen)>=55:break
            if vals[i]>-1e-8 and len(chosen)>=15:break
            if all(np.linalg.norm((samples[i]-samples[j])/U)>.045 for j in chosen):chosen.append(int(i))
        cuts=samples[chosen];worst=float(vals[order[0]]);where=samples[order[0]].tolist()
    else:cuts=np.empty((0,3));worst=None;where=None
    return {'visited':visited,'accepted':accepted,'pruned':pruned,'pending':len(stack),'maximum_depth':maxdepth,
            'deep_leaf_samples':len(samples),'worst_actual_sample_residual':worst,'worst_pair':where,
            'float_cover_complete':not stack and not leaves,'c':c,'seconds':time.monotonic()-start},cuts

def refine(K,Qs,eps,starts):
    # Include hard two-dimensional boundary minimizations explicitly.
    found=[]
    def f(p):return float(K.residual(p,Qs,eps)[0])
    constraints=[{'type':'ineq','fun':lambda p:p[0]*p[1]-p[2]**2},
                 {'type':'ineq','fun':lambda p:p[0]+p[1]-2*p[2]-DMIN}]
    for start in starts[:35]:
        r=minimize(f,start,method='SLSQP',bounds=[(0,U),(0,U),(-U,U)],constraints=constraints,
                   options={'ftol':1e-11,'maxiter':160})
        if r.success:found.append((f(r.x),r.x.tolist()))
    for u0,t0 in ((.5264,.9033),(1.1945,-.7353),(1.2388,-1.1359),(.8,0),(U,-U)):
        cons=[{'type':'ineq','fun':lambda p:p[0]*U-p[1]**2},
              {'type':'ineq','fun':lambda p:p[0]+U-2*p[1]-DMIN}]
        r=minimize(lambda p:f([p[0],U,p[1]]),[u0,t0],method='SLSQP',bounds=[(0,U),(-U,U)],constraints=cons,
                   options={'ftol':1e-12,'maxiter':200})
        if r.success:found.append((f([r.x[0],U,r.x[1]]),[float(r.x[0]),U,float(r.x[1])]))
    found.sort();return found

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--iterations',type=int,default=7);parser.add_argument('--start',type=int,default=0);args=parser.parse_args()
    HERE.mkdir(parents=True,exist_ok=True);K,data=frozen_kernel();points=initial_points();history=[];start=time.monotonic()
    if args.start:
        history=json.loads((HERE/'history.json').read_bytes())
        cuts=[]
        for row in history:
            cuts.extend(p for v,p in row['best_local_violations'])
            if row['float_bernstein']['worst_pair']:cuts.append(row['float_bernstein']['worst_pair'])
        for p in HERE.glob('candidate-*-c*.json'):
            for entry in json.loads(p.read_bytes()).get('unresolved_first',[]):cuts.append([sum(pair)/2 for pair in entry['box']])
        points=canonical([*points,*cuts])
    Qs=[cp.Variable((d,d),symmetric=True) for d in K.dims];eps=cp.Variable(nonneg=True)
    for iteration in range(args.start,args.start+args.iterations):
        base,matrices=K.coefficients(points)
        expression=base+eps/N+sum(A@cp.reshape(Q,(Q.size,),order='C') for A,Q in zip(matrices,Qs))
        constraints=[expression>=1e-9]+[Q-2e-7*np.eye(d)>>0 for Q,d in zip(Qs,K.dims)]
        problem=cp.Problem(cp.Minimize(eps+1e-9*sum(cp.trace(Q) for Q in Qs)),constraints)
        t=time.monotonic()
        problem.solve(solver='CLARABEL',max_iter=300,tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
        if eps.value is None:
            print('No SDP solution',problem.status,flush=True);break
        current=[(np.asarray(Q.value)+np.asarray(Q.value).T)/2 for Q in Qs];gap=float(eps.value)
        out=export(data,current,gap)
        (HERE/f'candidate-{iteration}.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
        cover,cuts=check_bounds(out,K,current,gap)
        numerical=refine(K,current,gap,[*cuts,*points[np.argsort(K.residual(points,current,gap))[:10]]])
        if numerical:
            cuts=canonical([*cuts,*[p for v,p in numerical[:25]]])
        row={'iteration':iteration,'SDP_status':problem.status,'gap':gap,'SDP_seconds':time.monotonic()-t-cover['seconds'],
             'points':len(points),'float_bernstein':cover,'best_local_violations':numerical[:8],
             'minimum_Q_eigenvalue_float':min(float(np.linalg.eigvalsh(Q).min()) for Q in current),
             'trace':sum(float(np.trace(Q)) for Q in current),'scope':'numerical discovery only'}
        history.append(row);(HERE/'history.json').write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(row),flush=True)
        if cover['float_cover_complete'] and (not numerical or numerical[0][0]>=-2e-5):break
        points=canonical([*points,*cuts]);np.save(HERE/'current_points.npy',points)
    print('Saved bounded discovery history; elapsed',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
