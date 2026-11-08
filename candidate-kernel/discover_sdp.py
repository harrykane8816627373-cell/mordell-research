"""Numerical candidate-adapted PSD kernel search; NOT a proof/acceptor.

Needs isolated numpy/scipy/cvxpy. Outputs are explicitly sampled numerical
evidence. Exact summation/nullspace identities are checked separately.
"""
from pathlib import Path
import argparse
import json
import math
import sys

import cvxpy as cp
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCRATCH = ROOT.parent/'kernel-candidate'


def target(n=10):
    raw = json.loads((ROOT/'validated-maxima/certificates.json').read_text())
    row = next(r for r in raw['rows'] if r['n'] == n)
    flat = np.array([float(q/(2**raw['center_exponent'])) for q in row['center_numerators']])
    flat = np.insert(flat, 1, 0)
    z = flat[::2]+1j*flat[1::2]
    rotation = z[8].conjugate()/abs(z[8]) if n==10 else 1
    x = z*rotation*math.sqrt((n-1)/n)
    # Discovery ONLY: enforce the D2 ansatz visible in the validated seed.
    if n==10:
        a = abs(x[8].real)
        b, c = abs(x[1].real), abs(x[1].imag)
        d, e = abs(x[0].real), abs(x[0].imag)
        x = np.array([a,-a]+[sx*b+1j*sy*c for sx in (-1,1) for sy in (-1,1)]
                     +[sx*d+1j*sy*e for sx in (-1,1) for sy in (-1,1)])
    x *= math.sqrt(n-1)/np.linalg.norm(x)
    us = abs(x)**2
    nodes=[]
    for u in sorted(us):
        if not nodes or u-nodes[-1]>1e-8:nodes.append(u)
    nodes=np.array(nodes)
    contacts = []
    for i in range(n):
        for j in range(i):
            u,v = sorted((us[i], us[j]))
            t = (x[i]*x[j].conjugate()).real
            if not any(np.linalg.norm(np.array((u,v,t))-p)<1e-9 for p in contacts):
                contacts.append(np.array([u,v,t]))
    return x, nodes, np.array(contacts)


class Kernel:
    def __init__(self, x, nodes, s=4, modes=9):
        self.s, self.modes = s, modes
        self.nodes = nodes
        self.n=len(x);self.mass=float(self.n-1)
        self.logtarget = sum(math.log(abs(a-b)**2) for i,a in enumerate(x) for b in x[:i])
        self.moments = np.array([np.mean(abs(x)**(2*j)) for j in range(1,s)])
        self.transforms = []
        self.raw_moments = []
        for k in range(1,modes+1):
            M = np.array([np.sum(x**k*(abs(x)**2-1)**j) for j in range(s)])
            self.raw_moments.append(M)
            if np.linalg.norm(M)<1e-8:
                E = np.eye(s)
            else:
                _,sing,vh = np.linalg.svd(np.vstack([M.real,M.imag]),full_matrices=True)
                rank=sum(q>1e-8 for q in sing)
                E = vh[rank:,:].T
            self.transforms.append(E)
        self.dims = [s]+[E.shape[1] for E in self.transforms]

    def feature(self, u, k):
        u = np.asarray(u)
        if k == 0:
            q = np.prod(u[...,None]-self.nodes,axis=-1)
            v = np.stack([u**j-self.moments[j-1] for j in range(1,self.s)]+[10*q*q],axis=-1)
            qp = sum(np.prod(u[...,None]-np.delete(self.nodes,i),axis=-1) for i in range(len(self.nodes)))
            derivative = np.stack([j*u**(j-1) for j in range(1,self.s)]+[20*q*qp],axis=-1)
            return v, derivative
        raw = np.stack([(u-1)**j for j in range(self.s)],axis=-1)
        der = np.stack([np.zeros_like(u)]+[j*(u-1)**(j-1) for j in range(1,self.s)],axis=-1)
        return raw@self.transforms[k-1], der@self.transforms[k-1]

    def coefficients(self, pts, gradients=False):
        pts = np.atleast_2d(pts)
        u,v,t = pts.T
        zero,one = np.zeros_like(u),np.ones_like(u)
        H=[one,t];Hu=[zero,zero];Hv=[zero,zero];Ht=[zero,one]
        for k in range(2,self.modes+1):
            H.append(2*t*H[-1]-u*v*H[-2])
            Hu.append(2*t*Hu[-1]-v*H[-3]-u*v*Hu[-2])
            Hv.append(2*t*Hv[-1]-u*H[-3]-u*v*Hv[-2])
            Ht.append(2*H[-2]+2*t*Ht[-1]-u*v*Ht[-2])
        matrices = []; derivs=[]
        for k in range(self.modes+1):
            fu,du = self.feature(u,k);fv,dv = self.feature(v,k)
            outer = np.einsum('ni,nj->nij',fu,fv)
            cross = (outer+outer.transpose(0,2,1))/2
            uu=np.einsum('ni,nj->nij',fu,fu);vv=np.einsum('ni,nj->nij',fv,fv)
            denom=2*(self.n-1)
            base = -H[k][:,None,None]*cross-(u**k)[:,None,None]*uu/denom-(v**k)[:,None,None]*vv/denom
            matrices.append(base.reshape(len(pts),-1))
            if gradients:
                ou=np.einsum('ni,nj->nij',du,fv);ov=np.einsum('ni,nj->nij',fu,dv)
                cu=(ou+ou.transpose(0,2,1))/2;cv=(ov+ov.transpose(0,2,1))/2
                uuder=np.einsum('ni,nj->nij',du,fu)+np.einsum('ni,nj->nij',fu,du)
                vvder=np.einsum('ni,nj->nij',dv,fv)+np.einsum('ni,nj->nij',fv,dv)
                du_diag=u**k
                dv_diag=v**k
                dbaseu=-Hu[k][:,None,None]*cross-H[k][:,None,None]*cu-du_diag[:,None,None]*uuder/denom
                dbasev=-Hv[k][:,None,None]*cross-H[k][:,None,None]*cv-dv_diag[:,None,None]*vvder/denom
                if k:
                    dbaseu-=(k*u**(k-1))[:,None,None]*uu/denom
                    dbasev-=(k*v**(k-1))[:,None,None]*vv/denom
                dbaset=-Ht[k][:,None,None]*cross
                derivs.append(np.stack([dbaseu,dbasev,dbaset],axis=1).reshape(len(pts),3,-1))
        A=self.n/2;m=self.n-1
        base=A*(u+v)/m+2*(self.logtarget-A*self.mass)/(self.n*m)
        residualbase=base-np.log(u+v-2*t)
        if gradients:
            d=u+v-2*t
            dbase=np.stack([A/m-1/d,A/m-1/d,2/d],axis=-1)
            return residualbase, matrices, dbase, derivs
        return residualbase,matrices

    def residual(self,pts,Qs,eps):
        base,matrices=self.coefficients(pts)
        return base+eps/(self.n*(self.n-1)/2)+sum(A@Q.reshape(-1) for A,Q in zip(matrices,Qs))


def schwarz_pair_margin(points):
    points=np.atleast_2d(points)
    u,v,t=points.T
    q=9-9*(u+v)/8-t/4
    delta=5*q-28
    d=u+v-2*t
    b=(81*u+v+18*t)*(u+81*v+18*t)/4096
    with np.errstate(divide='ignore',invalid='ignore'):
        raw=2*np.log(delta+2)+np.log(d)+8*np.log(b)+28*np.log(q)-2*np.log(2-delta)-math.log((97/80)*9**9)-21*math.log(7)
    raw=np.where(delta>=2,1,np.minimum(raw,1))
    return np.where(np.isfinite(raw),raw,-1e6)


def domain_filter(points,strong=False):
    if not strong:return np.ones(len(points),dtype=bool)
    u,v,t=np.asarray(points).T
    return ((9*(u+v)+2*t<=664/25+1e-12)&(schwarz_pair_margin(points)>=-1e-12))


def sample_domain(count,seed,strong=False,n=10):
    rng=np.random.default_rng(seed);U=15403/10000 if strong else 2*(n-1)**2/n**2
    dmin=4/n
    pts=[]
    while sum(len(p) for p in pts)<count:
        u=rng.uniform(0,U,count);v=rng.uniform(0,U,count)
        t=rng.uniform(-1,1,count)*np.sqrt(u*v)
        mask=u+v-2*t>=dmin
        p=np.stack([u[mask],v[mask],t[mask]],axis=-1)
        pts.append(p[domain_filter(p,strong)])
    pts=np.concatenate(pts)[:count]
    # Necessary-domain boundary and origin pairs included explicitly.
    for u in np.linspace(0,U,25):
        for v in np.linspace(0,U,25):
            for t in (-math.sqrt(u*v),math.sqrt(u*v),(u+v-dmin)/2):
                if abs(t)<=math.sqrt(u*v)+1e-13 and u+v-2*t>=dmin-1e-13:
                    pts=np.vstack((pts,[u,v,t]))
    return pts[domain_filter(pts,strong)]


def attempt(s=4,modes=9,count=1800,iterations=3,contact_gradients=True,strong=False,n=10):
    if strong and n!=10:raise ValueError('strong Schwarz domain currently n10 only')
    x,nodes,contacts=target(n);kernel=Kernel(x,nodes,s,modes)
    points=np.vstack((sample_domain(count,1238,strong,n),contacts))
    U=15403/10000 if strong else 2*(n-1)**2/n**2
    dmin=4/n;N=n*(n-1)/2
    Qs=[cp.Variable((r,r),PSD=True) for r in kernel.dims]
    eps=cp.Variable(nonneg=True)
    history=[]
    matrices=None
    for iteration in range(iterations):
        base,matrices=kernel.coefficients(points)
        expression=base+eps/N+sum(A@cp.reshape(Q,(Q.size,),order='C') for A,Q in zip(matrices,Qs))
        constraints=[expression>=-1e-10]
        if contact_gradients:
            _,_,dbase,Ds=kernel.coefficients(contacts,True)
            grad=dbase+sum(cp.reshape(D.reshape(len(contacts)*3,-1)@cp.reshape(Q,(Q.size,),order='C'),(len(contacts),3),order='C') for D,Q in zip(Ds,Qs))
            for i,(u,v,t) in enumerate(contacts):
                if u*v-t*t>1e-8:
                    constraints.extend([grad[i,:]<=1e-6,grad[i,:]>=-1e-6])
                else:
                    constraints.extend([cp.abs(grad[i,0]+t/(2*u)*grad[i,2])<=1e-6,
                                        cp.abs(grad[i,1]+t/(2*v)*grad[i,2])<=1e-6,
                                        grad[i,2]>=-1e-6])
        problem=cp.Problem(cp.Minimize(eps+1e-8*sum(cp.trace(Q) for Q in Qs)),constraints)
        try:
            problem.solve(solver='CLARABEL',max_iter=250,tol_gap_abs=1e-8,tol_feas=1e-8,tol_gap_rel=1e-8)
        except Exception as error:
            history.append({'iteration':iteration,'status':'solver_error','error':str(error)})
            break
        status=problem.status
        if eps.value is None:
            history.append({'iteration':iteration,'status':status})
            break
        current=[np.array(Q.value) for Q in Qs];gap=float(eps.value)
        check=sample_domain(10000,2501+iteration,strong,n)
        values=kernel.residual(check,current,gap)
        indices=np.argsort(values)[:16]
        violations=[(float(values[i]),check[i].tolist()) for i in indices]
        for start in contacts.tolist()+[check[i].tolist() for i in indices[:8]]:
            domain_constraints=[{'type':'ineq','fun':lambda p:p[0]*p[1]-p[2]**2},
                                {'type':'ineq','fun':lambda p:p[0]+p[1]-2*p[2]-dmin}]
            if strong:
                domain_constraints.extend([{'type':'ineq','fun':lambda p:664/25-9*(p[0]+p[1])-2*p[2]},
                                           {'type':'ineq','fun':lambda p:float(schwarz_pair_margin(p)[0])}])
            result=minimize(lambda p: float(kernel.residual(p,current,gap)[0]),start,
                            method='SLSQP',bounds=[(0,U),(0,U),(-U,U)],
                            constraints=domain_constraints,
                            options={'ftol':1e-11,'maxiter':250})
            if result.success:
                violations.append((float(result.fun),result.x.tolist()))
        violations.sort()
        worst=violations[0]
        history.append({'iteration':iteration,'status':status,'candidate_gap':gap,
                        'worst_numerical_domain_residual':worst[0],'worst_pair':worst[1],
                        'minimum_Q_eigenvalue':min(float(np.linalg.eigvalsh(Q).min()) for Q in current),
                        'trace_sum':sum(float(np.trace(Q)) for Q in current),
                        'samples':len(points)})
        print(history[-1],flush=True)
        if worst[0]>=-1e-7:
            break
        points=np.vstack((points,[p for value,p in violations[:16]]))
    out={'date':'2026-11-08','status':'numerical_discovery_only_not_a_global_certificate',
         'mass':n-1,'n':n,'candidate_log_delta':kernel.logtarget,
         'radial_nodes':nodes.tolist(),'feature_order':s,'harmonics_through':modes,
         'contact_classes':contacts.tolist(),'contact_class_count':len(contacts),
         'contact_gradients_imposed':contact_gradients,'history':history,
         'strong_schwarz_domain':strong,
         'scope':'floating SDP and numerical violation searches only; no continuous-domain coverage proved'}
    if eps.value is not None:
        out['Q_blocks']=[Q.value.tolist() for Q in Qs]
        out['candidate_gap']=float(eps.value)
        out['nullspace_transforms']=[E.tolist() for E in kernel.transforms]
    SCRATCH.mkdir(exist_ok=True)
    path=SCRATCH/f'sdp_n{n}_s{s}_k{modes}_grad{int(contact_gradients)}_strong{int(strong)}.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print('Saved',path,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--features',type=int,default=4)
    parser.add_argument('--modes',type=int,default=9);parser.add_argument('--samples',type=int,default=1800)
    parser.add_argument('--iterations',type=int,default=3);parser.add_argument('--no-contact-gradients',action='store_true')
    parser.add_argument('--strong-domain',action='store_true')
    parser.add_argument('--n',type=int,default=10)
    args=parser.parse_args();attempt(args.features,args.modes,args.samples,args.iterations,not args.no_contact_gradients,args.strong_domain,args.n)
