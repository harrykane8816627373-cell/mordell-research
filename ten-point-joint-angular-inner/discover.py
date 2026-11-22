import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from tubes import EPS,S,Tube

HERE=Path(__file__).resolve().parent
def discover(mode,repo,out,eps):
    tube=Tube(repo)
    if eps: EPS[mode]=eps
    data=tube.data(mode)
    e=tube.eng
    pts=np.array([complex(*[float(F(t.lo+t.hi,2*e.scale)) for t in z]) for z in data[0]])
    outer=pts[2:]*np.sqrt(float(S)/sum(abs(pts[2:])**2))
    indices=np.triu_indices(10,1)
    basis=np.tile(np.vstack((np.eye(3),-np.ones((1,3)))),(2,1))
    def proposal(center):
        phase=(1+1j*float(center))/(1-1j*float(center))
        def loss(q):
            w=q[0];z=outer*phase*np.exp(1j*(basis@q[1:]))
            xs=np.r_[np.sqrt(w),-np.sqrt(w),z]
            d=abs(xs[:,None]-xs[None,:])**2
            np.fill_diagonal(d,1)
            g=2*np.sum((xs[:,None]*xs[None,:].conj()).imag/d,axis=1)
            gw=1/w-2*np.real(1/(z*z-w)).sum()-90/(float(S)+2*w)
            return 45*np.log((float(S)+2*w)/9)-np.log(d[indices]).sum(),-np.r_[gw,basis.T@g[2:]]
        fit=minimize(loss,np.r_[.12,np.zeros(3)],jac=True,method='L-BFGS-B',
            bounds=[(.107,.133)]+[(-.15,.15)]*3,options={'gtol':3e-12,'ftol':1e-15,'maxiter':60})
        return [F(round(float(a)*10**9),10**9) for a in np.r_[fit.x[0],np.tan(fit.x[1:]/2)]]
    stack=[(-F(1),F(1),0)]
    leaves=[]
    while stack:
        lo,hi,d=stack.pop()
        mid=(lo+hi)/2
        vals=proposal(mid)
        res=tube.accept(mode,data,lo,hi,vals)
        if res:
            leaves.append([str(lo),str(hi),*map(str,vals),res])
        elif d==18:
            raise ArithmeticError(f'phase {lo} {hi} unresolved, w0={vals[0]}, previous leaves={len(leaves)}')
        else:stack.extend([(mid,hi,d+1),(lo,mid,d+1)])
    out.write_text(json.dumps(leaves,separators=(',',':'))+'\n',encoding='ascii',newline='\n')
    print(mode,eps,len(leaves),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('mode',choices=EPS)
    p.add_argument('--repo',type=Path,default=HERE.parent)
    p.add_argument('--out',type=Path)
    p.add_argument('--eps',type=F)
    a=p.parse_args()
    discover(a.mode,a.repo.resolve(),a.out or HERE/(a.mode+'.json'),a.eps)
