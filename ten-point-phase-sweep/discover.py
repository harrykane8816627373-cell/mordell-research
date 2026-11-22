"""Find rational references and phase cells for the exact replay."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

from verify_all import AngularBox, HERE, MODES


def discover(mode, output):
    engine=AngularBox()
    data=engine.setup_mode(mode)
    positions=np.array([complex(*[float(F(t.lo+t.hi,2*engine.scale)) for t in z]) for z in data[0]])
    outer=positions[2:]
    indices=np.triu_indices(10,1)
    basis=np.vstack((np.eye(3),-np.ones((1,3))))
    directions=np.tile(basis,(2,1))

    def slopes(center):
        phase=(1+1j*float(center))/(1-1j*float(center))
        def loss(q):
            points=np.r_[positions[:2],outer*phase*np.exp(1j*(directions@q))]
            d=abs(points[:,None]-points[None,:])**2
            np.fill_diagonal(d,1)
            g=2*np.sum((points[:,None]*points[None,:].conj()).imag/d,axis=1)
            return -np.log(d[indices]).sum(), -directions.T@g[2:]
        fit=minimize(loss,np.zeros(3),jac=True,method='BFGS',options={'gtol':1e-10})
        return [F(round(float(a)*10**8),10**8) for a in np.tan(fit.x/2)]

    stack=[(-F(1),F(1),0)]
    leaves=[]
    while stack:
        lo,hi,depth=stack.pop()
        middle=(lo+hi)/2
        angles=slopes(middle)
        verdict=engine.accept(mode,data,lo,hi,angles)
        if verdict is not None:
            leaves.append({'interval':[str(lo),str(hi)],'corrections':list(map(str,angles)), 'result':verdict})
        elif depth==20:
            raise ArithmeticError(f'phase interval {float(lo)}, {float(hi)} is still unresolved')
        else:
            stack.extend([(middle,hi,depth+1),(lo,middle,depth+1)])
    output.write_text(json.dumps(leaves,indent=2)+'\n',encoding='ascii',newline='\n')
    print(mode,len(leaves),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=MODES)
    parser.add_argument('--tolerance',type=F)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.tolerance is not None:
        MODES[args.mode]=args.tolerance
    discover(args.mode,args.output or HERE/f'{args.mode}.json')
