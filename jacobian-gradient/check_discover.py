"""Independent numeric comparison with the old discovery derivative formulas.

Run in the existing SDP environment; imports no Jacobian runtime.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json
import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    analytic=load(HERE/'analytic.py','exact_generated_gradient')
    old=load(ROOT/'candidate-kernel/discover_sdp.py','original_numerical_formulas')
    data=json.loads((ROOT/'twelve-point-kernel-refinement/kernel.json').read_bytes())
    kernel=old.Kernel.__new__(old.Kernel)
    kernel.n=12;kernel.mass=11;kernel.s=data['feature_order'];kernel.modes=data['harmonics_through']
    kernel.nodes=np.array([float.fromhex(q) for q in data['radial_nodes']])
    kernel.moments=np.array([float.fromhex(q) for q in data['radial_moment_centers']])
    kernel.transforms=[np.array([[float.fromhex(q) for q in row] for row in E]) for E in data['transforms']]
    kernel.logtarget=float.fromhex(data['B'])
    Qs=[np.array([[float.fromhex(q) for q in row] for row in Q]) for Q in data['Q_blocks']]
    residual=analytic.Residual()
    points=[[.25,2/3,-.2],[7/6,1.5,-.75],[.5264,121/72,.9033],[1.238802,1.680555,-1.135873]]
    rng=np.random.default_rng(20261108)
    for _ in range(80):
        u,v=rng.uniform(.05,121/72,2);t=rng.uniform(-.95,.95)*np.sqrt(u*v)
        if u+v-2*t>=1/3:points.append([u,v,t])
    pts=np.array(points)
    _,_,base,D=kernel.coefficients(pts,True)
    original=base+sum(M@Q.reshape(-1) for M,Q in zip(D,Qs))
    generated=np.array([residual.jac(p) for p in pts])
    errors=np.max(np.abs(original-generated),axis=1)
    if max(errors)>1e-7:raise ArithmeticError('existing derivatives do not agree with generated exact derivative evaluation')
    starts=points[:4];results=[];U=121/72
    constraints=[{'type':'ineq','fun':lambda p:p[0]*p[1]-p[2]**2,
                  'jac':lambda p:np.array([p[1],p[0],-2*p[2]])},
                 {'type':'ineq','fun':lambda p:p[0]+p[1]-2*p[2]-1/3,
                  'jac':lambda p:np.array([1,1,-2])}]
    for start in starts:
        r=minimize(residual.value,start,jac=residual.jac,method='SLSQP',bounds=[(0,U),(0,U),(-U,U)],
            constraints=constraints,options={'ftol':1e-12,'maxiter':200})
        results.append({'start':start,'success':bool(r.success),'point':r.x.tolist(),
                        'h_minus_log_distance':float(r.fun),'evaluations':r.nfev,'iterations':r.nit})
    report={'passed':True,'scope':'numeric diagnostic only; no exhaustive minimum or new bound',
        'source_sha256':hashlib.sha256((ROOT/'candidate-kernel/discover_sdp.py').read_bytes()).hexdigest(),
        'points_compared':len(points),'maximum_gradient_disagreement':float(max(errors)),
        'old_gradient_formula_bug_detected':False,
        'old_attempt_slsqp_supplied_analytic_jacobian':False,
        'analytic_gradient_slsqp_examples':results}
    (HERE/'numeric_crosscheck.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
