from pathlib import Path
import sys,itertools
ROOT=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
from flint import arb,ctx
from arb_bernstein import coefficients
from physical_boundary_bernstein import correction_upper,exact_q_coefficients
ctx.prec=160;degrees=(17,17,7);checked=0
for box in [(0.,2.,0.,2.,-2.,2.),(.99,1.01,.98,1.02,-1.01,-.99),(0.,.001,.99,1.01,-.01,.01),(1.11,1.12,1.11,1.12,-1.12,-1.11)]:
    exact_q,exact_q2=exact_q_coefficients(box,degrees)
    for lam in [0.,.125,1.7,12.]:
        for kappa in [0.,1.,4.]:
            expected=coefficients({(1,1,0):arb(lam),(0,0,2):-arb(lam),(2,2,0):-arb(kappa),(1,1,2):2*arb(kappa),(0,0,4):-arb(kappa)},degrees,box)
            upper=correction_upper(box,degrees,lam,kappa)
            assert all(arb(a)>=b for a,b in zip(upper.flat,expected.flat))
            assert all((arb(lam)*q-arb(kappa)*q2).overlaps(b)for q,q2,b in zip(exact_q.flat,exact_q2.flat,expected.flat))
            checked+=len(upper.flat)
print('Verified physical-boundary correction upper coefficients:',checked)
