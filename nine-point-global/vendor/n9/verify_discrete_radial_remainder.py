"""Exact algebra for the remaining-radii bound; not a global M9 proof."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
import sympy as s
from flint import arb,ctx,fmpq
ctx.prec=160
r,rho=s.symbols('r rho',real=True)
f=r*r*(r-rho)**2
assert s.expand(f-rho**2*r**2/4-r*r*(rho/2-r)*(3*rho/2-r))==0
assert s.expand(f-rho**2*(r-rho)**2/4-(r-rho)**2*(r-rho/2)*(r+rho/2))==0
xx=s.symbols('x:7')
assert s.expand(7*sum(x*x for x in xx)-sum(xx)**2-sum((xx[i]-xx[j])**2 for i in range(7)for j in range(i+1,7)))==0
rho=(225-arb(18465).sqrt())/90;L=arb(14372)/405-arb(20)/81*arb(18465).sqrt()
# The previously observed outer boundary violation, treated as exact binary input.
u=arb(.06347517146819519);v=arb(1.7736153273751138)
dist2=[(8-u-v-k*rho)**2 for k in range(8)]
minimum=min(arb(a.lower())for a in dist2)
excess=u*u*(u-rho)**2+v*v*(v-rho)**2+rho*rho*minimum/28-L
assert excess>0
out={'piecewise_polynomial_identities':True,'seven_variable_cauchy_identity':True,'old_boundary_point_excluded':True,'old_boundary_excess':str(excess),'global_M9_proved':False}
(ROOT/'discrete_radial_remainder_certificate.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
