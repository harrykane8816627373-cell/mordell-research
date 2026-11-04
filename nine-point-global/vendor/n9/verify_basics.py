"""Exact basic identities, candidate stationarity, spectrum, and source polynomial."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parent
if not __debug__: raise RuntimeError('Run without -O: assertions are proof checks.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
import sympy as s
start=time.monotonic();t,u=s.symbols('t u');n=9;lam=s.Integer(4)
a=s.sqrt(2)/2;roots=[1,a*(1+s.I),s.I,a*(-1+s.I),-1,a*(-1-s.I),-s.I,a*(1-s.I)]
z=[s.expand(3*s.sqrt(2)/4*x)for x in roots]+[s.Integer(0)]
assert s.simplify(sum(z))==0
assert s.simplify(sum(x*s.conjugate(x)for x in z))==9
from sympy.polys.matrices import DomainMatrix
F=s.QQ.algebraic_field(s.sqrt(2),s.I)
zz=[F.from_sympy(x)for x in z];zc=[F.from_sympy(s.conjugate(x))for x in z]
for i,x in enumerate(zz):assert sum((F.one/(x-y)for j,y in enumerate(zz)if i!=j),F.zero)==4*zc[i]
def build_A(points):
    mat=[[F.zero for j in range(9)]for i in range(9)]
    for i in range(9):
        for j in range(i):
            val=F.one/(points[i]-points[j])**2
            mat[i][i]+=val;mat[j][j]+=val;mat[i][j]=-val;mat[j][i]=-val
    return DomainMatrix(mat,(9,9),F)
AA=build_A(zz).matmul(build_A(zc));coeff=AA.charpoly()
char=s.factor(sum(F.to_sympy(c)*t**(9-j)for j,c in enumerate(coeff)))
print('A Abar characteristic',char,flush=True)
eigs=s.roots(char,t);assert sum(eigs.values())==9
assert eigs.get(s.Integer(16),0)==1 and all(ev<16 for ev in eigs if ev!=16)
P=s.Poly(384*t**3-1024*t**2+765*t-741,t)
H=128*(t**3+(u-t*t)**2/(9-t))-360*u+2025
assert s.cancel((9-t)*H.subs(u,12)-3*P.as_expr())==0
threshold=t*t+s.Rational(45,32)*(9-t)-12
assert threshold.subs(t,1)==s.Rational(1,4)
assert s.diff(threshold,t).subs(t,1)>0 and s.diff(threshold,t,2)>0
assert P.count_roots(1,9)==1
alpha=[(lo,hi)for (lo,hi),mult in P.intervals(eps=s.Rational(1,10**20))if lo>1][0]
U=s.Rational(1.9182004928588867)
assert U>s.Rational(8,9)*alpha[1]
assert s.Rational(9,32)*s.Rational(8,9)==s.Rational(1,4)
assert 128*s.Rational(9,8)**3*s.Rational(8,729)==2
assert 360*s.Rational(9,8)**2*s.Rational(8,729)==5
assert 2025*s.Rational(8,729)==s.Rational(200,9)
aa=s.Rational(71,40);vv=s.Rational(256,27)
pair_threshold=2*aa**2+s.Rational(5,4)*(8-2*aa)
pair_test=2*(2*aa**3+(vv-2*aa**2)**2/(8-2*aa))-5*vv+s.Rational(200,9)
assert pair_threshold>vv and pair_test>1
ratio_coeff=[1,-3,-8,79,-181,-426,2941,-5791,5528,-2876,795,-88,-4,1]
R=s.Poly(sum(c*t**k for k,c in enumerate(ratio_coeff)),t)
assert R.count_roots(1,s.oo)==1
ints=R.intervals(eps=s.Rational(1,10**24));assert len(ints)>=3
lo,hi=ints[2][0];assert s.Rational(6600,1000)<lo<hi<s.Rational(6601,1000)
assert s.Rational(6601,1000)<4+2*s.sqrt(2)
C=s.Rational(9**36,8**28)
rec={'global_M9_determined':False,'candidate_stationary_exact':True,'candidate_strict_local_maximum_mod_rotation':True,
 'candidate_value':str(C),'candidate_value_approx':float(C),'A_Abar_characteristic':str(char),'A_Abar_eigenvalues':{str(k):v for k,v in eigs.items()},
 'all_maxima_radius_squared_polynomial':str(P.as_expr()),'alpha_interval':list(map(str,alpha)),
 'rescaled_radius_squared_bound':str(U),'rescaled_min_distance_squared':'1/4',
 'joint_radius_moment_scaling_verified':True,'excluded_equal_squared_radii':'71/40','excluded_pair_FV':str(pair_test),
 'cantrell_ratio_polynomial':str(R.as_expr()),'cantrell_third_real_root_interval':[str(lo),str(hi)],'cantrell_ratio_squared_approx':float((lo+hi)/2),
 'source_coordinate_identity_scope':'Verified separately by derive_contact_coordinates.py and certify_contact_configuration.py; see FRIEDMAN_N9_EXCLUSION_PROOF.txt.',
 'source_url':'https://groups.google.com/g/sci.math.research/c/UFEX1TnyO_0','seconds':time.monotonic()-start}
(ROOT/'basic_exact_certificate.json').write_text(json.dumps(rec,indent=2))
print(json.dumps(rec,indent=2),flush=True)
