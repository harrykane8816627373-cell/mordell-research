"""Independent exact check of the serialized kernel, polynomial, and contacts."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data'
sys.path.insert(0,str(ROOT.parent.parent/'friedman_full_audit_20261023/deps'))
import sympy as s
from sympy.polys.domains import QQ
rawbytes=(OUT/'global_pair_bound_exact_contacts.json').read_bytes();raw=json.loads(rawbytes)
aexpr=s.cos(2*s.pi/7);F=QQ.algebraic_field(aexpr);a=F.from_sympy(aexpr)
R=F.poly_ring('u','v','t');u,v,t=R.gens
zero=F.zero;one=F.one

def fld(q):
    ans=zero
    for c in q:ans=ans*a+F.convert(s.Rational(c))
    return ans
Q=[[[[zero for j in range(4)]for i in range(4)]for k in range(7)]for ell in range(4)]
for e in raw['matrix_entries']:
    k,i,j=e['label']
    for ell,c in enumerate(e['value']):Q[ell][k][i][j]=Q[ell][k][j][i]=fld(c)

def feats(x,k):return [x-R.convert(s.Rational(7,8)),x*(x-1),x*x*(x-1),x**3*(x-1)]if k==0 else [R.one,x-1,(x-1)**2,(x-1)**3]
H=[R.one,t]
for k in range(2,7):H.append(2*t*H[-1]-u*v*H[-2])
P=[]
for ell in range(4):
    kernel=R.zero;diag_u=R.zero;diag_v=R.zero
    for k in range(7):
        fu,fv=feats(u,k),feats(v,k)
        for i in range(4):
            for j in range(4):
                q=R.convert(Q[ell][k][i][j],F)
                kernel+=H[k]*q*fu[i]*fv[j]
                diag_u+=u**k*q*fu[i]*fu[j];diag_v+=v**k*q*fv[i]*fv[j]
    base=R.convert(s.Rational(4,7))*(u+v)-1 if ell==0 else R.convert(s.Rational(1,4))
    P.append(base-kernel-R.convert(s.Rational(1,14))*(diag_u+diag_v))
encoded=[]
for ell in range(4):
    pp=R.zero
    for e in raw['polynomial_terms']:
        pp+=R.convert(fld(e['value'][ell]),F)*u**e['powers'][0]*v**e['powers'][1]*t**e['powers'][2]
    encoded.append(pp);assert pp==P[ell],('kernel definition',ell)
cs=[a,2*a*a-one,-2*a*a-a+one/2]
contacts=[(zero,one,zero)]+[(one,one,c)for c in cs]
def ev(poly,p):
    ans=zero
    for ex,coeff in poly.to_dict().items():
        term=coeff
        for xx,k in zip(p,ex):term*=xx**k
        ans+=term
    return ans
for j,p in enumerate(contacts):
    d=p[0]+p[1]-2*p[2]
    for ell in range(4):
        target=one if j>0 and ell==j else zero
        assert ev(P[ell],p)==target,('contact value',j,ell)
        axes=[1,2]if j==0 else [0,1,2]
        for axis in axes:
            targ=(one if axis<2 else -2*one)/d if ell==0 else zero
            assert ev(P[ell].diff((u,v,t)[axis]),p)==targ,('contact derivative',j,ell,axis)
result={'coefficient_sha256':hashlib.sha256(rawbytes).hexdigest(),'kernel_definition_exact':True,'contact_values_and_gradients_exact':True,'global_majorant_proved':False}
(OUT/'global_pair_bound_exact_definition_certificate.json').write_text(json.dumps(result,indent=2))
print(result,flush=True)
