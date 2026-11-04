"""Independently rebuild serialized auxiliary polynomials from the kernel matrices."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
import sympy as s
from sympy.polys.domains import QQ
blob=(ROOT/'global_pair_bound_exact_contacts.json').read_bytes();raw=json.loads(blob)
F=QQ.algebraic_field(s.sqrt(2));a=F.from_sympy(s.sqrt(2));zero=F.zero;one=F.one
R=F.poly_ring('u','v','t');u,v,t=R.gens
D=raw['radial_degree'];N=raw['matrix_dimension'];K=raw.get('angular_maximum',7)
assert N==D+1 and D>=2 and K>=7
def fld(co):
    q=zero
    for c in co:q=q*a+F.convert(s.Rational(c))
    return q
Q=[[[[zero for j in range(N)]for i in range(N)]for k in range(K+1)]for ell in range(3)]
expected={(k,i,j)for k in range(K+1)for i in range(N)for j in range(i,N)}
seen=set()
for item in raw['matrix_entries']:
    label=tuple(item['label']);assert label in expected and label not in seen;seen.add(label)
    assert len(item['value'])==3;k,i,j=label
    for ell,co in enumerate(item['value']):Q[ell][k][i][j]=Q[ell][k][j][i]=fld(co)
assert seen==expected
def feats(x,k):
    if k==0:return [x-R.convert(s.Rational(8,9))]+[x**j*(x-1)for j in range(1,D+1)]
    if k%8==0:return [(x-1)**j for j in range(1,D+2)]
    return [(x-1)**j for j in range(D+1)]
H=[R.one,t]
for k in range(2,K+1):H.append(2*t*H[-1]-u*v*H[-2])
P=[]
for ell in range(3):
    kernel=R.zero;du=R.zero;dv=R.zero
    for k in range(K+1):
        fu,fv=feats(u,k),feats(v,k)
        for i in range(N):
            for j in range(N):
                c=R.convert(Q[ell][k][i][j],F)
                kernel+=H[k]*c*fu[i]*fv[j]
                du+=u**k*c*fu[i]*fu[j];dv+=v**k*c*fv[i]*fv[j]
    base=R.convert(s.Rational(9,16))*(u+v)-1 if ell==0 else (R.convert(s.Rational(2,3))if ell==2 else R.zero)
    P.append(base-kernel-R.convert(s.Rational(1,16))*(du+dv))
encoded=[R.zero]*3;seen=set()
for item in raw['polynomial_terms']:
    powers=tuple(item['powers']);assert len(powers)==3 and powers not in seen;seen.add(powers)
    assert all(isinstance(e,int)and e>=0 for e in powers)and len(item['value'])==3
    monomial=u**powers[0]*v**powers[1]*t**powers[2]
    for ell,co in enumerate(item['value']):encoded[ell]+=R.convert(fld(co),F)*monomial
assert encoded==P
contacts=[(zero,one,zero)]+[(one,one,c)for c in [a/2,zero,-a/2,-one]]
logvalues=[(0,0,0),(0,1,0),(0,0,1),(0,-1,1),(0,0,2)]
def ev(poly,point):
    result=zero
    for ex,c in poly.to_dict().items():
        term=c
        for xx,k in zip(point,ex):term*=xx**k
        result+=term
    return result
for j,point in enumerate(contacts):
    d=point[0]+point[1]-2*point[2]
    directions=[(0,1,0),(0,0,1)]if j==0 else ([(1,0,s.Rational(-1,2)),(0,1,s.Rational(-1,2))]if j==4 else[(1,0,0),(0,1,0),(0,0,1)])
    for ell in range(3):
        assert ev(P[ell],point)==F.convert(logvalues[j][ell]),('contact value',j,ell)
        for direction in directions:
            actual=sum((F.convert(q)*ev(P[ell].diff(var),point)for q,var in zip(direction,(u,v,t))),zero)
            target=F.convert(direction[0]+direction[1]-2*direction[2])/d if ell==0 else zero
            assert actual==target,('contact derivative',j,ell,direction)
out={'coefficient_sha256':hashlib.sha256(blob).hexdigest(),'kernel_definition_exact':True,
     'all_serialized_entries_checked':True,'contact_values_and_gradients_exact':True,
     'matrix_count':K+1,'matrix_dimension':N,'global_majorant_proved':False}
(ROOT/'global_pair_bound_exact_definition_certificate.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
