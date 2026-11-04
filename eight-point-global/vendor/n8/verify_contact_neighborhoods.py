"""Rigorous local convexity certificates for the exact auxiliary contacts.

Full-domain coverage is handled separately by verify_global_cover.py.
"""
from pathlib import Path
import sys,json,time,math,hashlib,collections
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data'
sys.path.insert(0,str(ROOT.parent.parent/'friedman_full_audit_20261023/deps'))
from flint import arb,ctx,fmpq
ctx.prec=128
source=OUT/'global_pair_bound_exact_contacts.json';data=source.read_bytes();sha=hashlib.sha256(data).hexdigest();raw=json.loads(data)
a=(2*arb.pi()/7).cos();cc=[a,2*a*a-1,-2*a*a-a+arb(1)/2];logs=[(2-2*c).log()for c in cc]
def fld(coefs):
    val=arb(0)
    for q in coefs:val=val*a+arb(fmpq(q))
    return val
def scalar(value):return fld(value[0])+sum((fld(c)*l for c,l in zip(value[1:],logs)),arb(0))
encoded={tuple(e['powers']):e['value']for e in raw['polynomial_terms']}
assert all(encoded[(j,i,k)]==value for (i,j,k),value in encoded.items())
P={tuple(e['powers']):scalar(e['value'])for e in raw['polynomial_terms']}
def diff(poly,i):
    result={}
    for e,c in poly.items():
        if e[i]:ee=list(e);ee[i]-=1;result[tuple(ee)]=c*e[i]
    return result
G=[diff(P,i)for i in range(3)];H=[[diff(G[i],j)for j in range(3)]for i in range(3)]
def tree(poly,axis=0):
    if axis==3:return next(iter(poly.values()),arb(0))
    deg=max((e[axis]for e in poly),default=0)
    return [tree({e:c for e,c in poly.items()if e[axis]==j},axis+1)for j in range(deg+1)]
def evaluate(tr,xx,axis=0):
    if axis==3:return tr
    val=arb(0)
    for sub in reversed(tr):val=val*xx[axis]+evaluate(sub,xx,axis+1)
    return val
PT=tree(P);GT=list(map(tree,G));HT=[[tree(h)for h in row]for row in H]
dd=[1,1,-2]
def shifted(poly,point):
    out={}
    for e,c in poly.items():
        for i in range(e[0]+1):
            for j in range(e[1]+1):
                for k in range(e[2]+1):
                    ix=(i,j,k);v=c
                    for ei,ii,x in zip(e,ix,point):v*=math.comb(ei,ii)*x**(ei-ii)
                    out[ix]=out.get(ix,arb(0))+v
    return out
def positive_definite(M):
    size=len(M);L=[[arb(int(i==j))for j in range(size)]for i in range(size)];D=[]
    for j in range(size):
        d=M[j][j]-sum((L[j][i]**2*D[i]for i in range(j)),arb(0))
        if not d>0:return False,[]
        D.append(d)
        for k in range(j+1,size):L[k][j]=(M[k][j]-sum((L[k][i]*L[j][i]*D[i]for i in range(j)),arb(0)))/d
    return True,D

# Each contact has an EXACT value and the stated zero derivatives; their
# algebraic identities are independently checked by reconstruct_pair_bound_exact.py.
contacts=[[arb(0),arb(1),arb(0)]]+[[arb(1),arb(1),c]for c in cc]
neighborhoods=[];local_records=[]
for j,p in enumerate(contacts):
    d=sum((x*y for x,y in zip(p,dd)),arb(0))
    if j==0:
        inward=evaluate(GT[0],p)-1/d
        assert inward>0
    hs=[[tree(shifted(H[i][k],p))for k in range(3)]for i in range(3)]
    for power in range(4,22):
        eps=2.0**(-power);zz=[arb(0,eps)]*3
        di=d+sum((z*w for z,w in zip(zz,dd)),arb(0))
        if not di>0:continue
        matrix=[[evaluate(hs[i][k],zz)+dd[i]*dd[k]/di**2 for k in range(3)]for i in range(3)]
        ok,piv=positive_definite(matrix)
        if ok:break
    assert ok
    neighborhoods.append((p,eps))
    local_records.append({'contact':j,'radius':eps,'Hessian_LDL_pivots':[str(x)for x in piv]})
    print('LOCAL CERTIFIED',j,'radius',eps,flush=True)
neighborhoods.append(([arb(1),arb(0),arb(0)],neighborhoods[0][1]))
(OUT/'global_pair_bound_local_certificate.json').write_text(json.dumps({'coefficient_sha256':sha,'local_positivity_proved':True,'records':local_records,'global_majorant_proved':False},indent=2))

print("All contact neighborhoods certified; full-domain coverage is verified separately.",flush=True)
