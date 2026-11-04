"""Rigorous PSD and contact-neighborhood checks. No global coverage claim."""
from pathlib import Path
import sys,json,time,math,hashlib
ROOT=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
from flint import arb,ctx,fmpq
ctx.prec=192;start=time.monotonic()
source=ROOT/'global_pair_bound_exact_contacts.json';blob=source.read_bytes();raw=json.loads(blob)
history=ROOT/'auxiliary_history';history.mkdir(exist_ok=True)
snapshot=history/(hashlib.sha256(blob).hexdigest()+'.json')
if not snapshot.exists():snapshot.write_bytes(blob)
a=arb(2).sqrt();logs=[(2-a).log(),arb(2).log()]
def fld(co):
    out=arb(0)
    for q in co:out=out*a+arb(fmpq(q))
    return out
def scalar(value):return fld(value[0])+sum((fld(q)*l for q,l in zip(value[1:],logs)),arb(0))
def pd(M):
    n=len(M);L=[[arb(int(i==j))for j in range(n)]for i in range(n)];D=[]
    for j in range(n):
        d=M[j][j]-sum((L[j][k]*L[j][k]*D[k]for k in range(j)),arb(0))
        if not d>0:return False,D
        D.append(d)
        for i in range(j+1,n):L[i][j]=(M[i][j]-sum((L[i][k]*L[j][k]*D[k]for k in range(j)),arb(0)))/d
    return True,D
dim=raw['matrix_dimension'];K=raw.get('angular_maximum',7);Q=[[[arb(0)for j in range(dim)]for i in range(dim)]for k in range(K+1)]
for item in raw['matrix_entries']:
    k,i,j=item['label'];Q[k][i][j]=Q[k][j][i]=scalar(item['value'])
qrecords=[]
for k,M in enumerate(Q):
    ok,piv=pd(M);assert ok,('Matrix not certified positive definite',k,list(map(str,piv)))
    qrecords.append(list(map(str,piv)))
print('All',K+1,'matrices rigorously positive definite',flush=True)
P={tuple(e['powers']):scalar(e['value'])for e in raw['polynomial_terms']}
encoded={tuple(e['powers']):e['value']for e in raw['polynomial_terms']}
assert all(encoded[(j,i,k)]==value for (i,j,k),value in encoded.items())
def diff(poly,i):
    out={}
    for e,c in poly.items():
        if e[i]:ee=list(e);ee[i]-=1;out[tuple(ee)]=c*e[i]
    return out
G=[diff(P,i)for i in range(3)];H=[[diff(G[i],j)for j in range(3)]for i in range(3)]
def tree(poly,axis=0):
    if axis==3:return next(iter(poly.values()),arb(0))
    degree=max((e[axis]for e in poly),default=0)
    return [tree({e:c for e,c in poly.items()if e[axis]==j},axis+1)for j in range(degree+1)]
def evaluate(tr,x,axis=0):
    if axis==3:return tr
    out=arb(0)
    for sub in reversed(tr):out=out*x[axis]+evaluate(sub,x,axis+1)
    return out
def shift(poly,p):
    out={}
    for e,c in poly.items():
        for i in range(e[0]+1):
            for j in range(e[1]+1):
                for k in range(e[2]+1):
                    ex=(i,j,k);v=c
                    for ei,ii,x in zip(e,ex,p):v*=math.comb(ei,ii)*x**(ei-ii)
                    out[ex]=out.get(ex,arb(0))+v
    return out
GT=list(map(tree,G));dd=[1,1,-2]
contacts=[[arb(0),arb(1),arb(0)]]+[[arb(1),arb(1),c]for c in [a/2,arb(0),-a/2,arb(-1)]]
records=[];Hq=[[0,1,0],[1,0,0],[0,0,-2]]
for j,p in enumerate(contacts):
    d=sum((x*y for x,y in zip(p,dd)),arb(0));tau=arb(0)
    if j==0:
        inward=evaluate(GT[0],p)-1/d;assert inward>0
        tau=inward
    if j==4:
        tau=(evaluate(GT[2],p)+2/d)/2;assert tau>0
    # Congruence before interval evaluation preserves the exact u/v cancellations.
    # The matrix is invertible, so positivity is equivalent to positivity of the
    # original Hessian; the surrounding cube remains in the original coordinates.
    if j==0:B=[[1,0,0],[0,1,0],[0,0,1]]
    elif j==4:
        # The last direction (1,1,-1) is tangent to uv-t²=0 at (1,1,-1).
        # This integer shear has determinant -2 and cancels large normal terms
        # before interval evaluation, while preserving the original cube.
        B=[[1,1,1],[-1,1,1],[0,0,-1]]
    else:B=[[1,1,0],[-1,1,0],[0,0,1]]
    transformed=[]
    for i in range(3):
        row=[]
        for k in range(3):
            poly={}
            for aa in range(3):
                for bb in range(3):
                    factor=B[aa][i]*B[bb][k]
                    if factor:
                        for ex,value in H[aa][bb].items():poly[ex]=poly.get(ex,arb(0))+factor*value
            row.append(tree(shift(poly,p)))
        transformed.append(row)
    HS=transformed
    db=[sum(B[a][i]*dd[a]for a in range(3))for i in range(3)]
    qb=[[sum(B[a][i]*Hq[a][b]*B[b][k]for a in range(3)for b in range(3))for k in range(3)]for i in range(3)]
    def corrected(eps,kappa):
        zz=[arb(0,eps)]*3;x=[c+z for c,z in zip(p,zz)]
        di=d+sum((z*w for z,w in zip(zz,dd)),arb(0))
        if not di>0:return None,False
        matrix=[[evaluate(HS[i][k],zz)+db[i]*db[k]/di**2 for k in range(3)]for i in range(3)]
        if j in [0,4]:
            # Multiplication gives a rigorous enclosure even for a zero-centered ball.
            # This python-flint version can return NaN for such a ball raised with **2.
            q=x[0]*x[1]-x[2]*x[2];grad=[x[1],x[0],-2*x[2]]
            grad=[sum((B[a][i]*grad[a]for a in range(3)),arb(0))for i in range(3)]
            if kappa and not tau-kappa*q>0:return None,False
            matrix=[[matrix[i][k]-tau*qb[i][k]+2*kappa*(grad[i]*grad[k]+q*qb[i][k])for k in range(3)]for i in range(3)]
        return matrix,True
    kappas=[0]if j not in [0,4] else[0]+[2**k for k in range(16)]
    certified=False
    for kappa in kappas:
        at_point,valid=corrected(0,kappa)
        if not valid or not pd(at_point)[0]:continue
        for power in range(4,41):
            eps=2.**(-power);mat,valid=corrected(eps,kappa)
            if not valid:continue
            ok,piv=pd(mat)
            if ok:certified=True;break
        if certified:break
    if not certified:
        diagnostic,_=corrected(0,0)
        failure={'coefficient_sha256':hashlib.sha256(blob).hexdigest(),'all_matrices_positive_definite':True,
                 'matrix_LDL_pivots':qrecords,'local_positivity_proved':False,'failed_contact':j,
                 'congruent_corrected_point_Hessian':[[str(v)for v in row]for row in diagnostic],'congruence_matrix':B,
                 'already_certified_contacts':records,'global_majorant_proved':False}
        negative=diagnostic[0][0]+diagnostic[1][1]-2*diagnostic[0][1]
        if j in [1,2,3] and negative<0:
            failure['strict_negative_direction']=[B[a][0]-B[a][1]for a in range(3)]
            failure['direction_in_congruent_basis']=[1,-1,0]
            failure['strict_negative_direction_second_derivative']=str(negative)
        (ROOT/'global_pair_bound_local_failure.json').write_text(json.dumps(failure,indent=2))
        raise RuntimeError(('Uncertified contact neighborhood',j))
    rec={'contact':j,'radius':eps,'Hessian_LDL_pivots':list(map(str,piv)),'kappa':kappa,'congruence_matrix':B}
    if j==0:rec['inward_derivative']=str(inward)
    if j in [0,4]:rec['tau']=str(tau);rec['method']='g-tau*(uv-t²)+kappa*(uv-t²)² has positive Hessian; tau-kappa*(uv-t²)>0.'
    records.append(rec);print('LOCAL CERTIFIED',j,'radius',eps,'kappa',kappa,flush=True)
rec={'coefficient_sha256':hashlib.sha256(blob).hexdigest(),'all_matrices_positive_definite':True,
     'matrix_LDL_pivots':qrecords,'local_positivity_proved':True,'records':records,'global_majorant_proved':False,'seconds':time.monotonic()-start}
(ROOT/'global_pair_bound_local_preconditioned_certificate.json').write_text(json.dumps(rec,indent=2))
print('Local certificates saved; no global coverage claim.',flush=True)
