"""Rigorous positive definiteness for the reconstructed auxiliary matrices only."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data'
sys.path.insert(0,str(ROOT.parent.parent/'friedman_full_audit_20261023/deps'))
from flint import arb,ctx,fmpq
ctx.prec=192
rawbytes=(OUT/'global_pair_bound_exact_contacts.json').read_bytes();raw=json.loads(rawbytes)
a=(2*arb.pi()/7).cos();cosines=[a,2*a*a-1,-2*a*a-a+arb(1)/2]
logs=[(2-2*c).log()for c in cosines]
def field(coefs):
    v=arb(0)
    for q in coefs:v=v*a+arb(fmpq(q))
    return v
def value(encoded):return field(encoded[0])+sum((field(c)*l for c,l in zip(encoded[1:],logs)),arb(0))
Q=[[[arb(0)for j in range(4)]for i in range(4)]for k in range(7)]
for entry in raw['matrix_entries']:
    k,i,j=entry['label'];Q[k][i][j]=Q[k][j][i]=value(entry['value'])
pivots=[]
for k,M in enumerate(Q):
    L=[[arb(int(i==j))for j in range(4)]for i in range(4)];D=[]
    for j in range(4):
        pivot=M[j][j]-sum((L[j][h]**2*D[h]for h in range(j)),arb(0))
        assert pivot>0,(k,j,str(pivot));D.append(pivot)
        for i in range(j+1,4):L[i][j]=(M[i][j]-sum((L[i][h]*L[j][h]*D[h]for h in range(j)),arb(0)))/D[j]
    pivots.append([str(d)for d in D])
polys=[{'powers':entry['powers'],'coefficient':str(value(entry['value'])),'approx':float(value(entry['value']).mid())}for entry in raw['polynomial_terms']]
result={'coefficient_sha256':hashlib.sha256(rawbytes).hexdigest(),'global_majorant_proved':False,'exact_auxiliary_matrices_positive_definite':True,
        'working_precision_bits':ctx.prec,'LDL_positive_pivot_balls':pivots,'polynomial_coefficient_balls':polys}
(OUT/'global_pair_bound_psd_certificate.json').write_text(json.dumps(result,indent=2))
print('All seven exact auxiliary matrices are positive definite.',flush=True)
print('Polynomial terms:',len(polys),'max total degree:',max(sum(x['powers'])for x in polys),flush=True)
print('Continuum majorant is NOT yet certified.',flush=True)
