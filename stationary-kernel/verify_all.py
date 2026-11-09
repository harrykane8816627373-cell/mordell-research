"""Exact Ward residual and zero-sum-kernel identities; no new global bound."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,random

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
ZERO=(F(0),F(0));ONE=(F(1),F(0))


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,s):return (a[0]*s,a[1]*s)
def conj(a):return (a[0],-a[1])
def norm(a):return a[0]**2+a[1]**2
def inv(a):
    d=norm(a);need(d>0,'nonzero complex denominator')
    return (a[0]/d,-a[1]/d)
def csum(values):
    result=ZERO
    for a in values:result=add(result,a)
    return result
def power(z,k):
    result=ONE
    for _ in range(k):result=mul(result,z)
    return result


def quadratic(M,A):
    return csum(scale(mul(conj(M[i]),M[j]),A[i][j])
                for i in range(5) for j in range(5))


def real_bilinear(a,A,b):
    return sum((a[i]*A[i][j]*b[j] for i in range(5) for j in range(5)),F(0))


def rational(rng):return F(rng.randrange(-9,10),rng.randrange(1,8))
def complex_rational(rng):return (rational(rng),rational(rng))


def correction_matrices(n,rng):
    A1=[[F(0)]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i,5):
            if i<2:A1[i][j]=A1[j][i]=rational(rng)
    matrices=[]
    for v in ([F(1,n),F(1),F(0),F(0),F(0)],
              [-F(n-3,n),F(1),F(0),F(0),F(0)],
              [-F(n-4,n),F(1),F(0),F(0),F(0)]):
        a=[rational(rng) for _ in range(5)]
        matrices.append([[v[i]*a[j]+a[i]*v[j] for j in range(5)] for i in range(5)])
    return [matrices[0],A1,matrices[1],matrices[2]]


def ward_residual(points):
    n=len(points);lam=F(n,2)
    residual=[sub(csum(inv(sub(z,w)) for j,w in enumerate(points) if j!=i),scale(conj(z),lam))
              for i,z in enumerate(points)]
    p=[csum(power(z,k) for z in points) for k in range(8)]
    need(csum(residual)==scale(conj(p[1]),-lam),'sum force residual equals minus lambda conjugate centroid')
    for k in range(1,7):
        left=csum(mul(power(z,k+1),f) for z,f in zip(points,residual))
        mixed=csum(scale(power(z,k),norm(z)) for z in points)
        convolution=csum(mul(p[r],p[k-r]) for r in range(1,k))
        right=sub(add(scale(p[k],F(n)-F(k+1,2)),scale(convolution,F(1,2))),scale(mixed,lam))
        need(left==right,'exact Ward identity with explicit nonzero force residual')
    return 6


def phi(u):return [(u-1)**j for j in range(5)]


def kernel(z,w,k,A):
    phase=power(mul(z,conj(w)),k)[0]
    return phase*real_bilinear(phi(norm(z)),A,phi(norm(w)))


def moment_sum_identity(points,matrices):
    count=0;n=len(points)
    for k,A in enumerate(matrices):
        M=[csum(scale(power(z,k),(norm(z)-1)**j) for z in points) for j in range(5)]
        direct=sum((kernel(z,w,k,A) for z in points for w in points),F(0))
        quad=quadratic(M,A)
        need(quad[1]==0 and direct==quad[0],'ordered kernel sum equals exact complex moment quadratic form')
        diagonal=[kernel(z,z,k,A) for z in points]
        hchange=sum((-kernel(z,w,k,A)-(diagonal[i]+diagonal[j])/F(2*(n-1))
                     for i,z in enumerate(points) for j,w in enumerate(points[:i])),F(0))
        need(hchange==-direct/2,'diagonal-corrected two-point change sums to minus half the ordered correction')
        count+=1
    return count


def null_moment_identity(n,matrices,rng):
    # Abstract exact vectors satisfying only the stated linear constraints.
    # No claim that every such vector has a physical point representation.
    M0=[(F(n),F(0)),(-F(1),F(0))]+[(rational(rng),F(0)) for _ in range(3)]
    M1=[ZERO,ZERO]+[complex_rational(rng) for _ in range(3)]
    z2=complex_rational(rng);z3=complex_rational(rng)
    M2=[z2,scale(z2,F(n-3,n))]+[complex_rational(rng) for _ in range(3)]
    M3=[z3,scale(z3,F(n-4,n))]+[complex_rational(rng) for _ in range(3)]
    for M,A in zip((M0,M1,M2,M3),matrices):
        need(quadratic(M,A)==ZERO,'linear-constraint moment vector annihilates its correction form')
    return 4


def matrix_parameter_dimensions():
    # Flatten exact rational symmetric matrices to independently check rank.
    dimensions=[]
    for which in (1,2,3,0):
        basis=[]
        if which==1:
            for i in range(2):
                for j in range(i,5):
                    A=[[F(0)]*5 for _ in range(5)];A[i][j]=A[j][i]=F(1)
                    basis.append([A[r][s] for r in range(5) for s in range(5)])
        else:
            v=([F(1,12),F(1),F(0),F(0),F(0)] if which==0 else
               [-F(3,4),F(1),F(0),F(0),F(0)] if which==2 else
               [-F(2,3),F(1),F(0),F(0),F(0)])
            for index in range(5):
                a=[F(int(i==index)) for i in range(5)]
                basis.append([v[r]*a[s]+a[r]*v[s] for r in range(5) for s in range(5)])
        rank=0;rows=[r[:] for r in basis]
        for col in range(25):
            pivot=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
            if pivot is None:continue
            rows[rank],rows[pivot]=rows[pivot],rows[rank]
            lead=rows[rank][col];rows[rank]=[x/lead for x in rows[rank]]
            for i in range(len(rows)):
                if i!=rank and rows[i][col]:
                    factor=rows[i][col];rows[i]=[a-factor*b for a,b in zip(rows[i],rows[rank])]
            rank+=1
        need(rank==len(basis),'free correction parameters are linearly independent')
        dimensions.append(rank)
    need(dimensions==[9,5,5,5],'nineteen Ward and five optional mass parameters')
    return dimensions


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args();manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed own input '+name)
    for name,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'fixed source dependency '+name)
    rng=random.Random(2026110819);ward=ordered=nulls=cases=0
    for n in range(2,10):
        for trial in range(2):
            points=[]
            while len(points)<n:
                z=complex_rational(rng)
                if z not in points:points.append(z)
            matrices=correction_matrices(n,rng)
            ward+=ward_residual(points);ordered+=moment_sum_identity(points,matrices)
            nulls+=null_moment_identity(n,matrices,rng);cases+=1
    for n in range(2,22):
        for trial in range(5):nulls+=null_moment_identity(n,correction_matrices(n,rng),rng)
    dimensions=matrix_parameter_dimensions()
    # An allowed A1 can be negative on a one-point Gram matrix; this prevents
    # falsely labelling the full stationary correction as a PSD kernel.
    bad=[[F(0)]*5 for _ in range(5)];bad[0][0]=-F(1)
    need(kernel(ONE,ONE,1,bad)==-1,'explicit correction Gram matrix is not PSD')
    report={'passed':True,'date':'2026-11-08','arithmetic':'Python standard-library exact Gaussian rationals',
        'scope':'existing full-force Ward identities and new zero-sum correction architecture',
        'Ward_identities_claimed_original_here':False,
        'arbitrary_nonstationary_point_sets_checked':cases,
        'exact_force_residual_identities_checked':ward,
        'ordered_sum_and_diagonal_correction_identities_checked':ordered,
        'abstract_constraint_vector_zero_form_checks':nulls,
        'n12_free_Ward_parameter_dimensions':dimensions[:3],
        'n12_total_free_Ward_parameters':sum(dimensions[:3]),
        'optional_mass_only_free_parameters':dimensions[3],
        'total_with_optional_mass_parameters':sum(dimensions),
        'n12_v2':['-3/4','1','0','0','0'],'n12_v3':['-2/3','1','0','0','0'],
        'optional_n12_v0':['1/12','1','0','0','0'],
        'zero_ordered_sum_proved_for_all_complete_stationary_configurations':True,
        'optional_radial_zero_sum_proved_for_all_mass_n_minus_one_configurations':True,
        'arbitrary_configurations_falsely_assumed_stationary':False,
        'full_correction_kernel_positive_semidefinite_claimed':False,
        'explicit_non_PSD_correction_example_checked':True,
        'new_global_numerical_bound_proved':False,
        'formal_proof_assistant_verification_claimed':False,
        'input_sha256':manifest['files'],'dependency_sha256':manifest['dependencies']}
    if args.record_report:
        (HERE/'verification_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
