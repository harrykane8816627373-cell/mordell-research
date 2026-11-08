"""Exact independent audits for the Nov8 ring population obstructions; -O safe."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from pathlib import Path
import argparse, hashlib, json, time

from four_ring_rigidity import certificate, classify_four, isolated_force_exponents

HERE=Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


class Poly:
    def __init__(self, value=0):
        self.d=({k:F(v) for k,v in value.items() if v} if isinstance(value,dict)
                else ({():F(value)} if value else {}))
    def __add__(self,other):
        other=other if isinstance(other,Poly) else Poly(other)
        d=dict(self.d)
        for k,v in other.d.items():d[k]=d.get(k,F(0))+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self,other):return self+(-other if isinstance(other,Poly) else -Poly(other))
    def __rsub__(self,other):return Poly(other)-self
    def __mul__(self,other):
        other=other if isinstance(other,Poly) else Poly(other)
        d=defaultdict(F)
        for a,x in self.d.items():
            for b,y in other.d.items():d[tuple(sorted(a+b))]+=x*y
        return Poly(dict(d))
    __rmul__=__mul__
    def __pow__(self,k):
        out=Poly(1)
        for _ in range(k):out=out*self
        return out


def symbols(*names):return [Poly({(name,):1}) for name in names]


def universal_identities():
    a,b,m,s,A,B,D,X,Y,T,U,V,h,n= symbols(
        'a','b','m','s','A','B','D','X','Y','T','U','V','h','n')
    checked=[]
    def zero(p,label):need(not p.d,label);checked.append(label)
    # An arbitrary background of equal largest rings only changes complex T.
    zero(A*(a+T)-T*A-a*A,'one exceptional smaller ring')
    for q in (0,1):
        R0=T*A*B;Re=(a+b+T)*D**q
        zero(A*B*Re-D**q*R0-(a+b)*A*B*D**q,
             'two smaller rings: unbalanced isolated exponents q='+str(q))
    Ra=-(a+T)*B;Rb=-(b+T)*A
    zero(A*Ra-B*Rb-(b-a)*A*B,'two smaller rings: unequal balanced populations')
    def balanced(root,shift):
        return -(a+shift)*(A+B),(2*a+shift)*root+shift*A*B
    Rx0,Rx1=balanced(X,U);Ry0,Ry1=balanced(Y,V)
    combo=(A+B)*(Rx1-Ry1)+(X+A*B)*Rx0-(Y+A*B)*Ry0
    zero(combo-a*(A+B)*(X-Y),'two largest rings forced to the same root constant')
    zero((4*a+h)*(4*a+h-1)-4*a*(4*a-1+2*h)-h*(h-1),
         'three-ring exceptional case: binary-center virial identity')
    # Balanced four-ring exceptional pattern (s,s,s,3s).
    S1,S2,S3=symbols('S1','S2','S3')
    # Q=t^3-S1*t^2+S2*t-S3; T is the full largest-ring shift.
    c0=(3*s+T)*D-T*S3
    zero(c0-((2*s)*D+s*S3)-(T+s)*(D-S3),
         'S1=0 branch at T=-s has 2D+S3=0')
    zero(c0-(s*D+2*s*S3)-(T+2*s)*(D-S3),
         'S2=0 branch at T=-2s has D+2S3=0')
    zero((6*s+h)*(6*s+h-1)-6*s*(6*s-1+2*h)-h*(h-1),
         'four-ring exceptional case: binary-center virial identity')
    zero(3*(6*s-1+2*h)-F(3,2)*(7*s-1+2*h)-F(3,2)*(5*s-1+2*h),
         'second branch: exact mean-radius ratio from virial')
    zero(11*(7*s-1+2*h)-15*(5*s-1+2*h)-2*(s+2-4*h),
         'second branch: uniform 11/15 upper bound')
    zero(n-F(1,2)*(s+1)-F(1,2)*(2*n-s-1),'full Ward moment coefficient')
    # Family (a,s,s,2s), a != s and a < 2s.
    S=A+B
    for q in (0,1):
        Rs=D*(s+T)*S;Rother=-(a+s+T)*S*X**q
        zero(X**q*Rs+D*Rother+a*D*S*X**q,
             'two half-population rings: incompatible residues q='+str(q))
    need(F(11,15)**3<F(1,2),'strict rational cubed AMGM comparison')
    q,=symbols('q')
    c0=(6*q+T)*D**2-T*S3
    zero(c0-(4*q*D**2+2*q*S3)-(T+2*q)*(D**2-S3),
         'doubling rule: S1=0 branch at T=-2q')
    zero(c0-(2*q*D**2+4*q*S3)-(T+4*q)*(D**2-S3),
         'doubling rule: S2=0 branch at T=-4q')
    zero((9*q+h)*(9*q+h-1)-9*q*(9*q-1+2*h)-h*(h-1),
         'doubling rule: binary-center virial identity')
    zero(9*(9*q-1+2*h)-3*(11*q-1+2*h)-6*(8*q-1+2*h),
         'doubling rule: exact small-radius mean')
    zero(3*(11*q-1+2*h)-4*(8*q-1+2*h)-(q+1-2*h),
         'doubling rule: uniform 3/4 bound')
    need(F(3,4)**3<F(1,2),'doubling rule strict rational cube comparison')
    return checked


def expand_binomials(populations):
    constants=symbols('A','B','C')
    p={0:Poly(1)}
    for a,c in zip(populations,constants):
        out={}
        for e,q in p.items():
            out[e+a]=out.get(e+a,Poly())+q
            out[e]=out.get(e,Poly())-c*q
        p=out
    return p,constants


def balanced_coefficient_audit():
    checked=0
    for a,b,c in combinations_with_replacement(range(2,31),3):
        if a==c:continue
        p,(A,B,C)=expand_binomials((a,b,c));M=a+b+c
        low,high=(C*(A+B),-(A+B)) if a==b else (B*C,-A)
        need(not (p[a]-low).d and not (p[M-a]-high).d,
             'two distinct unavoidable intermediate coefficients in balanced exterior product')
        need(a!=M-a and 0<a<M-a<M,'intermediate degrees are distinct')
        checked+=1
    return checked


def subset_sum_audit():
    audited=0;boundary=0;examples=[]
    for t in range(1,5):
        for outside in combinations_with_replacement(range(2,17),t):
            # Independent mask expansion, not the implementation's doubling list.
            sums=[sum(a for a,selected in zip(outside,mask) if selected)
                  for mask in product((0,1),repeat=t)]
            if len(set(sums))!=len(sums):continue
            M=sum(outside)
            for m in range(max(max(outside)+1,(M+1)//2),41):
                groups=defaultdict(list)
                for exponent in sums:groups[exponent % m].append(exponent)
                single=sorted(v[0] for v in groups.values() if len(v)==1)
                need(len(single)>=2,'two incompatible full-force constants in every audited dissociated profile')
                need(single==isolated_force_exponents(outside,m),'independent isolated force coefficients')
                if M==2*m:
                    need(m not in sums,'at the degree boundary the midpoint subset is absent')
                    boundary+=1
                audited+=1
    return {'cases_checked':audited,'exterior_ring_counts':[1,2,3,4],
            'exterior_population_range':[2,16],'largest_population_maximum':40,
            'degree_equals_two_m_cases':boundary,'general_proof_replaced_by_finite_enumeration':False}


def ward_audit():
    checked=0
    for s in range(2,129):
        product_coefficients=defaultdict(int)
        for j in range(s+1):
            product_coefficients[(s-j+1,j)]+=1
            product_coefficients[(s-j,j+1)]-=1
        product_coefficients={k:v for k,v in product_coefficients.items() if v}
        need(product_coefficients=={(s+1,0):1,(0,s+1):-1},'exact divided difference behind the full Ward identity')
        checked+=1
    return checked


def equal_triplet_endpoint_audit():
    checked=0
    for s in range(2,129):
        for m in range(s+1,3*s+6):
            if m in (2*s,3*s) or 2*m==3*s:continue
            need(all(e % m for e in (s,2*s,3*s)),
                 'constant endpoint is isolated for every audited nonresonant outlier')
            need(all((3*s-e) % m for e in (0,s,2*s)),
                 'leading endpoint is isolated for every audited nonresonant outlier')
            checked+=1
    return checked


def census_audit(stored):
    counts=Counter();profiles=0
    # A separate fixed-total partition loop checks the coverage of the catalog.
    for total in range(8,101):
        for a in range(2,total//4+1):
            for b in range(a,(total-a)//3+1):
                for c in range(b,(total-a-b)//2+1):
                    m=total-a-b-c
                    need(m>=c,'complete sorted four-ring population partition')
                    counts[classify_four((a,b,c,m))]+=1;profiles+=1
    need(profiles==stored['census']['profiles_checked'] and dict(sorted(counts.items()))==stored['census']['counts'],
         'independent finite-profile coverage')
    for row in stored['examples']:
        p=row['populations'];need(classify_four(p)==row['case'],'fixed example conclusion')
    for p in [(2,2,2,4),(2,2,4,6),(2,4,4,6),(2,4,6,8)]:
        need(classify_four(p)=='not_excluded_by_these_theorems','unresolved families remain unresolved')
    return {'profiles_checked':profiles,'counts':dict(sorted(counts.items()))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args();start=time.monotonic()
    raw=(HERE/'source_manifest.json').read_bytes();manifest=json.loads(raw)
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed input '+name)
    stored=json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    need(stored==certificate(),'recomputed exact certificates')
    algebra=universal_identities();balanced=balanced_coefficient_audit()
    subset=subset_sum_audit();ward=ward_audit();census=census_audit(stored)
    endpoints=equal_triplet_endpoint_audit()
    report={'passed':True,'date':'2026-11-08','source_manifest_sha256':hashlib.sha256(raw).hexdigest(),
            'input_sha256':manifest['files'],'universal_polynomial_identities_checked':algebra,
            'balanced_exterior_coefficient_profiles_checked':balanced,
            'dissociated_subset_sum_audit':subset,'exact_ward_divided_differences_checked':ward,
            'three_equal_population_endpoint_cases_checked':endpoints,
            'four_ring_population_census':census,
            'arbitrary_continuous_phases_allowed':True,'full_pointwise_force_used':True,
            'all_four_ring_populations_equal_claimed':False,'remaining_patterns_certified_to_exist':False,
            'general_n_sharp_values_determined':False,'floating_point_used':False,
            'analytic_lemmas_not_formalized':['positive-power modulus injectivity',
                'three-vector radial weighted closure','three-variable AMGM','half-circle product sign'],
            'elapsed_seconds':round(time.monotonic()-start,3)}
    if args.record_report:
        (HERE/'verification_report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
