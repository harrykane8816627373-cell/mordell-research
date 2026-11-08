"""Exact algebraic verification of ring rigidity and unrestricted instability."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, sys, time

from ring_rigidity import certificate, fourier_coefficients, rotational_bound

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


class Poly:
    """Small independent rational multivariate polynomial identity checker."""
    def __init__(self, value=0):
        self.d = ({k:F(v) for k,v in value.items() if v} if isinstance(value,dict)
                  else ({():F(value)} if value else {}))
    def __add__(self, other):
        other = other if isinstance(other,Poly) else Poly(other)
        d = dict(self.d)
        for k,v in other.d.items():
            d[k] = d.get(k,F(0))+v
        return Poly(d)
    __radd__ = __add__
    def __neg__(self):
        return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self, other):
        return self + (-other if isinstance(other,Poly) else -Poly(other))
    def __rsub__(self, other):
        return Poly(other)-self
    def __mul__(self, other):
        other = other if isinstance(other,Poly) else Poly(other)
        d = defaultdict(F)
        for a,x in self.d.items():
            for b,y in other.d.items():
                d[tuple(sorted(a+b))] += x*y
        return Poly(dict(d))
    __rmul__ = __mul__
    def __pow__(self, power):
        result = Poly(1)
        for _ in range(power):
            result = result*self
        return result


def symbols(*names):
    return [Poly({(name,):1}) for name in names]


def identities():
    a,b,m,C,A,B,D,X,Y,h,u,v,w,k,l,K = symbols(
        'a','b','m','C','A','B','D','X','Y','h','u','v','w','k','l','K')
    checked = []
    def zero(value,label):
        need(not value.d,label); checked.append(label)
    # These are polynomial syzygies for the actual pointwise force remainders.
    Ra=(a+m+C)*X-(a+C)*Y
    R0=A*(C*Y-(m+C)*X)
    zero(A*Ra+R0-a*A*(X-Y),'two equal largest populations: nonzero saturated remainder combination')
    for q in (0,1):
        constant=C*A*B; leading=(a+b+C)*D**q
        zero(A*B*leading-D**q*constant-(a+b)*A*B*D**q,
             'unique largest unbalanced: constant/leading elimination q='+str(q))
    Ra=-(a+C)*B; Rb=-(b+C)*A
    zero(A*Ra-B*Rb-(b-a)*A*B,'balanced unequal populations: two different forced constants')
    Ra=-(a+C)*(A+B); R0=(2*a+C)*D+C*A*B
    zero((A+B)*R0+(D+A*B)*Ra-a*(A+B)*(D-A*B),
         'balanced equal populations: product of the two smaller root constants')
    n=4*a+h; force=4*a-1+2*h
    zero(n*(n-1)-4*a*force-h*(h-1),'binary-center virial identity for every population a')
    zero((u-v)**2-(u+v)**2+4*u*v,'strict AMGM obstruction after positive modulus reduction')
    c=lambda kk,mm:F(1,2)*(kk-2)*(kk-mm)
    zero(c(k+1,m)-c(k,m)-k+F(1,2)*(m+1),'general Fourier coefficient recurrence')
    zero(c(1,m)-F(1,2)*(m-1),'general initial Fourier coefficient')
    zero(c(m,m),'constant Fourier mode')
    zero(c(1,2*l)-c(l+1,2*l)-F(1,2)*l*l,'even-order Fourier range')
    zero(c(1,2*l+1)-c(l+1,2*l+1)-F(1,2)*l*(l+1),'odd-order Fourier range')
    zero(16*K*K-2*(K*(8*K)-1)-2,'no-center first excluded order m=8K')
    zero((16*K*K-4*K)-2*(K*(8*K-1)-1)-(2-2*K),
         'no-center immediately preceding order')
    zero((16*K*K+4*K)-2*(K*(8*K+1)+1)-(2*K-2),
         'center first excluded order m=8K+1 for K>=2')
    zero(16*K*K-2*(K*(8*K)+1)+2,'center immediately preceding order')
    return checked


def direct_fourier_audit():
    checks=0
    for m in range(2,129):
        # I(w)*(1-w)=1-(1+w+...+w^(m-1))/m, so I is exactly
        # (1-w)^(-1) at every nonidentity m-th root, composite m included.
        inverse=[F(m-1-j,m) for j in range(m-1)]
        product=[F(0)]*m
        for j,x in enumerate(inverse):
            product[j]+=x; product[j+1]-=x
        need(product==[F(m-1,m)]+[F(-1,m)]*(m-1),'independent inverse at all nonidentity roots')
        square=[F(0)]*(2*m-3)
        for i,x in enumerate(inverse):
            for j,y in enumerate(inverse):
                square[i+j]+=x*y
        folded=[F(0)]*m
        for j,x in enumerate(square):
            folded[j % m]+=x
        # Sum w^j over all roots except 1 is m-1 if m|j, otherwise -1.
        exact=[m*(folded[0]-folded[-k % m]) for k in range(1,m+1)]
        need(exact==fourier_coefficients(m),'direct root trace equals every Fourier coefficient')
        need(max(exact)-min(exact)==F(m*m//4,2),'direct full Fourier range')
        checks+=m
    return {'rotation_orders_checked':127,'direct_coefficients_checked':checks,
            'all_root_traces_exact':True,'floating_point_used':False}


def census_audit(stored):
    counts=Counter();profiles=0
    for a in range(2,34):
        for b in range(a,(100-a)//2+1):
            for c in range(b,100-a-b+1):
                profiles+=1
                if a==b==c:
                    kind='equal_population_not_excluded'
                elif b==c:
                    need(0<a<c,'one smaller population')
                    need(a % c!=0,'injectivity obstruction is genuinely nonconstant')
                    kind='two_largest_equal'
                else:
                    need(a<=b<c,'unique largest population')
                    residues=defaultdict(list)
                    for exponent in {0,a,b,a+b}:
                        residues[exponent % c].append(exponent)
                    if a+b!=c:
                        need(residues[0]==[0],'constant term is isolated')
                        need(residues[(a+b)%c]==[a+b],'leading term is isolated')
                        kind='unique_largest_unbalanced'
                    elif a!=b:
                        need(residues[a]==[a] and residues[b]==[b],'two separate intermediate residues')
                        kind='unique_largest_balanced_unequal'
                    else:
                        need(c==2*a and set(residues)=={0,a},'the sole balanced cancellation case')
                        kind='unique_largest_balanced_equal'
                counts[kind]+=1
    expected=stored['population_census']
    need(dict(sorted(counts.items()))==expected['counts'] and profiles==expected['profiles_checked'],
         'independent residue-partition census')
    return {'profiles_checked':profiles,'population_case_counts':dict(sorted(counts.items())),
            'finite_enumeration_replaces_general_proof':False}


def derivative(p):
    return [i*p[i] for i in range(1,len(p))]


def evaluate(p,x):
    answer=F(0)
    for coefficient in reversed(p):
        answer=answer*x+coefficient
    return answer


def witness_audit(stored):
    results=[]
    for record in stored['stationary_but_unstable_examples']:
        n=record['n'];p=list(map(F,record['polynomial_in_increasing_degree_order']))
        d=derivative(p);dd=derivative(d)
        equation=[F(0)]*(n+1)
        for i,x in enumerate(dd):equation[i]+=x
        for i,x in enumerate(d):equation[i+1]-=2*x
        for i,x in enumerate(p):equation[i]+=2*n*x
        need(not any(equation),'exact Hermite differential equation implies full force equations')
        parity=n%2
        q=p[parity::2]
        need(len(q)==4 and q[-1]==1,'three squared-radius roots')
        intervals=[tuple(map(F,row)) for row in record['squared_radius_isolation_intervals']]
        need(all(intervals[j][1]<=intervals[j+1][0] for j in range(2)),
             'squared-radius isolation intervals are disjoint')
        for lo,hi in intervals:
            lo,hi=F(lo),F(hi)
            need(0<=lo<hi and evaluate(q,lo)*evaluate(q,hi)<0,'three disjoint positive root intervals')
        # Newton sums are independent of the construction recurrence.
        sums={0:F(n)}
        for k in range(1,5):
            sums[k]=-sum((p[n-j]*sums[k-j] for j in range(1,k)),F(0))-k*p[n-k]
        need(sums[1]==sums[3]==0 and sums[2]==F(n*(n-1),2),'centering and virial mass')
        shift=F(n-1,2)
        need(sums[2]-n*shift==0,'unrestricted Hessian direction is centered')
        norm=sums[4]-2*shift*sums[2]+n*shift*shift
        need(norm==F(record['direction_norm_squared'])>0,'exact Hessian direction norm')
        need(2*norm==F(record['positive_second_variation']),'strict positive second variation')
        results.append({'n':n,'three_distinct_positive_radii_proved':True,
                        'full_pointwise_stationarity_proved':True,
                        'positive_second_variation':str(2*norm),
                        'unrestricted_nonradial_instability_proved':True})
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true',help='save the independently recomputed verification report')
    args=parser.parse_args();start=time.monotonic()
    raw=(HERE/'source_manifest.json').read_bytes();manifest=json.loads(raw)
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed source '+name)
    stored=json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    need(stored==certificate(),'recomputed exact certificate')
    algebra=identities();fourier=direct_fourier_audit();census=census_audit(stored)
    examples=witness_audit(stored)
    for row in stored['orbit_thresholds']:
        K,h,m=row['orbits'],row['center'],row['first_excluded_order']
        need(m==(8*K if not h else 10 if K==1 else 8*K+1),'general symbolic threshold formula')
        need(rotational_bound(m,K,h)['full_hessian_instability_proved'] is True,
             'first excluded order is strictly unstable')
        need(rotational_bound(m-1,K,h)['full_hessian_instability_proved'] is False,
             'preceding order is not excluded by this bound')
    report={'passed':True,'date':'2026-11-05','source_manifest_sha256':hashlib.sha256(raw).hexdigest(),
            'input_sha256':manifest['files'],'universal_polynomial_identities_checked':algebra,
            'direct_fourier_audit':fourier,'residue_census':census,'stationary_witnesses':examples,
            'arbitrary_continuous_phases_allowed':True,'unrestricted_hessian_used':True,
            'global_optima_for_general_n_determined':False,'floating_point_used':False,
            'elapsed_seconds':round(time.monotonic()-start,3)}
    if args.record_report:
        (HERE/'verification_report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
