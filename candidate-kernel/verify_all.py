"""Exact PSD-kernel identities and rational rejection of two sampled SDP fits.

This acceptor accepts identities and *rejections*, never a global sharp bound.
Only Python's standard library is required.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import hashlib,importlib.util,json

HERE=Path(__file__).resolve().parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def bilinear(a,Q,b):return sum((a[i]*Q[i][j]*b[j] for i in range(len(a)) for j in range(len(b))),F(0))
def cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cpower(z,k):
    out=(F(1),F(0))
    for _ in range(k):out=cmul(out,z)
    return out
def norm(z):return z[0]**2+z[1]**2


def hex_fraction(x):
    # A hexadecimal binary64 literal specifies an exact dyadic rational.
    return F.from_float(float.fromhex(x))


def decode(x):
    return [decode(q) for q in x] if isinstance(x,list) else hex_fraction(x)


def features(u,k,data):
    s=data['feature_order']
    if not k:
        nodes=decode(data['radial_nodes']);means=decode(data['radial_moment_centers'])
        return [u**j-means[j-1] for j in range(1,s)]+[10*prod(u-a for a in nodes)**2]
    E=decode(data['transforms'][k-1]);raw=[(u-1)**j for j in range(s)]
    return [sum((raw[i]*E[i][j] for i in range(s)),F(0)) for j in range(len(E[0]))]


def encoded_kernel(u,v,t,data):
    modes=data['harmonics_through'];H=[F(1),t]
    for k in range(2,modes+1):H.append(2*t*H[-1]-u*v*H[-2])
    Qs=decode(data['Q_blocks'])
    return sum((H[k]*bilinear(features(u,k,data),Q,features(v,k,data)) for k,Q in enumerate(Qs)),F(0))


def log_lower(d,terms=120):
    need(d>1,'positive atanh-series terms')
    a=(d-1)/(d+1)
    return 2*sum((a**(2*j+1)/(2*j+1) for j in range(terms)),F(0))


def reject_sampled_fit(n):
    data=json.loads((HERE/f'data/n{n}_exact_interpretation.json').read_bytes())
    discovery=json.loads((HERE/f'data/n{n}_discovery.json').read_bytes())
    need(data['n']==n and data['mass']==n-1,'fixed target normalization')
    need(decode(data['Q_blocks'])==[[[F.from_float(q) for q in row] for row in Q] for Q in discovery['Q_blocks']],
         'encoded coefficient blocks reproduce frozen discovery exactly')
    need(decode(data['transforms'])==[[[F.from_float(q) for q in row] for row in E] for E in discovery['nullspace_transforms']],
         'encoded nullspace transforms reproduce frozen discovery exactly')
    need(hex_fraction(data['B'])==F.from_float(discovery['candidate_log_delta']+discovery['candidate_gap']),
         'encoded B reproduces the frozen floating fitted target')
    u,v,t=map(F,data['test_pair']);d=u+v-2*t
    need(0<=u and 0<=v and t*t<u*v and d>F(4,n),'pair inside physical and strict separation domain')
    if n==10:
        q=9-F(9,8)*(u+v)-t/4;delta=5*q-28
        b=(81*u+v+18*t)*(u+81*v+18*t)/4096
        need(u<F(15403,10000) and v<F(15403,10000) and q>F(142,25),'inside new strict linear necessary domain')
        if 0<delta<2:
            need((delta+2)**2*d*b**8*q**28 >= (2-delta)**2*F(97,80)*9**9*7**21,
                 'pair passes nonlinear necessary cut too')
        domain='physical, strict separation, improved radial/complement caps and nonlinear Schwarz cut'
    else:
        need(n==12 and max(u,v)<F(2*(n-1)**2,n*n),'inside strict subset radial bound')
        domain='physical, strict separation and strict subset radial cap'
    B=hex_fraction(data['B']);N=n*(n-1)//2;b0=(B-N)/n
    def p(a):return (F(n,2)*a+b0-encoded_kernel(a,a,a,data)/2)/(n-1)
    h=p(u)+p(v)-encoded_kernel(u,v,t,data)
    lower=log_lower(d)
    # Coarse rational gaps are enough; no floating logarithm is used.
    gap=F(1,500) if n==10 else F(1,2000)
    need(lower-h>gap,'exact point defeats claimed continuous majorant')
    return {'n':n,'test_pair':[str(u),str(v),str(t)],'domain_checks':domain,
            'log_distance_minus_encoded_h_strict_lower':str(gap),
            'atanh_positive_terms':120,'encoded_fit_is_a_global_majorant':False,
            'candidate_global_optimality_refuted':False}


def fixture_test():
    # An exact centered D2 configuration, mass 9, with distinct layer radii.
    layers=[[(F(3,10),F(0)),(-F(3,10),F(0))],
            [(sx*F(9,10),sy*F(9,20)) for sx in (-1,1) for sy in (-1,1)],
            [(sx*F(3,10),sy*F(21,20)) for sx in (-1,1) for sy in (-1,1)]]
    nodes=[norm(layer[0]) for layer in layers];counts=[2,4,4]
    need(nodes==[F(9,100),F(81,80),F(477,400)],'exact fixture nodes')
    T2=sum(count*a*a for count,a in zip(counts,nodes))
    c=[sum((cpower(z,2)[0] for z in layer),F(0)) for layer in layers]
    need(c==[F(9,50),F(243,100),-F(81,20)],'nonzero exact layer second moments')
    Q2=[[dot(c,c)*(i==j)-c[i]*c[j] for j in range(3)] for i in range(3)]
    def radial(u,k):
        if k==0:return [u-F(9,10),u*u-T2/10,prod(u-a for a in nodes)**2]
        return [prod((u-nodes[j])/(nodes[r]-nodes[j]) for j in range(3) if j!=r) for r in range(3)]
    def K(z,w):
        u,v=norm(z),norm(w);t=dot(z,w);H=[F(1),t]
        for k in (2,3):H.append(2*t*H[-1]-u*v*H[-2])
        out=F(0)
        for k in range(4):
            a,b=radial(u,k),radial(v,k)
            out+=H[k]*(bilinear(a,Q2,b) if k==2 else dot(a,b))
        return out
    def test(points,should_vanish):
        need(sum(norm(z) for z in points)==9,'fixture mass')
        need(sum(z[0] for z in points)==sum(z[1] for z in points)==0,'fixture centered')
        total=sum((K(z,w) for z in points for w in points),F(0))
        penalty=F(0)
        for k in range(4):
            M=[(sum((cpower(z,k)[0]*radial(norm(z),k)[j] for z in points),F(0)),
                sum((cpower(z,k)[1]*radial(norm(z),k)[j] for z in points),F(0))) for j in range(3)]
            if k==2:
                # Exact Cauchy--Schwarz identity gives a PSD rank-two block.
                direct=bilinear([a for a,b in M],Q2,[a for a,b in M])+bilinear([b for a,b in M],Q2,[b for a,b in M])
                squares=sum((norm((c[i]*M[j][0]-c[j]*M[i][0],c[i]*M[j][1]-c[j]*M[i][1]))
                             for i in range(3) for j in range(i+1,3)),F(0))
                need(direct==squares and squares>=0,'rank-two PSD identity')
                penalty+=squares
            else:penalty+=sum((norm(m) for m in M),F(0))
        need(total==penalty and total>=0,'ordered pair sum equals sum of squared moments')
        need((total==0)==should_vanish,'candidate vanishing versus rotated-layer detection')
        B=F(42);b0=(B-45)/10
        def p(z):return (5*norm(z)+b0-K(z,z)/2)/9
        hsum=sum((p(z)+p(w)-K(z,w) for i,z in enumerate(points) for w in points[:i]),F(0))
        need(hsum==B-total/2,'diagonal-corrected two-point sum identity')
        return str(total)
    points=sum(layers,[]);base=test(points,True)
    turned=layers[0]+[cmul(z,(F(3,5),F(4,5))) for z in layers[1]]+layers[2]
    perturbed=test(turned,False)
    return {'n':10,'mass':9,'rational_fixture_candidate_kernel_sum':base,
            'one_layer_rotated_kernel_sum':perturbed,
            'positive_semidefinite_block_identity_checked':True,
            'two_point_summation_identity_checked':True,
            'fixture_is_the_actual_candidate_claimed':False}


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed kernel input '+name)
    spec=importlib.util.spec_from_file_location('kernel_exact_gram_check',HERE/'verify_psd.py')
    psd=importlib.util.module_from_spec(spec);spec.loader.exec_module(psd)
    gram=psd.verify()
    need(gram['passed'] is True,'exact Gram block acceptance')
    repaired_gap=F(1,500)-F(231,10**6)
    need(repaired_gap>0,'original failure survives exact PSD repair')
    report={'passed':True,'arithmetic':'exact rational arithmetic; positive logarithm series lower bounds',
            'rational_PSD_fixture':fixture_test(),
            'rejected_sampled_fits':[reject_sampled_fit(n) for n in (10,12)],
            'exact_Gram_block_checks':gram,
            'repaired_n10_fit_still_fails_same_pair_strict_gap_lower':str(repaired_gap),
            'continuous_global_majorant_accepted':False,'global_sharp_bound_proved':False}
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
