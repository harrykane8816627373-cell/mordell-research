from fractions import Fraction as Q
from math import comb
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('candidate', ROOT / 'verify_mordell.py')
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)

# Small formal polynomial algebra over Q. Monomials are tuples of variable names.
def add(*polys):
    out = {}
    for poly in polys:
        for key, value in poly.items():
            out[key] = out.get(key, Q(0)) + value
    return {key: value for key, value in out.items() if value}

def mul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(sorted(ka + kb))
            out[key] = out.get(key, Q(0)) + va * vb
    return {key: value for key, value in out.items() if value}

def scale(a, c):
    return {key: value * c for key, value in a.items() if value * c}

def const(c):
    return {(): Q(c)} if c else {}

def var(name):
    return {(name,): Q(1)}

def powpoly(a, n):
    out = const(1)
    for _ in range(n):
        out = mul(out, a)
    return out

def check(actual, expected, name):
    if actual != expected:
        raise ValueError((name, add(actual, scale(expected, -1))))

T, U, V, q, p, w2, m, mb = map(var, ['T','U','V','q','p','w2','m','mb'])
p2, pb2, p3, pb3, p4, pb4 = map(var, ['p2','pb2','p3','pb3','p4','pb4'])
moments = {(0,0):const(5), (1,0):{}, (0,1):{}, (1,1):const(5),
           (2,0):p2, (0,2):pb2, (3,0):p3, (0,3):pb3,
           (4,0):p4, (0,4):pb4, (2,1):{}, (1,2):{},
           (3,1):scale(p2,Q(7,4)), (1,3):scale(pb2,Q(7,4)),
           (4,1):scale(p3,Q(3,2)), (1,4):scale(pb3,Q(3,2)),
           (2,2):T, (3,3):U, (4,4):V, (4,2):m, (2,4):mb,
           (3,2):var('M32'), (2,3):var('M23'),
           (4,3):var('M43'), (3,4):var('M34')}

def reduce_products(poly):
    out = {}
    for key, value in poly.items():
        key = list(key)
        for a, b, replacement in [('p2','pb2','q'),('p3','pb3','p'),('p4','pb4','w2')]:
            while a in key and b in key:
                key.remove(a); key.remove(b); key.append(replacement)
        key = tuple(sorted(key))
        out[key] = out.get(key,Q(0)) + value
    return {key:value for key,value in out.items() if value}

distance_moments = {}
for k in range(1,5):
    terms = []
    for a in range(k+1):
        for b in range(k+1):
            terms.append(scale(mul(moments[k-a,k-b],moments[a,b]),
                               Q((-1)**(a+b)*comb(k,a)*comb(k,b),2)))
    distance_moments[k] = reduce_products(add(*terms))
expected = {
    1:const(25),
    2:add(scale(T,5),const(50),q),
    3:add(scale(U,5),scale(T,45),scale(q,Q(21,2)),scale(p,-1)),
    4:add(scale(V,5),scale(U,80),scale(mul(T,T),18),scale(q,49),scale(p,-12),w2,
          scale(mul(m,pb2),6),scale(mul(mb,p2),6))}
for k in range(1,5):
    check(distance_moments[k], expected[k], f'distance moment {k}')

# Derive the A/B coefficients from the unsimplified formulas in (23)-(24).
C, K = candidate.C, candidate.K
x = var('x')
x2, x3, x4 = (powpoly(x,n) for n in [2,3,4])

def derive(s, t, u, v):
    radial = add(scale(v,5),scale(u,30),scale(mul(t,t),18),scale(t,-275),const(750))
    f = add(scale(s,-Q(1,2)+Q(15,2)*C),scale(add(u,const(-5)),5*C),scale(radial,-K))
    d0 = add(scale(v,-60),scale(u,144),scale(t,-81))
    return f,d0

sA=scale(x2,20); tA=add(const(5),sA)
uA=add(const(5),scale(x2,60),scale(x3,-60))
vA=add(const(5),scale(x2,120),scale(x3,-240),scale(x4,260))
fA,dA=derive(sA,tA,uA,vA)
check(fA,add(scale(x2,-10+450*C-500*K),scale(x3,-300*C+3000*K),scale(x4,-8500*K)),'branch A f')
check(dA,add(const(15),scale(x2,-180),scale(x3,5760),scale(x4,-15600)),'branch A D')
sB=x; tB=add(const(5),x)
uB=scale(powpoly(tB,2),Q(1,5)); vB=scale(powpoly(tB,3),Q(1,25))
fB,dB=derive(sB,tB,uB,vB)
check(dB,scale(mul(tB,add(scale(powpoly(tB,2),-12),scale(tB,144),const(-405))),Q(1,5)), 'branch B D')

# Evaluate independently derived polynomials at rational points and compare exact
# point intervals with the candidate, including its outward square-root interval.
def evalpoly(poly, xv):
    return sum((value*xv**len(key) for key,value in poly.items()),Q(0))

point_checks=0
for branch,fp,dp,left,right in [('A',fA,dA,Q(1,20),Q(1,4)),('B',fB,dB,Q(5,4),Q(5,3))]:
    for i in range(41):
        xv=left+(right-left)*i/40
        sv=20*xv*xv if branch=='A' else xv
        for j in range(41):
            qv=(25-15*sv)*j/40
            rad=qv*(20-12*sv-Q(4,5)*qv)*(evalpoly(dp,xv)+5*qv)
            cut=max(Q(0),5+sv-(Q(4,7)+Q(9,8)/sv)*qv)
            own=(candidate.I(evalpoly(fp,xv)-candidate.AQ*qv+candidate.BQ*qv*qv-K*cut*cut)
                 + candidate.DQ*candidate.I(rad).sqrt())
            actual=candidate.upper_expression(branch,candidate.I(xv),candidate.I(qv))
            if (own.lo,own.hi)!=(actual.lo,actual.hi):
                raise ValueError((branch,xv,qv,'encoding mismatch'))
            point_checks+=1

result={'status':'PASSED','distance_moment_identities':[1,2,3,4],
        'scalar_polynomial_derivations':['A f','A D','B f','B D'],
        'exact_point_encoding_checks':point_checks,
        'note':'These checks audit algebra and encoding; they do not replace the analytic proof or interval coverage.'}
with open(ROOT / 'algebra_identity_result.json','w',encoding='utf-8') as out:
    json.dump(result,out,indent=2)
print(json.dumps(result,indent=2))
