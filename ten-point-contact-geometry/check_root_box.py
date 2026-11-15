"""Independent exact replay of the five-dimensional Jacobian root enclosure."""
import json
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def plus(*polys):
    out={}
    for p in polys:
        for e,c in p.items():
            out[e]=out.get(e,0)+c
    return {e:c for e,c in out.items() if c}


def scaled(p,n):
    return {e:c*n for e,c in p.items() if c*n}


def mul(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            key=tuple(i+j for i,j in zip(e,f))
            out[key]=out.get(key,0)+c*d
    return {e:c for e,c in out.items() if c}


def system():
    one={(0,)*5:1}
    a,b,c,x,y=({tuple(int(i==j) for i in range(5)):1} for j in range(5))
    A,B,C=[one,a],[scaled(one,2),scaled(x,2)],[scaled(one,2),scaled(y,2)]
    for k in range(2,5):
        A.append(mul(a,A[k-1]))
        for L,r,s in ((B,b,x),(C,c,y)):
            L.append(plus(scaled(mul(s,L[k-1]),2),scaled(mul(mul(r,r),L[k-2]),-1)))
    S=[plus(A[k],B[k],C[k]) for k in range(5)]
    T=[plus(mul(a,A[k]),mul(b,B[k]),mul(c,C[k])) for k in range(5)]
    return [plus(scaled(S[k],-(2*k+1)),scaled(T[k],-10),
                 *(scaled(mul(S[j],S[k-j]),2) for j in range(k+1))) for k in range(5)]


def rat(s):
    return F(int(s['num']),int(s['den']))


def interval(s):
    return rat(s['lower']),rat(s['upper'])


def check_box(s,names):
    need(tuple(s['variables'])==names,'box variable order')
    box=[interval(i) for i in s['intervals']]
    need(len(box)==len(names) and all(l<=u for l,u in box),'valid rational box')
    return box


def isum(items):
    pairs=list(items)
    return sum((a for a,b in pairs),F(0)),sum((b for a,b in pairs),F(0))


def imul(x,y):
    p=[a*b for a in x for b in y]
    return min(p),max(p)


def ipow(x,k):
    if not k:
        return F(1),F(1)
    a,b=x
    if k%2:
        return a**k,b**k
    return (F(0) if a<=0<=b else min(a**k,b**k)),max(a**k,b**k)


def ival(p,box):
    def term(e,c):
        x=(F(c),F(c))
        for f,y in zip(e,box):
            x=imul(x,ipow(y,f))
        return x
    return isum(term(e,c) for e,c in p.items())


def dif(p,k):
    return {tuple(v-int(i==k) for i,v in enumerate(e)):e[k]*c
            for e,c in p.items() if e[k]}


def check():
    raw=json.loads((HERE/'root_box_jacobian.json').read_text())
    need(raw['operation_id']=='polynomial.system.real_root_box.certify','operation ID')
    out=raw['output']
    names=tuple('abcxy')
    data=out['polynomial_map']
    need(tuple(data['input_variables'])==names,'system variable order')
    polys=system()
    need(len(data['output_polynomials'])==5,'five equations')
    for p,q in zip(polys,data['output_polynomials']):
        need(q['domain']=='QQ' and tuple(q['variables'])==names,'component variable order')
        terms=q['polynomial']['terms']
        got={tuple(t['exponents']):rat(t['coefficient']) for t in terms}
        need(len(got)==len(terms) and p==got,'independent sparse polynomial construction')
    box=check_box(out['box'],names)
    a,b,c,x,y=box
    need(0<a[0]<a[1]<b[0]<b[1]<c[0]<c[1]<F(17,10),'nondegenerate ordered radii')
    need(0<max(abs(i) for i in x)<b[0] and x[1]<0,'nondegenerate middle rectangle')
    need(0<max(abs(i) for i in y)<c[0] and y[0]>0,'nondegenerate outer rectangle')
    center=[(l+u)/2 for l,u in box]
    need(center==[F(n,2**20) for n in (118037,1084060,1216217,-797427,810046)],'dyadic center')
    need(all(u-l==F(1,2**19) for l,u in box),'dyadic width')
    cert=out['conclusion']
    need(cert['status']=='CERTIFIED_UNIQUE_NONSINGULAR','unique nonsingular output')
    ev=cert['evidence']
    need(tuple(ev['center']['variables'])==names and
         [rat(z) for z in ev['center']['values']]==center,'recorded center')
    value=[rat(z) for z in ev['value_at_center']]
    point=[(x,x) for x in center]
    need([ival(p,point) for p in polys]==[(x,x) for x in value],'exact point values')
    J=[[dif(p,k) for k in range(5)] for p in polys]
    def matrix(s):
        need(s['row_count']==s['column_count']==len(s['entries'])==5 and
             all(len(r)==5 for r in s['entries']),'5x5 matrix')
        return [list(map(rat,r)) for r in s['entries']]
    midJ=matrix(ev['jacobian_at_center'])
    need([[ival(p,point) for p in row] for row in J]==[[(x,x) for x in row] for row in midJ],
         'point derivative')
    inverse=matrix(ev['preconditioner'])
    for a0,b0 in ((midJ,inverse),(inverse,midJ)):
        need(all(sum(a0[i][k]*b0[k][j] for k in range(5))==int(i==j)
                 for i in range(5) for j in range(5)),'two-sided inverse')
    enclosure=[[interval(e) for e in row] for row in ev['jacobian_enclosure']['entries']]
    need([[ival(p,box) for p in row] for row in J]==enclosure,'interval derivative')
    rows=[]
    image=[]
    for i in range(5):
        row=[]
        for j in range(5):
            l,u=isum(imul((inverse[i][k],)*2,enclosure[k][j]) for k in range(5))
            row.append((int(i==j)-u,int(i==j)-l))
        rows.append(row)
        l,u=isum(imul(row[j],(box[j][0]-center[j],box[j][1]-center[j])) for j in range(5))
        q=center[i]-sum(inverse[i][j]*value[j] for j in range(5))
        image.append((q+l,q+u))
    need(image==check_box(ev['krawczyk_image'],names),'exact Krawczyk image')
    need(all(l<a0<=b0<u for (a0,b0),(l,u) in zip(image,box)),'strict interior')
    norm=max(sum(max(abs(l),abs(u)) for l,u in row) for row in rows)
    need(norm<F(21,10000),'strict contraction with numerical margin')

    # Equality recovery: exactly one integer population is consistent with mass 9.
    feasible=[]
    closest=F(100)
    for p in range(11):
        for q in range(11-p):
            counts=(p,q,10-p-q)
            l,u=isum(imul((F(m),)*2,r) for m,r in zip(counts,box))
            if l<=9<=u:
                feasible.append(counts)
            else:
                closest=min(closest,abs(l-9),abs(u-9))
    need(feasible==[(2,4,4)],'mass and support force 2/4/4 populations')
    need(closest>F(1,8),'population alternatives separated by at least 1/8 in mass')
    return box,norm,closest


if __name__=='__main__':
    _,norm,closest=check()
    print('Exact five-variable root box replay passed')
    print('Uniform ||I-C Df||_infinity <',float(norm))
    print('Mass 9 forces population (2, 4, 4); other populations have margin >',float(closest))
