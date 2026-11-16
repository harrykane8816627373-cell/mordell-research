"""Exact certificate for the sum of the two conditioned inner rows."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INNER = (F(37,500),F(43,250))
OUTER = (F(169,200),F(164,125))
DOT = (F(-27,200),F(-1,10))
BUDGET = F(2101,10**6)
D = {(1,0,0):F(1),(0,1,0):F(1),(0,0,1):F(-2)}


def load(path):
    spec = importlib.util.spec_from_file_location(path.replace('/', '_'), ROOT/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def geometry(box):
    (u0,u),(v0,v),(t0,t1)=box
    return (u+v-2*t0<F(11,20) and u0+v0-2*t1>F(41,100)
            and 10000*t1*t1>5929*u*v)


def add_d(P,U):
    for x,p in U.items():
        for y,q in D.items():
            k=tuple(a+b for a,b in zip(x,y))
            P[k]=P.get(k,F(0))+p*q


def outer_poly(P,L):
    A,B,a,b,c=L
    Q=P.copy()
    add_d(Q,{(0,0,0):-A,(1,0,0):-B,(0,1,0):-a,(0,0,1):-b})
    Q[1,0,0]=Q.get((1,0,0),F(0))-c
    Q[0,0,1]=Q.get((0,0,1),F(0))+c
    return Q


def pair_poly(L,R):
    A,B,a,b,c=L
    E,G,e,f,g=R
    correction={}
    add_d(correction,{(0,0,0):8*(A+E)+9*(a+e)-BUDGET,
                      (1,0,0):8*B-a-e-b+5*c,(0,1,0):8*G-a-e-f+5*g,(0,0,1):-b-f})
    for k,x in (((1,0,0),-c),((0,1,0),-g),((0,0,1),c+g)):
        correction[k]=correction.get(k,F(0))+x
    return correction


def value(poly, point):
    return sum(c*point[0]**i*point[1]**j*point[2]**k for (i,j,k),c in poly.items())


def consistent(outer, polynomial, L, R):
    u,v,t=F(1,9),F(2,15),F(-2,15)
    d=u+v-2*t
    A,B,a,b,c=L
    E,G,e,f,g=R
    a_row=8*(A+B*u)+a*(9-u-v)+b*(-u-t)+c*(5*u-(u-t)/d)
    b_row=8*(E+G*v)+e*(9-u-v)+f*(-v-t)+g*(5*v-(v-t)/d)
    need(value(pair_poly(L,R),(u,v,t))==d*(a_row+b_row-BUDGET),'conditioned pair correction identity')
    u,v,t=F(1,7),F(111,100),F(1,4)
    d=u+v-2*t
    expected=value(polynomial,(u,v,t))-d*(A+B*u+a*v+b*t)-c*(u-t)
    need(value(outer,(u,v,t))==expected,'outer support correction identity')


def quadratic_lanes(Q, shape, box, pack, upper):
    """Return lower Bernstein padding for a quadratic, with no bias."""
    q={}
    for k,v in Q.items():
        need(sum(k)<=2,'quadratic correction')
        for ell in product(*(range(i+1) for i in k)):
            value=v
            for j,(a,b) in enumerate(box):
                m,n=k[j],ell[j]
                value*=((2 if m==2 and n==1 else 1)*a**(m-n)*(b-a)**n)
            q[ell]=q.get(ell,F(0))+value
    result=0;amplitude=0
    degrees=tuple(d-1 for d in shape)
    for k,v in q.items():
        cells=[j for j,i in enumerate(k) if i]
        if not cells:
            mask,bound=pack.ones,1
        elif len(cells)==1:
            j=cells[0];n=degrees[j]
            mask,bound=(pack.index[j],n) if k[j]==1 else (pack.choose2[j],n*(n-1))
        else:
            j,l=cells
            mask,bound=pack.cross[(j,l)],degrees[j]*degrees[l]
        coeff=upper.S*v//bound
        result+=coeff*mask
        amplitude+=abs(coeff)*bound
    # The base has every lane between 2^59 and 2^61 at its root. Midpoint
    # splitting preserves that interval. This padding cannot under/overflow it
    # into adjacent lanes or violate the comparator's bound 2^62.
    need(amplitude<1<<58,'safe signed correction in every lane')
    return result


def count(action,stats):
    stats[action]=stats.get(action,0)+1


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--record-report',action='store_true');args=ap.parse_args()
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for directory,key in ((HERE,'files'),(ROOT,'dependencies')):
        for file,digest in manifest[key].items():
            need(sha256((directory/file).read_bytes()).hexdigest()==digest,'fixed '+file)
    banded=load('ten-point-banded-kernel/verify_all.py')
    row=load('ten-point-ward-radii/verify_all.py')
    domain=load('ten-point-counted-radii/verify_all.py')
    refined=load('ten-point-banded-kernel-refined/verify_all.py')
    upper=banded.upper
    profile=json.loads((ROOT/'ten-point-refined-radii/verification_report.json').read_bytes())
    need(profile['passed'] and profile['applies_to_every_global_maximizer'] and profile['exact_inner_count']==2,
         'certified two inner points')
    need(tuple(map(F,profile['inner_squared_radius_open_interval']))==INNER
         and tuple(map(F,profile['outer_squared_radius_open_interval']))==OUTER,'fixed radial bands')
    midpoint=F(239,1000)**2
    need(4*F(profile['midpoint']['inner_midpoint_squared_strict_upper'])==midpoint,
         'fixed bound on the squared sum of the inner pair')
    report=json.loads((ROOT/'ten-point-banded-kernel-refined/verification_report.json').read_bytes())
    need(report['passed'] and F(report['c'])==refined.C==F(3,40000),'confirmed old outer-pair deficit')
    need(tuple(map(F,report['inner_pair']['inner_dot_product_open_interval']))==DOT
         and F(report['inner_pair']['angle_cosine_strict_upper'])==F(-2,3)
         and F(report['inner_pair']['inner_squared_distance_strict_upper'])==F(3,5),'fixed inherited inner geometry')
    data=json.loads((ROOT/'ten-point-banded-kernel-refined/kernel.json').read_bytes())
    need(upper.hx(data['B'])-upper.log.log_lower(upper.witness())+28*refined.C<BUDGET,'confirmed budget for 17 edges')
    need(F(12,25)**2>INNER[1]*OUTER[1],'outer dot box covers all feasible projections')
    A,force=upper.corrections(data)
    base,shape=upper.polynomial(data,A,force)
    weights=banded.parsed(data)
    O=banded.root_polynomial(base,shape,weights,0,1)
    I=banded.root_polynomial(base,shape,weights,0,0)
    pack=upper.Packed(shape,True);upper.old.selftest()
    certificate=json.loads((HERE/'cover.json').read_bytes())
    marks=[tuple(map(F,b)) for b in certificate['marks']]
    need(marks and marks[0][0]==INNER[0] and marks[-1][1]==INNER[1]
         and all(a<b for a,b in marks)
         and all(b==c for (_,b),(c,_) in zip(marks,marks[1:])), 'all inner mark radii covered exactly')
    support={};outside={}
    for i,j,co,tree in certificate['outer']:
        need((i,j) not in support and 0<=i<len(marks) and j>=0,'unique support key')
        L=tuple(map(F,co));support[i,j]=L
        hull=(marks[i],OUTER,(-F(12,25),F(12,25)))
        P=outer_poly(O,L);consistent(P,O,L,L)
        stack=[(hull,upper.initial(P,shape,hull),0)]
        for action in tree:
            need(stack and action in 'AP012','valid outer support node')
            box,B,depth=stack.pop()
            count(action,outside)
            if action=='A':need(pack.dominates(B,box,F(0)),'exact outer-neighbor support lower bound')
            elif action=='P':need(domain.prunable(row,upper,box),'exact outer-neighbor prune')
            else:stack.extend(banded.split(pack,B,box,depth,int(action)))
        need(not stack,'all outer-neighbor boxes resolved')
    exact_keys=[(i,j) for i in range(len(marks)) for j in range(i,len(marks))]
    need([(i,j) for i,j,_ in certificate['inner']]==exact_keys,'all unordered inner slab pairs covered')
    inside={}
    for i,j,tree in certificate['inner']:
        hull=(marks[i],marks[j],DOT)
        init=upper.initial(I,shape,hull)
        need(((init+pack.high-(1<<59)*pack.ones)&pack.high)==pack.high,'initial base lanes above 2^59')
        need(((init+pack.high-(1<<61)*pack.ones)&pack.high)==0,'initial base lanes below 2^61')
        stack=[(hull,init,0)]
        for action in tree:
            need(stack,'valid inner pair node')
            box,B,depth=stack.pop()
            if isinstance(action,list):
                x,y=action
                left,right=support[i,x],support[j,y]
                consistent(outer_poly(O,left),O,left,right)
                correction=quadratic_lanes(pair_poly(left,right),shape,box,pack,upper)
                need(pack.dominates(B+correction,box,F(0)),'two inner rows exceed the allowed joint budget')
                count('A',inside)
            elif action=='G':
                need(geometry(box),'verified stricter inner angle and distance')
                count(action,inside)
            elif action=='P':
                (u0,u),(v0,v),(t0,t1)=box
                need(u0>v or t1*t1>u*v or u0+v0-2*t1>F(3,5) or u0+v0+2*t0>=midpoint
                     or 9*t0*t0<4*u0*v0 or domain.prunable(row,upper,box),'valid exact inherited inner pair prune')
                count(action,inside)
            else:
                need(action in '012','valid inner split')
                stack.extend(banded.split(pack,B,box,depth,int(action)))
                count(action,inside)
        need(not stack,'all inner boxes resolved')
    output={'passed':True,'applies_to_every_global_maximizer':True,'normalization_mass':9,
            'inner_squared_distance_open_interval':['41/100','11/20'],'inner_angle_cosine_strict_upper':'-77/100',
            'outer_supports':len(support),'outer_support_counts':outside,'inner_pair_counts':inside,'unresolved':0}
    text=json.dumps(output,indent=2)+'\n'
    if args.record_report:(HERE/'verification_report.json').write_bytes(text.encode())
    print(text,end='')


if __name__=='__main__':main()
