"""Exact +/- ring-phase optimization by 2-adic population classes.

All arithmetic is integer/Fraction. This solves the finite sign family
at prescribed radii and populations, not unrestricted planar Mn.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod


def v2(m):
    if type(m) is not int or m < 1:
        raise ValueError('positive integer population required')
    return (m & -m).bit_length()-1


def polar_data(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n >= 2 required')
    k=0
    while 3*(k+1)*(k+2)<=n:
        k+=1
    n0=3*k*(k+1);r=n-n0
    rings=tuple([(6*j,6*j*j) for j in range(1,k+1)]+([(r,n+n0)] if r else []))
    if sum(m for m,b in rings)!=n or sum(m*b for m,b in rings)!=n*n:
        raise ArithmeticError('polar normalization failed')
    return n0,r,rings


def ring_product(n, rings, signs, center_singleton=False):
    """Exact squared-distance product at squared radii b_i/n."""
    if type(n) is not int or n<2 or sum(m for m,b in rings)!=n:
        raise ValueError('integer n >= 2 must equal the total population')
    if any(type(m) is not int or m<1 or type(b) is not int or b<=0 for m,b in rings):
        raise ValueError('positive integer populations and radius numerators required')
    if len(rings)!=len(signs) or any(s not in (-1,1) for s in signs):
        raise ValueError('one +/-1 sign per ring required')
    N=n*(n-1)//2
    numerator=prod(m**m*b**(m*(m-1)//2) for m,b in rings)
    for i,j in combinations(range(len(rings)),2):
        mi,a=rings[i];mj,b=rings[j]
        g=gcd(mi,mj);L=mi*mj//g
        if L%2:
            raise ValueError('this rational sign formulation requires even cross L')
        cross=abs(signs[j]**(mi//g)*b**(L//2)-signs[i]**(mj//g)*a**(L//2))
        if not cross:
            raise ValueError('collision')
        numerator*=cross**(2*g)
    value=F(numerator,n**N)
    singletons=[b for m,b in rings if m==1]
    if center_singleton:
        if sum(m*b for m,b in rings)!=n*n:
            raise ValueError('centroid correction requires initial total moment n')
        if len(singletons)>1:
            raise ValueError('centroid correction here permits at most one singleton ring')
        if singletons:
            value*=F(n**3,n**3-singletons[0])**N
    return value


def phase_graph(rings):
    """Build rational unary factors and equal-v2 edge factors."""
    rings=tuple(rings)
    if not rings or any(type(b) is not int or b<=0 for m,b in rings):
        raise ValueError('nonempty positive integer radius numerators required')
    if any(a>=b for (mi,a),(mj,b) in zip(rings,rings[1:])):
        raise ValueError('strictly increasing squared radii required')
    exponents=[v2(m) for m,b in rings]
    groups={e:tuple(i for i,a in enumerate(exponents) if a==e) for e in sorted(set(exponents))}
    unary=[F(1) for _ in rings]
    edge={}
    for i,j in combinations(range(len(rings)),2):
        mi,a=rings[i];mj,b=rings[j];g=gcd(mi,mj);L=mi*mj//g
        if L%2:
            raise ValueError('even cross least-common-multiples required')
        A,B=a**(L//2),b**(L//2)
        factor=F(B+A,B-A)**(2*g)
        if exponents[i]==exponents[j]:
            edge[i,j]=factor
        else:
            higher=i if exponents[i]>exponents[j] else j
            unary[higher]*=factor
    return groups,tuple(unary),edge


def graph_gain(signs, unary, edges):
    gain=prod(unary[i] for i,s in enumerate(signs) if s==-1)
    gain*=prod(f for (i,j),f in edges.items() if signs[i]!=signs[j])
    return gain


def optimize_signs(rings, max_states=1000000):
    """Globally optimize the entire finite sign family by independent blocks.

    The smallest v2 class has no unary factors, so fix its first sign +1
    using complement symmetry. No incomplete search is returned on a
    budget refusal.
    """
    groups,unary,edges=phase_graph(rings)
    first=min(groups)
    planned=sum(1<<(len(ids)-int(e==first)) for e,ids in groups.items())
    if max_states is not None and planned>max_states:
        raise ValueError('exact phase blocks need %s states, above budget %s'%(planned,max_states))
    signs=[1]*len(rings);blocks=[]
    for e,ids in groups.items():
        free=ids[1:] if e==first else ids
        best=F(0);best_signs=None
        for mask in range(1<<len(free)):
            candidate={i:1 for i in ids}
            for bit,i in enumerate(free):
                candidate[i]=-1 if mask>>bit&1 else 1
            gain=prod(unary[i] for i in ids if candidate[i]==-1)
            gain*=prod(f for (i,j),f in edges.items() if i in candidate and j in candidate and candidate[i]!=candidate[j])
            if gain>best:
                best=gain;best_signs=candidate
        for i,s in best_signs.items():signs[i]=s
        blocks.append({'v2':e,'indices':ids,'states':1<<len(free),'gain':best})
    gain=prod(block['gain'] for block in blocks)
    if gain!=graph_gain(signs,unary,edges):
        raise ArithmeticError('block product mismatch')
    return {'signs':tuple(signs),'gain':gain,'states':planned,'blocks':blocks}


def optimized_polar(n,max_states=1000000):
    n0,r,rings=polar_data(n)
    optimized=optimize_signs(rings,max_states)
    value=ring_product(n,rings,optimized['signs'],True)
    minus=(-1,)*len(rings);flip=(1,)+minus[1:]
    previous=max(ring_product(n,rings,s,True) for s in (minus,flip))
    raw=ring_product(n,rings,(1,)*len(rings))
    centered_base=ring_product(n,rings,(1,)*len(rings),True)
    if value!=centered_base*optimized['gain'] or value<previous:
        raise ArithmeticError('exact family maximum consistency failed')
    return {'n':n,'n0':n0,'remainder':r,'rings':rings,'signs':optimized['signs'],
            'states':optimized['states'],'class_sizes':[len(b['indices']) for b in optimized['blocks']],
            'value':value,'previous_two_sign_value':previous,'original_polar_value':raw,
            'scope':'Exact maximum only over +/- phases of prescribed polar rings, with the same singleton centering.'}


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('n',type=int)
    args=p.parse_args();result=optimized_polar(args.n)
    print(json.dumps(result,default=str,ensure_ascii=False,indent=2))
