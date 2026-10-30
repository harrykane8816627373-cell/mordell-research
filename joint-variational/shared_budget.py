"""Shifted Gram duals share ONE Dirichlet-plus-pair energy budget.

Measured moments and any orthogonal-energy input are caller hypotheses.
The radial adapter supplies exact Gram/cross terms with support checks.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys
from affine_pair import pi_bounds,sqrt_bounds,require,ROOT
sys.path.insert(0,str(ROOT/'radial-angular'))
from mixed_energy import solve


def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))


def checked_gram(H,m,c,Hf):
    H=[list(map(F,row)) for row in H];m=list(map(F,m));c=list(map(F,c));Hf=F(Hf)
    d=len(H)
    if not d or len(m)!=d or len(c)!=d or any(len(row)!=d for row in H):raise ValueError('matching nonempty Gram dimensions')
    if Hf<0 or any(H[i][j]!=H[j][i] for i in range(d) for j in range(d)):raise ValueError('symmetric Gram and nonnegative norm')
    hm,hc=solve(H,m),solve(H,c)
    J=dot(m,hm);p=dot(c,hm);h=Hf-dot(c,hc)
    if h<0:raise ValueError('negative residual norm: incompatible augmented Gram')
    return H,m,c,Hf,J,p,h


def radial_packet(n,a,T,measurements):
    a,T=F(a),F(T);m=list(map(F,measurements));degree=len(m)
    if type(n)is not int or n<2 or not 0<a<=1 or T<1 or not 1<=degree<=12:
        raise ValueError('n>=2, 0<a<=1, T>=1, radial degree 1..12')
    margin=T-F(n-1,n)
    if margin<0 or margin*margin<F(1,n):raise ValueError('T must cover all globally allowed smeared disks')
    if m[0]!=0:raise ValueError('matched second moment forces the first radial measurement to zero')
    S=T*T
    H=[[F(2*j*k,(j+k)*(j+k-1))*S**(j+k) for k in range(1,degree+1)] for j in range(1,degree+1)]
    c=[]
    for j in range(1,degree+1):
        c.append(2*j*sum(((-1)**l*comb(j-1,l)*S**(j-1-l)*
            (a**(2*l+2)/F(2*l+3)+(T**(2*l+1)-a**(2*l+1))/F(2*l+1)) for l in range(j)),F(0)))
    Hf=(1+1/(a*a))/2
    checked_gram(H,m,c,Hf)
    return {'n':n,'a':a,'T':T,'H':H,'m':m,'c':c,'Hf':Hf}


def geometry_constants(n,K,a):
    a=F(a)
    if type(n)is not int or n<2 or type(K)is not int or K<1 or 2*K>=n or not 0<a<=1:
        raise ValueError('n>=2, integer 1<=K<n/2, 0<a<=1')
    root=sqrt_bounds(n,10**12)[1]
    V=K*root+F(K,4)/a-F(n,4)*(1-a+a*a/3)+F(1,32)/a**3
    return n-2*K,V


def add_centered_holes(packet,radii):
    """Append overlapping central bubble tests, conditional on EMPTY smeared disks.

    This adapter requires an unaugmented radial-polynomial packet. It does
    not infer empty disks from radii, point separation, or moment data.
    """
    radii=list(map(F,radii));degree=len(packet['m']);a,T=packet['a'],packet['T'];S=T*T
    if not radii or len(radii)!=len(set(radii)) or any(not 0<b<=1 for b in radii):
        raise ValueError('distinct central energy-hole radii in (0,1]')
    expected=radial_packet(packet['n'],a,T,packet['m'])
    if packet!=expected:raise ValueError('unaugmented radial packet required')
    H=[row[:] for row in packet['H']];cross=[]
    for b in radii:
        cross.append([4*j*sum(((-1)**l*comb(j-1,l)*S**(j-1-l)*b**(2*l+4)/F(2*l+4)
                              for l in range(j)),F(0)) for j in range(1,degree+1)])
    for j in range(degree):H[j].extend(row[j] for row in cross)
    for i,b in enumerate(radii):H.append(cross[i]+[min(b,d)**4 for d in radii])
    c=packet['c']+[2*min(a,b)**3/(3*a)+2*max(F(0),b-a) for b in radii]
    m=packet['m']+[-F(packet['n'],2)*b**4 for b in radii]
    checked_gram(H,m,c,packet['Hf'])
    return dict(packet,H=H,m=m,c=c,centered_energy_hole_radii=radii)


def exact_certificate(packet,K,line,orthogonal_energy_lower=0,disjoint_pair_lower=0):
    H,m,c,Hf,J,p,h=checked_gram(packet['H'],packet['m'],packet['c'],packet['Hf'])
    A,C=F(orthogonal_energy_lower),F(disjoint_pair_lower)
    if min(A,C)<0:raise ValueError('nonnegative independent energy inputs')
    L,V=geometry_constants(packet['n'],K,packet['a'])
    alpha,beta=F(line['alpha']),F(line['beta'])
    if alpha<0 or beta<=0:raise ValueError('positive certified affine-pair parameters')
    pl,pu=pi_bounds();linear=beta*(p/2-2*V)
    lower=A+C+J+alpha*L+min(pl*linear,pu*linear)-pu*pu*beta*beta*h/16
    return {'n':packet['n'],'maximum_layers':K,'projection_energy':J,'transport_projection':p,
            'residual_norm_squared':h,'adjacent_edge_count_lower':L,'base_circumference_upper':V,
            'orthogonal_energy_lower':A,'disjoint_pair_lower':C,'shared_gap_lower':lower,
            'scope':'Conditional on actual measurements, certified affine pair line, valid full Gram/cross terms, and energy components orthogonal/disjoint as specified. This is one joint lower bound, not a sum of overlapping old certificates.'}


def interval_certificate(packet,K,line,intervals,coefficients=None,orthogonal_energy_lower=0):
    H,m,c,Hf,J,p,h=checked_gram(packet['H'],packet['m'],packet['c'],packet['Hf'])
    bounds=[tuple(map(F,row)) for row in intervals]
    if len(bounds)!=len(H) or any(len(row)!=2 or row[0]>row[1] for row in bounds):raise ValueError('ordered matching measurement intervals')
    alpha,beta=F(line['alpha']),F(line['beta']);A=F(orthogonal_energy_lower)
    if alpha<0 or beta<=0 or A<0:raise ValueError('positive pair parameters and nonnegative orthogonal input')
    pl,pu=pi_bounds();pm=(pl+pu)/2
    if coefficients is None:
        mid=[(lo+hi)/2+pm*beta*x/4 for (lo,hi),x in zip(bounds,c)]
        v=solve(H,mid)
    else:v=list(map(F,coefficients))
    if len(v)!=len(H):raise ValueError('matching dual coefficients')
    L,V=geometry_constants(packet['n'],K,packet['a'])
    pairing=sum((min(x*lo,x*hi) for x,(lo,hi) in zip(v,bounds)),F(0))
    norm=sum((v[i]*H[i][j]*v[j] for i in range(len(H)) for j in range(len(H))),F(0))
    linear=beta*(dot(v,c)/2-2*V)
    lower=A+alpha*L+2*pairing-norm+min(pl*linear,pu*linear)-pu*pu*beta*beta*Hf/16
    return {'shared_gap_lower':lower,'fixed_dual_vector':v,'measurement_intervals':bounds,
            'scope':'One fixed shifted-Gram dual covers every actual measurement in the entire supplied box; midpoint only selects coefficients.'}


def credit_ledger(n,claims,reserved_edges=()):
    if type(n)is not int or n<2:raise ValueError('integer n>=2')
    def edge(raw):
        if len(raw)!=2 or any(type(x)is not int for x in raw):raise ValueError('integer edge endpoints')
        i,j=sorted(raw)
        if not 0<=i<j<n:raise ValueError('distinct endpoints inside n')
        return i,j
    used={edge(e):F(1) for e in reserved_edges};cost=F(0)
    for claim in claims:
        weight,lower=F(claim['weight']),F(claim['per_pair_lower'])
        if not 0<=weight<=1 or lower<0:raise ValueError('nonnegative fractional credit <=1')
        edges=[edge(e) for e in claim['edges']]
        if len(edges)!=len(set(edges)):raise ValueError('duplicate edge inside a claim')
        for e in edges:
            used[e]=used.get(e,F(0))+weight
            if used[e]>1:raise ValueError('pair credited more than once or overlaps a reserved layer edge')
        cost+=weight*lower*len(edges)
    return {'disjoint_pair_lower':cost,'edge_credit_totals':[(e,w) for e,w in sorted(used.items())],
            'scope':'Each claimed per-pair lower and any geometric edge classification must be proved separately; this ledger only validates resource allocation.'}
