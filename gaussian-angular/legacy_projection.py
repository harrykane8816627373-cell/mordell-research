"""Full cross terms with old holomorphic/first mixed angular tests.

The interval dual remains valid without treating midpoint Gram entries as exact.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
from gaussian_energy import (add,sub,scale,mul,square,outward,require,cpow,norm,
    exp_negative,bessel_sum,center_power,measurement,probe_gram)

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'radial-angular'))
from mixed_energy import solve


def old_mode(points,k,S):
    points=[tuple(map(F,p)) for p in points]
    n=len(points);S=F(S)
    if type(k)is not int or not 1<=k<=16 or n<2 or S<1 or any(len(p)!=2 for p in points):raise ValueError('k=1..16, n>=2 planar points, S>=1')
    qmax=max(norm(p) for p in points);margin=S-qmax-F(1,n)
    if margin<0 or margin*margin<4*qmax/n:raise ValueError('support fails to cover actual smear')
    p=[sum((cpow(z,k)[j] for z in points),F(0)) for j in range(2)]
    Q=[sum((norm(z)*cpow(z,k)[j] for z in points),F(0)) for j in range(2)]
    bubble=[(S-F(k+1,2*n))*p[j]-Q[j] for j in range(2)]
    diag=[k*S**k,k*S**k,S**(k+2)/F(k+2),S**(k+2)/F(k+2)]
    values=p+bubble
    return {'mode':k,'support_squared_q':S,'diagonal':diag,'measurements':values,
            'old_projection':sum((v*v/h for v,h in zip(values,diag)),F(0))}


def cross_mode(probe,k,S):
    S=F(S);t=probe['t'];rho2=probe['center_squared_radius']
    x=rho2/(4*t);y=S/(4*t);z=S*rho2/(16*t*t)
    power=center_power(probe,k)
    pure=scale(mul(exp_negative(x+y),bessel_sum(z,k)),2*k*S**k/(4*t)**k)
    inside=scale(mul(exp_negative(x),bessel_sum(z,k,y)),2*(k+1)*S**(k+1)/(4*t)**k)
    bubble=sub(inside,scale(pure,S/k))
    return [scale(pure,power[0]),scale(pure,power[1]),
            scale(bubble,power[0]),scale(bubble,power[1])]


def assembled_packet(old,probes,measurements):
    d=4;count=len(probes)
    if count<1 or len(measurements)!=count:raise ValueError('matching nonempty Gaussian data')
    H=[[(F(0),F(0)) for _ in range(d+count)] for _ in range(d+count)]
    for i,value in enumerate(old['diagonal']):H[i][i]=(value,value)
    for j,probe in enumerate(probes):
        for i,c in enumerate(cross_mode(probe,old['mode'],old['support_squared_q'])):
            H[i][d+j]=H[d+j][i]=c
        for l,other in enumerate(probes):H[d+j][d+l]=probe_gram(probe,other)
    m=[(v,v) for v in old['measurements']]+list(measurements)
    return H,m


def choose_dual(H,m,mesh=10**9):
    mid=[[sum(v)/2 for v in row] for row in H];rhs=[sum(v)/2 for v in m]
    v=solve(mid,rhs)
    out=[]
    for x in v:
        y=x*mesh;out.append(F((y.numerator*2+y.denominator)//(2*y.denominator),mesh))
    return out


def dual_lower(H,m,v):
    d=len(m);v=list(map(F,v))
    if len(H)!=d or len(v)!=d or any(len(row)!=d for row in H):raise ValueError('matching dual dimensions')
    if any(lo>hi for lo,hi in m) or any(lo>hi for row in H for lo,hi in row):raise ValueError('ordered intervals')
    pairing=sum((min(a*lo,a*hi) for a,(lo,hi) in zip(v,m)),F(0))
    penalty=sum((max(v[i]*v[j]*H[i][j][0],v[i]*v[j]*H[i][j][1]) for i in range(d) for j in range(d)),F(0))
    return 2*pairing-penalty


def packet_certificate(old,probes,values):
    H,m=assembled_packet(old,probes,values);v=choose_dual(H,m)
    lower=dual_lower(H,m,v)
    return {'energy_lower':lower,'old_projection':old['old_projection'],
            'improvement_lower':lower-old['old_projection'],'fixed_dual_vector':v,
            'Gram_intervals':H,'measurement_intervals':m,
            'scope':'Full old-mode/Gaussian Gram, actual disk means or independently proved intervals. No stationarity or optimality is inferred.'}
