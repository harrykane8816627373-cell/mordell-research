"""Rigorous upper Bernstein coefficients of lambda*(uv-t²), lambda>=0.

Subtracting this polynomial makes a stronger lower-bound test on the physical
domain uv-t²>=0. The choice of lambda may be exploratory; positivity of the
resulting Bernstein lower coefficients remains a rigorous sufficient test.
"""
import numpy as np
import math
from flint import arb,fmpq

def power_bounds(a,b,n,power):
    values=[]
    for k in range(n+1):
        value=sum((math.comb(power,j)*a**(power-j)*b**j*arb(fmpq(math.comb(k,j),math.comb(n,j)))for j in range(min(k,power)+1)),arb(0))
        values.append(value)
    return (np.array([np.nextafter(float(v.lower()),-np.inf)for v in values]),
            np.array([np.nextafter(float(v.upper()),np.inf)for v in values]))

def exact_q_coefficients(box,degrees):
    def values(axis,power):
        n=degrees[axis];a=arb(box[2*axis]);b=arb(box[2*axis+1])-a
        return [sum((math.comb(power,j)*a**(power-j)*b**j*arb(fmpq(math.comb(k,j),math.comb(n,j)))for j in range(min(k,power)+1)),arb(0))for k in range(n+1)]
    u,v,t2=values(0,1),values(1,1),values(2,2)
    u2,v2,t4=values(0,2),values(1,2),values(2,4)
    shape=tuple(n+1 for n in degrees);q=np.empty(shape,dtype=object);q2=np.empty(shape,dtype=object)
    for i in range(shape[0]):
        for j in range(shape[1]):
            uv=u[i]*v[j];uv2=u2[i]*v2[j]
            for k in range(shape[2]):
                q[i,j,k]=uv-t2[k]
                q2[i,j,k]=uv2-2*uv*t2[k]+t4[k]
    return q,q2

def correction_upper(box,degrees,lam,kappa=0.):
    assert np.isfinite(lam)and lam>=0
    assert np.isfinite(kappa)and kappa>=0
    ul=[];vl=[]
    for axis,out in [(0,ul),(1,vl)]:
        n=degrees[axis];a=arb(box[2*axis]);b=arb(box[2*axis+1])-a
        for k in range(n+1):out.append(np.nextafter(float((a+b*arb(fmpq(k,n))).upper()),np.inf))
    assert min(ul)>=0 and min(vl)>=0
    n=degrees[2];a=arb(box[4]);b=arb(box[5])-a
    t2=[]
    for k in range(n+1):
        value=a*a+2*a*b*arb(fmpq(k,n))+b*b*arb(fmpq(k*(k-1),n*(n-1)))
        t2.append(np.nextafter(float(value.lower()),-np.inf))
    uv=np.nextafter(np.array(ul)[:,None,None]*np.array(vl)[None,:,None],np.inf)
    q=np.nextafter(uv-np.array(t2)[None,None,:],np.inf)
    result=np.nextafter(lam*q,np.inf)
    if not kappa:return result
    # q²=u²v²-2uvt²+t⁴. Bound each tensor coefficient outward.
    powers=[]
    for axis in [0,1]:
        aa=arb(box[2*axis]);bb=arb(box[2*axis+1])-aa
        one_lo,_=power_bounds(aa,bb,degrees[axis],1)
        two_lo,_=power_bounds(aa,bb,degrees[axis],2)
        powers.append((np.maximum(0,one_lo),np.maximum(0,two_lo)))
    uvlo=np.maximum(0,np.nextafter(powers[0][0][:,None,None]*powers[1][0][None,:,None],-np.inf))
    _,t2hi=power_bounds(a,b,n,2)
    t4lo,_=power_bounds(a,b,n,4)
    uv_t2_hi=np.nextafter(np.where(t2hi[None,None,:]>=0,uv,uvlo)*t2hi[None,None,:],np.inf)
    q2lo=np.nextafter(powers[0][1][:,None,None]*powers[1][1][None,:,None],-np.inf)
    q2lo=np.nextafter(q2lo-np.nextafter(2*uv_t2_hi,np.inf),-np.inf)
    q2lo=np.nextafter(q2lo+t4lo[None,None,:],-np.inf)
    return np.nextafter(result-np.nextafter(kappa*q2lo,-np.inf),np.inf)
