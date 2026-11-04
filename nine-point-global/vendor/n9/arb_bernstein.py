"""Directed-ball tensor Bernstein coefficients on an arbitrary exact binary box.

For degree N, the coefficient of B_j^N after substituting x=a+b*t into x^i is
sum_{k<=min(i,j)} C(i,k)*a^(i-k)*b^k*C(j,k)/C(N,k).
Apply this linear change of basis separately on the three tensor axes.
NumPy only permutes arrays of Arb objects; all arithmetic is in Arb.
"""
import math
import numpy as np
from flint import arb,arb_mat,fmpq

def transform(n,lo,hi):
    a=arb(lo);b=arb(hi)-a
    return arb_mat([[sum((math.comb(i,k)*a**(i-k)*b**k*arb(fmpq(math.comb(j,k),math.comb(n,k)))
                         for k in range(min(i,j)+1)),arb(0))for i in range(n+1)]for j in range(n+1)])

def coefficients(poly,degrees,box):
    shape=tuple(d+1 for d in degrees);arr=np.empty(shape,dtype=object);arr.fill(arb(0))
    for e,c in poly.items():arr[e]=c
    for axis,n in enumerate(degrees):
        moved=np.moveaxis(arr,axis,0);old_shape=moved.shape;flat=moved.reshape(n+1,-1)
        transformed=transform(n,box[2*axis],box[2*axis+1])*arb_mat(flat.tolist())
        values=np.array([[transformed[i,j]for j in range(flat.shape[1])]for i in range(n+1)],dtype=object)
        arr=np.moveaxis(values.reshape(old_shape),0,axis)
    return arr

def cubic_log_upper_power(box):
    center=[(box[2*k]+box[2*k+1])/2 for k in range(3)]
    d0=arb(max(center[0]+center[1]-2*center[2],.25))
    delta={(0,0,0):arb(box[0])+arb(box[2])-2*arb(box[4])-d0,
           (1,0,0):arb(box[1])-arb(box[0]),(0,1,0):arb(box[3])-arb(box[2]),
           (0,0,1):-2*(arb(box[5])-arb(box[4]))}
    out={(0,0,0):d0.log()};power={(0,0,0):arb(1)}
    for m in range(1,4):
        product={}
        for e,c in power.items():
            for f,a in delta.items():
                ex=tuple(i+j for i,j in zip(e,f));product[ex]=product.get(ex,arb(0))+c*a
        power=product;factor=arb((-1)**(m+1))/(m*d0**m)
        for ex,c in power.items():out[ex]=out.get(ex,arb(0))+factor*c
    return out

def residual_coefficients(poly,degrees,box,order=3):
    arr=coefficients(poly,degrees,box)
    if order==3:
        upper=coefficients(cubic_log_upper_power(box),degrees,(0.,1.,0.,1.,0.,1.))
        for ix in np.ndindex(arr.shape):arr[ix]=arr[ix]-upper[ix]
        return arr
    assert order==1
    center=[(box[2*k]+box[2*k+1])/2 for k in range(3)]
    d0=arb(max(center[0]+center[1]-2*center[2],.25))
    grids=[[arb(box[2*k])+(arb(box[2*k+1])-arb(box[2*k]))*arb(fmpq(j,n))
            for j in range(n+1)]for k,n in enumerate(degrees)]
    constant=d0.log()-1
    for ix in np.ndindex(arr.shape):
        linear=(grids[0][ix[0]]+grids[1][ix[1]]-2*grids[2][ix[2]])/d0
        arr[ix]=arr[ix]-constant-linear
    return arr

def strictly_positive(poly,degrees,box):
    arr=residual_coefficients(poly,degrees,box)
    return all(x>0 for x in arr.flat)

def coefficients_and_cubic_residual(poly,degrees,box):
    original=coefficients(poly,degrees,box)
    upper=coefficients(cubic_log_upper_power(box),degrees,(0.,1.,0.,1.,0.,1.))
    residual=np.empty(original.shape,dtype=object)
    for ix in np.ndindex(original.shape):residual[ix]=original[ix]-upper[ix]
    return original,residual
