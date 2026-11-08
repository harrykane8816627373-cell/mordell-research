"""Bounded FLOAT Bernstein exploration, NOT an accepted continuous bound."""
from pathlib import Path
from math import comb,log
import argparse,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

def mul(a,b):return np.polynomial.polynomial.polymul(a,b)
def matrix_features(data,k):
 s=data['feature_order']
 if k==0:
  nodes=np.array([float.fromhex(a) for a in data['radial_nodes']]);means=np.array([float.fromhex(a) for a in data['radial_moment_centers']]);q=np.polynomial.polynomial.polyfromroots(nodes)
  M=np.zeros((2*len(nodes)+1,s))
  for j in range(1,s):M[j,j-1]=1;M[0,j-1]=-means[j-1]
  M[:,-1]=10*mul(q,q)
  return M
 raw=np.zeros((s,s))
 for j in range(s):
  for i in range(j+1):raw[i,j]=comb(j,i)*(-1)**(j-i)
 E=np.array([[float.fromhex(a) for a in row] for row in data['transforms'][k-1]])
 return raw@E

def coefficients(data):
 modes=data['harmonics_through'];degree=modes+2*(data['feature_order']-1)
 C=np.zeros((degree+1,degree+1,modes+1));H=[{(0,0,0):1.},{(0,0,1):1.}]
 for k in range(2,modes+1):
  q={}
  for (a,b,c),v in H[-1].items():q[(a,b,c+1)]=q.get((a,b,c+1),0)+2*v
  for (a,b,c),v in H[-2].items():q[(a+1,b+1,c)]=q.get((a+1,b+1,c),0)-v
  H.append(q)
 n=data['n'];m=n-1;B=float.fromhex(data['B']);b0=(B-n*m/2)/n
 C[0,0,0]=2*b0/m;C[1,0,0]=C[0,1,0]=n/(2*m)
 for k,rawQ in enumerate(data['Q_blocks']):
  Q=np.array([[float.fromhex(a) for a in row] for row in rawQ]);F=matrix_features(data,k);G=F@Q@F.T
  for i in range(G.shape[0]):
   for j in range(G.shape[1]):
    for (a,b,c),v in H[k].items():C[a+i,b+j,c]-=v*G[i,j]
    C[i+j+k,0,0]-=G[i,j]/(2*m);C[0,i+j+k,0]-=G[i,j]/(2*m)
 return C

def transform(deg,lo,hi):
 T=np.zeros((deg+1,deg+1));D=np.zeros_like(T)
 for i in range(deg+1):
  for j in range(i+1):T[j,i]=comb(i,j)*lo**(i-j)*(hi-lo)**j;D[i,j]=comb(i,j)/comb(deg,j)
 return D@T

def bernstein(C,box):
 out=C
 for axis in range(3):
  T=transform(C.shape[axis]-1,*box[axis]);out=np.moveaxis(np.tensordot(T,out,axes=(1,axis)),0,axis)
 return out

def split(B,axis):
 a=np.moveaxis(B,axis,0).copy();left=np.empty_like(a);right=np.empty_like(a);n=a.shape[0]-1
 left[0]=a[0];right[n]=a[n]
 for j in range(1,n+1):
  a[:n-j+1]=(a[:n-j+1]+a[1:n-j+2])/2;left[j]=a[0];right[n-j]=a[n-j]
 return np.moveaxis(left,0,axis),np.moveaxis(right,0,axis)

def dcoeff(box,shape):
 u=np.linspace(*box[0],shape[0])[:,None,None];v=np.linspace(*box[1],shape[1])[None,:,None];t=np.linspace(*box[2],shape[2])[None,None,:]
 return u+v-2*t

