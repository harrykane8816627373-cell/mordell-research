"""Outward float upper coefficients of the GLOBAL cubic Taylor upper bound to log.

The fourth derivative of log is negative, so its degree-three Taylor polynomial
at any d0>0 lies above log(d) for every d>0, on either side of d0.
All coefficient arithmetic is first performed in Arb. Nonnegative exact
Bernstein conversion weights are then enclosed with outward float operations.
"""
import math
import numpy as np
from flint import arb,fmpq
from arb_bernstein import cubic_log_upper_power

class CubicLogUpper:
    def __init__(self,degrees):
        assert min(degrees)>=3
        self.degrees=degrees;self.shape=tuple(n+1 for n in degrees);self.weights={}
        intervals=[]
        for n in degrees:
            rows=[]
            for k in range(4):
                exact=[arb(fmpq(math.comb(j,k),math.comb(n,k)))if j>=k else arb(0)for j in range(n+1)]
                lo=np.maximum(0,np.array([np.nextafter(float(v.lower()),-np.inf)for v in exact]))
                hi=np.array([np.nextafter(float(v.upper()),np.inf)for v in exact]);rows.append((lo,hi))
            intervals.append(rows)
        for i in range(4):
            for j in range(4-i):
                for k in range(4-i-j):
                    (ul,uh),(vl,vh),(tl,th)=intervals[0][i],intervals[1][j],intervals[2][k]
                    lo=np.maximum(0,np.nextafter(np.nextafter(ul[:,None,None]*vl[None,:,None],-np.inf)*tl[None,None,:],-np.inf))
                    hi=np.nextafter(np.nextafter(uh[:,None,None]*vh[None,:,None],np.inf)*th[None,None,:],np.inf)
                    self.weights[(i,j,k)]=(lo,hi)
    def upper(self,box):
        result=np.zeros(self.shape)
        for ex,c in cubic_log_upper_power(box).items():
            hi=np.nextafter(float(c.upper()),np.inf);wl,wh=self.weights[ex]
            term=np.nextafter(hi*(wh if hi>=0 else wl),np.inf)
            result=np.nextafter(result+term,np.inf)
        return result
