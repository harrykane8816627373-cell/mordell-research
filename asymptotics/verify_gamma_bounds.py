"""Outward rational enclosures for the Mordell linear-coefficient bounds.

Uses only integer/Fraction arithmetic. The triangular lattice eta product is
imported mathematically from Betermin-Sandier, arXiv:1404.4485v5, Section7.3.
This verifies constants; it does not prove the Coulomb-gas asymptotic theorem.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
D=70
S=10**D

def floor_f(x):return x.numerator//x.denominator

def ceil_f(x):return -floor_f(-x)

class I:
 def __init__(self,lo,hi=None):
  lo=F(lo);hi=lo if hi is None else F(hi)
  if lo>hi:raise ArithmeticError('Reversed interval')
  self.lo=F(floor_f(lo*S),S);self.hi=F(ceil_f(hi*S),S)
 @staticmethod
 def cv(x):return x if isinstance(x,I) else I(x)
 def __add__(self,x):
  x=self.cv(x);return I(self.lo+x.lo,self.hi+x.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,x):return self+-self.cv(x)
 def __rsub__(self,x):return self.cv(x)+-self
 def __mul__(self,x):
  x=self.cv(x);v=[a*b for a in (self.lo,self.hi) for b in (x.lo,x.hi)]
  return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(self,x):
  x=self.cv(x)
  if x.lo<=0<=x.hi:raise ArithmeticError('Interval denominator contains zero')
  return self*I(1/x.hi,1/x.lo)
 def __rtruediv__(self,x):return self.cv(x)/self
 def power(self,k):
  y=I(1)
  for _ in range(k):y=y*self
  return y

def atan_recip(q,terms=100):
 """Alternating series, with the first omitted term bounding the error."""
 x=F(1,q);v=F(0)
 for k in range(terms):v+=(-1)**k*x**(2*k+1)/(2*k+1)
 err=x**(2*terms+1)/(2*terms+1)
 return I(v-err,v+err)

def log_i(x,terms=180):
 x=I.cv(x)
 if x.lo<=0:raise ArithmeticError('Logarithm argument not positive')
 t=(x-1)/(x+1);a=max(abs(t.lo),abs(t.hi))
 if a>=1:raise ArithmeticError('Logarithm series does not contract')
 t2=t*t;power=t;v=I(0)
 for k in range(terms):v=v+power/F(2*k+1);power=power*t2
 # The exact remainder is bounded using the interval bound on t^(2K+1).
 # Including fixed-point rounding keeps power.lo/power.hi outward.
 tail=2*max(abs(power.lo),abs(power.hi))/((2*terms+1)*(1-a*a))
 return 2*v+I(-tail,tail)

def exp_i(x,terms=180):
 x=I.cv(x)
 if x.lo<0:raise ArithmeticError('This positive-series routine requires x>=0')
 term=I(1);v=I(1)
 for k in range(1,terms+1):term=term*x/k;v=v+term
 ratio=x.hi/F(terms+2)
 if ratio>=1:raise ArithmeticError('Exponential tail ratio does not contract')
 tail=term.hi*x.hi/F(terms+1)/(1-ratio)
 return I(v.lo,v.hi+tail)

def sqrt_i(x):
 x=I.cv(x)
 if x.lo<0:raise ArithmeticError('Negative square-root argument')
 low=isqrt(floor_f(x.lo*S*S));high=isqrt(ceil_f(x.hi*S*S))+1
 return I(F(low,S),F(high,S))

def fixed_decimal(x,places,up=False):
 a=(ceil_f if up else floor_f)(x*10**places)
 sign='-' if a<0 else '';a=abs(a);whole,part=divmod(a,10**places)
 return sign+str(whole)+'.'+str(part).zfill(places)

def output_interval(x,places=30):
 return {'lower':fixed_decimal(x.lo,places),'upper':fixed_decimal(x.hi,places,True)}

def main():
 pi=16*atan_recip(5)-4*atan_recip(239)
 root3=I(F(isqrt(3*S*S),S),F(isqrt(3*S*S)+1,S))
 a=pi*root3
 rho=1/exp_i(a)
 if not F(4,1000)<rho.lo<rho.hi<F(5,1000):raise ArithmeticError('Unexpected eta nome')
 K=16;partial=I(0)
 for k in range(1,K+1):partial+=log_i(1-((-1)**k)*rho.power(k))
 tail=2*rho.hi**(K+1)/(1-rho.hi)**2
 gamma_tri=log_i(pi)/2+log_i(3)/4-a/12+2*partial+I(-tail,tail)
 gamma_upper=log_i(pi)/2
 gamma_disk=I(F(3,4))-log_i(2)/2
 if not gamma_tri.hi<gamma_disk.lo<gamma_disk.hi<gamma_upper.lo:raise ArithmeticError('Disk-smearing bound ordering failed')
 # Betermin-Sandier Remark7.6, importing the RSZ lower bound for W_min.
 root27=sqrt_i(27);root2pi=sqrt_i(2*pi);root_sum=sqrt_i(2*pi+root27)
 rsz_a=2*root2pi/root27*(root_sum+root2pi)
 rsz_b=(root_sum-root2pi)/(root_sum+root2pi)
 gamma_rsz=gamma_upper+rsz_b*log_i(1-1/exp_i(rsz_a))/2
 if not gamma_disk.hi<gamma_rsz.lo<gamma_rsz.hi<gamma_upper.lo:raise ArithmeticError('RSZ improvement check failed')
 if not F(4021788952,10**10)<gamma_tri.lo<gamma_tri.hi<F(4021788953,10**10):raise ArithmeticError('Triangular constant check failed')
 if not F(5723649429,10**10)<gamma_upper.lo<gamma_upper.hi<F(5723649430,10**10):raise ArithmeticError('Upper constant check failed')
 if not gamma_tri.hi<gamma_upper.lo:raise ArithmeticError('Bound ordering check failed')
 for bound in (gamma_tri,gamma_disk,gamma_upper,gamma_rsz):
  if not F(fixed_decimal(bound.lo,30))<bound.lo<=bound.hi<F(fixed_decimal(bound.hi,30,True)):raise ArithmeticError('Printed decimal endpoints were not strictly outward')
 report={'status':'PASS','arithmetic':'Outward fixed-denominator Fraction intervals; no floating-point acceptance',
  'precision_decimal_digits':D,'strictly_outward_printed_decimals':True,'eta_product_terms':K,'log_series_terms':180,'exp_series_terms':180,
  'triangular_lower_endpoint':output_interval(gamma_tri),'upper_endpoint_half_log_pi':output_interval(gamma_upper),
  'stronger_disk_smearing_upper_endpoint':output_interval(gamma_disk),
  'gamma_disk_exact':'3/4 - log(2)/2',
  'stronger_RSZ_upper_endpoint':output_interval(gamma_rsz),
  'RSZ_a':output_interval(rsz_a),'RSZ_b':output_interval(rsz_b),
  'pi':output_interval(pi),'sqrt3':output_interval(root3),'nome_exp_minus_pi_sqrt3':output_interval(rho),
  'gamma_tri_exact':'3*log(Gamma(1/3)) + log(3)/2 - 2*log(2) - 3*log(pi)/2',
  'source':'https://arxiv.org/html/1404.4485v5#S7.SS3',
  'interpretation':'The exact triangular endpoint lies in its certified interval and gamma is at least this endpoint. Gamma is at most the elementary disk-smearing endpoint 3/4-log(2)/2, strictly below the RSZ endpoint and half(log(pi)). The exact universal gamma is not computed.'}
 Path(__file__).with_name('gamma_bounds_certificate.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
