"""Exact coefficient identities and outward rational finite-bound constants.

Analytic assertions (superharmonic averaging, equilibrium measure, and the
external renormalized-energy theorem) are proved/stated in README.md.
This script verifies their scalar normalization and the eight-point threshold.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from verify_gamma_bounds import I,log_i,exp_i,output_interval

def require(c,m):
 if not c:raise ArithmeticError(m)

def add(a,b):
 out=[F(0)]*max(len(a),len(b))
 for k,v in enumerate(a):out[k]+=v
 for k,v in enumerate(b):out[k]+=v
 return out

def scale(a,s):return [s*v for v in a]
def mul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for k,v in enumerate(a):
  for j,w in enumerate(b):out[k+j]+=v*w
 return out

def same(a,b):
 return add(a,scale(b,-1))==[F(0)]*max(len(a),len(b))

def integrate(a):return sum((v/F(k+1) for k,v in enumerate(a)),F(0))

LOG2=log_i(2)
def logq(x):
 x=F(x);require(x>0,'Positive logarithm argument required');k=0
 while x>=2:x/=2;k+=1
 while x<1:x*=2;k-=1
 return log_i(x)+k*LOG2

def main():
 # n(n-1)/2 as a polynomial in n. Compare coefficients, not sample n.
 N=[F(0),F(-1,2),F(1,2)]
 require(same(add(N,[0,F(3,4),F(-3,4)]),scale(N,F(-1,2))), 'Finite constant coefficient')
 require(same([0,F(-1,2),F(1,2)],N),'Finite log2 coefficient')
 require(same(add(N,[0,F(1,2)]),[0,0,F(1,2)]),'Finite logn coefficient')
 # Exact uniform-disk radial density 2t: second moment and self energy.
 require(integrate([0,0,0,2])==F(1,2),'Disk moment')
 require(integrate([0,1,0,-1])==F(1,4),'Disk self energy at radius one')
 # t*d/dt [h(t)-(1-t)^2] = -(1-t)^2 on (0,1].
 require(same([-1,2,-1],scale(mul([1,-1],[1,-1]),-1)),'Positive-part kernel derivative')
 # The midpoint lens subdisk: (rho/r)^2 * minimum h = (1-u)^4/16.
 p=scale(mul([1,-1],[1,-1]),F(1,4))
 require(same(mul(p,p),scale(mul(mul([1,-1],[1,-1]),mul([1,-1],[1,-1])),F(1,16))),'Lens coefficient')
 gap=40*LOG2-7*logq(7)-14
 ratio=exp_i(gap)
 require(gap.lo>F(1,10) and gap.hi<F(11,100),'Eight-point logarithmic gap')
 require(ratio.lo>F(111017,100000) and ratio.hi<F(111018,100000),'Eight-point ratio enclosure')
 # Reconstruct the unsimplified source formula and compare to closed form.
 rows=[]
 for n in [2,3,5,6,7,8,10,20,100]:
  Nn=F(n*(n-1),2);ln=logq(n);lnprev=logq(n-1)
  A=Nn*(1+ln-lnprev);i0=F(3,4)-LOG2/2;gd=F(3,4)-LOG2/2
  original=A-n*n*i0+F(n,2)*ln+n*gd
  closed=F(n,2)*ln+Nn*(LOG2+ln-lnprev)-Nn/2
  require(not (original.hi<closed.lo or closed.hi<original.lo),'Source/closed interval consistency')
  bridge=A-F(n*n,2)
  require(bridge.lo>F(-1,4)-F(1,12*(n-1)) and bridge.hi<F(-1,4)-F(1,12*n),'Bridge error enclosure')
  rows.append({'n':n,'log_upper':output_interval(closed),'bridge_correction':output_interval(bridge)})
 report={'status':'PASS','arithmetic':'Integer/Fraction polynomial identities and outward rational intervals; no floating-point acceptance',
  'uniform_disk_moment_coefficient':'1/2','uniform_disk_self_energy_constant':'1/4','lens_lower_bound_coefficient':'1/16',
  'finite_bound':'U_n = n^(n/2) * (2*n/(n-1))^(n*(n-1)/2) * exp(-n*(n-1)/4)',
  'eight_point_delta_exact':'40*log(2) - 7*log(7) - 14','eight_point_log_gap':output_interval(gap),
  'eight_point_U_over_candidate':output_interval(ratio),'finite_scalar_spot_checks':rows,
  'global_sharp_eight_point_bound_proved':False,
  'scope':'Universal coefficient identities are checked as polynomials; displayed finite values are scalar checks. The analytic energy inequalities and external asymptotic theorem require the accompanying mathematical proof.'}
 Path(__file__).with_name('finite_bound_certificate.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
 print(json.dumps({k:v for k,v in report.items() if k!='finite_scalar_spot_checks'},indent=2))
if __name__=='__main__':main()
