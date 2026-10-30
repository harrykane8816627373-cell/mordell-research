"""Stationary mixed-moment closure and exact angular energy certificates.

Moments supplied to the stationary interface are Coulomb x moments.
Stationarity is a hypothesis, never inferred from the moments or small errors.
"""
from fractions import Fraction as F
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'joint-bounds'))
from subset_moment import log_interval


def cadd(z,w):return z[0]+w[0],z[1]+w[1]
def csub(z,w):return z[0]-w[0],z[1]-w[1]
def cmul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def cscale(z,a):return a*z[0],a*z[1]
def norm(z):return z[0]*z[0]+z[1]*z[1]


def index(n,k):
    if type(n)is not int or n<2 or type(k)is not int or not 1<=k<=64:
        raise ValueError('integer n>=2 and angular degree 1..64 required')


def centered_moments(n,moments):
    values=[tuple(map(F,z)) for z in moments]
    if not values or len(values)>64 or any(len(z)!=2 for z in values):
        raise ValueError('1..64 complex rational power sums required')
    index(n,len(values))
    if values[0]!=(0,0):
        raise ValueError('stationary Coulomb moments must be centered')
    return values


def convolution(moments,k):
    return tuple(sum((cmul(moments[j-1],moments[k-j-1])[a]
                      for j in range(1,k)),F(0)) for a in range(2))


def compact_measurements(n,support_squared,moments):
    """n times actual averages of x^k(S-|x|^2)_+, at h=1.

    Requires a stationary normalized configuration and support of every
    smeared disk inside |w|^2<=S. Returns measurements, not residual tests
    of whether an arbitrary input configuration is stationary.
    """
    values=centered_moments(n,moments)
    S=F(support_squared)
    if S<=0:raise ValueError('positive squared support required')
    return [csub(cscale(p,S-2),cscale(convolution(values,k),F(1,n)))
            for k,p in enumerate(values,1)]


def stationary_energy(n,support_squared,moments):
    S=F(support_squared)
    values=centered_moments(n,moments)
    a=compact_measurements(n,S,values)
    rows=[]
    for k,(p,v) in enumerate(zip(values,a),1):
        old=norm(p)/(k*S**k)
        extra=(k+2)*norm(v)/S**(k+2)
        rows.append({'k':k,'holomorphic_energy':old,'stationary_extra':extra})
    return {'n':n,'support_squared':S,'rows':rows,
            'angular_energy_lower':sum((r['holomorphic_energy']+r['stationary_extra'] for r in rows),F(0)),
            'scope':'Conditional on actual Coulomb moments, stationarity and the supplied support; no optimum inferred.'}


def low_mode_coefficient(n,k,support_squared):
    """Original-z coefficient for k=2,3, where the convolution vanishes."""
    index(n,k)
    S=F(support_squared)
    if k not in (2,3) or S<=0:raise ValueError('k=2 or 3 and S>0 required')
    return F(n-1,n)**k*(1/(k*S**k)+(k+2)*(S-2)**2/S**(k+2))


def residual_mode_lower(k,S,power_modulus_lower,weighted_residual_upper,n):
    """Robust k=2/3 bound for centered, possibly nonstationary inputs.

    weighted_residual_upper bounds |sum x_i^(k+1) R_i|, with
    R_i=sum_{j!=i}1/(x_i-x_j)-(n/2)*conj(x_i). Uses h=1.
    """
    index(n,k)
    S,p,e=map(F,(S,power_modulus_lower,weighted_residual_upper))
    if k not in (2,3) or S<=0 or p<0 or e<0:
        raise ValueError('k=2,3, positive S and nonnegative bounds required')
    extra=max(F(0),abs(S-2)*p-2*e/n)
    return p*p/(k*S**k)+(k+2)*extra*extra/S**(k+2)


def eliminated_convolution_coefficient(n,k,S):
    """Original-z coefficient of |sum P_l P_(k-l)|^2, minimizing over P_k."""
    index(n,k)
    S=F(S)
    if S<=0:raise ValueError('positive squared support required')
    return F(n-1,n)**k*(k+2)/(n*n*S**k*(S*S+k*(k+2)*(S-2)**2))


def eight_point_certificate():
    S=F(151,50)
    # 2*(7/8+1/sqrt(8))^2 = 57/32+7/(4sqrt(2)).
    margin=S-F(57,32)
    if not margin>0 or not 2*margin*margin>F(49,16):
        raise ArithmeticError('eight-point global smeared support')
    l2,h2=log_interval(2);l7,h7=log_interval(7)
    gap_upper=40*h2-7*l7-14
    c2,c3=[low_mode_coefficient(8,k,S) for k in (2,3)]
    cap=F(64,49)
    if not gap_upper<F(523,5000)<c2*cap or not c3>F(27,100)*c2:
        raise ArithmeticError('eight-point improved ellipse')
    return {'n':8,'support_squared':S,'gap_upper':gap_upper,
            'simple_gap_strict_upper':F(523,5000),
            'p2_coefficient':c2,'p3_coefficient':c3,
            'p2_fourth_power_coefficient':eliminated_convolution_coefficient(8,4,S),
            'p2_squared_p3_squared_coefficient':4*eliminated_convolution_coefficient(8,5,S),
            'ellipse':'|p2(z)|^2+(27/100)|p3(z)|^2 < 64/49',
            'p2_modulus_strict_upper':F(8,7),
            'scope':'Every global eight-point optimizer; not the full sharp product bound or the p2*p3!=0 classification.'}
