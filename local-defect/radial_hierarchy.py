"""Matched-moment radial Gram hierarchy, orthogonal to angular modes."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'radial-angular'))
from mixed_energy import solve


def radial_gram(degree,S):
    if type(degree)is not int or not 1<=degree<=12 or F(S)<=0:
        raise ValueError('degree 1..12 and positive squared support required')
    S=F(S)
    return [[F(2*j*k,(j+k)*(j+k-1))*S**(j+k)
             for k in range(1,degree+1)] for j in range(1,degree+1)]


def radial_data(n,radii_squared,degree,S):
    if type(n)is not int or n<2:raise ValueError('integer n>=2 required')
    H=radial_gram(degree,S);S=F(S)
    q=list(map(F,radii_squared));r2=F(2,n)
    if len(q)!=n or min(q)<0 or sum(q)!=n-1 or S<2:
        raise ValueError('Coulomb squared radii with sum n-1 and S>=2 required')
    margin=S-max(q)-r2
    if margin<0 or margin*margin<4*max(q)*r2:
        raise ValueError('S must cover every full smeared disk')
    moments=[sum((sum((F(comb(m,l)**2,l+1)*r2**l*v**(m-l)
                           for l in range(m+1)),F(0)) for v in q),F(0))
             for m in range(degree+1)]
    differences=[moments[m]-F(n*2**m,m+1) for m in range(degree+1)]
    if differences[0]!=0 or differences[1]!=0:raise ArithmeticError('mass/second moment matching failed')
    measurements=[sum(((-1)**m*comb(j,m)*S**(j-m)*differences[m]
                       for m in range(j+1)),F(0)) for j in range(1,degree+1)]
    solution=solve(H,measurements)
    energy=sum((a*b for a,b in zip(measurements,solution)),F(0))
    return {'n':n,'degree':degree,'support_squared':S,'smeared_moment_sums':moments,
        'radial_measurements':measurements,'radial_energy_lower':energy,
        'scope':'Lower bound on n^2 D(radial nu), at h=1 and stated support. Can add orthogonal angular energy; not an independent extra copy of radial CDF or hole energy.'}


def eight_radial_coefficient(S=F(151,50)):
    S=F(S)
    if S<=0:raise ValueError('S>0 required')
    return F(9,2)*F(49,64)**2/S**4
