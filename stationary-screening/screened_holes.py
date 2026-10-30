"""Exact joint Gram of harmonic moments, compact modes and disjoint holes.

All geometry and moments use Coulomb x scale; h=1 throughout.
The support and empty-smoothed-ball hypotheses are not inferred from input.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
from stationary_moments import centered_moments,compact_measurements,cmul,norm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'radial-angular'))
from mixed_energy import solve


def geometry(S,holes):
    S=F(S)
    if S<2:raise ValueError('support must include the background disk')
    rows=[]
    for x,y,b in holes:
        x,y,b=map(F,(x,y,b))
        if b<=0:raise ValueError('positive energy-hole radius required')
        q=x*x+y*y;remaining=2-q-b*b
        if remaining<0 or remaining*remaining<4*q*b*b:
            raise ValueError('every energy hole must lie in the background disk')
        for a,c,d in rows:
            if (x-a)**2+(y-c)**2<(b+d)**2:
                raise ValueError('hole interiors overlap')
        rows.append((x,y,b))
    if sum((b**4 for x,y,b in rows),F(0))>=S*S:
        raise ValueError('degenerate radial conditioning')
    return S,rows


def conditioned_energy(n,S,holes,moments=None):
    """Joint necessary energy; moments=None uses only second-moment matching.

    Supplying actual centered moments invokes the stationary mixed closure
    for every supplied mode. Even [p1=0] adds the stationary dipole constraint.
    All given holes must be empty of the SMEARED measure; for point holes
    their radii must first be reduced by sqrt(2/n).
    """
    if type(n)is not int or n<2:raise ValueError('integer n>=2 required')
    S,holes=geometry(S,holes)
    values=[] if moments is None else centered_moments(n,moments)
    compact=[] if moments is None else compact_measurements(n,S,values)
    diagonal=[S*S]
    data=[F(0)]
    for k,v in enumerate(compact,1):
        diagonal.extend([S**(k+2)/F(k+2)]*2)
        data.extend(v)
    size=len(diagonal)
    C=[[diagonal[i] if i==j else F(0) for j in range(size)] for i in range(size)]
    s=[F(0)]*size
    T=F(0)
    for x,y,b in holes:
        t=b**4;T+=t
        v=[F(1)];power=(F(1),F(0))
        for k in range(1,len(values)+1):
            power=cmul(power,(x,y))
            v.extend(((k+1)*power[0],(k+1)*power[1]))
        for i in range(size):
            s[i]+=t*v[i]
            for j in range(size):C[i][j]-=t*v[i]*v[j]
    residual=[a+F(n,4)*b for a,b in zip(data,s)]
    solution=solve(C,residual)
    extra=sum((a*b for a,b in zip(residual,solution)),F(0))
    basic=n*n*T/16
    harmonic=sum((norm(p)/(k*S**k) for k,p in enumerate(values,1)),F(0))
    return {'n':n,'support_squared':S,'holes':holes,'screening_modes':len(values),
            'basic_hole_energy':basic,'conditioning_gain':extra,
            'harmonic_energy':harmonic,'joint_energy_lower':basic+extra+harmonic,
            'Schur_matrix':C,'Schur_data':residual,
            'scope':'Second-moment matching, disjoint empty smoothed balls and proven support are required; modes additionally require stationarity and actual Coulomb moments. Full radial energy is not independently added.'}


def one_hole_dipole(n,S,center,radius):
    if len(center)!=2:raise ValueError('rational center pair required')
    return conditioned_energy(n,S,[(center[0],center[1],radius)],[(0,0)])
