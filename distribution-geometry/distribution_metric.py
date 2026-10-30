"""Explicit full-plane W1 distribution bounds from the Coulomb deficit.

Classical Coulomb transport framework; disk constants and normalization are
proved in DISTRIBUTION_METRIC.md. All accepted arithmetic is rational.
W1 bounds are necessary, not proofs of any lattice or global optimality.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'joint-bounds'))
from subset_moment import log_interval,outward_decimal


def require_n(n):
    if type(n)is not int or n<2:
        raise ValueError('integer n>=2 required')


def sqrt_bounds(q,scale=10**12):
    q=F(q)
    if q<0 or type(scale)is not int or scale<1:
        raise ValueError('nonnegative rational and positive root scale required')
    v=isqrt(q.numerator*scale*scale//q.denominator)
    exact=v*v*q.denominator==q.numerator*scale*scale
    return F(v,scale),F(v if exact else v+1,scale)


@lru_cache(maxsize=16)
def pi_interval(terms=14):
    if type(terms)is not int or terms<1:
        raise ValueError('positive integer arctangent length required')
    def atan(q):
        value=sum(((-1)**j*q**(2*j+1)/(2*j+1) for j in range(terms)),F(0))
        nxt=(-1)**terms*q**(2*terms+1)/(2*terms+1)
        return min(value,value+nxt),max(value,value+nxt)
    a,b=atan(F(1,5))
    c,d=atan(F(1,239))
    return 16*a-4*d,16*b-4*c


def smear_penalty_upper(h,terms=28):
    h=F(h)
    if h<=0:
        raise ValueError('positive smear multiplier required')
    lo,hi=log_interval(h,terms)
    return max(F(0),(h*h-1)/2-lo)


def quantization_lower(n,scale=10**12):
    """Universal W1 lower bound for ANY n equally weighted atoms vs mu0."""
    require_n(n)
    return F(2,3)*sqrt_bounds(F(2,n),scale)[0]


def transport_bound(n,proven_gap_upper,centre_radius_upper,h=1,
                    original_z_scale=True,scale=10**12):
    """W1(empirical points, uniform sqrt(2) disk), under proved hypotheses.

    The caller must prove delta<=gap and every COULOMB center |x_i|<=rho.
    Merely passing positive numbers does not prove those assumptions.
    R includes BOTH the full background disk and every smoothed point disk.
    """
    require_n(n)
    gap,rho,h=map(F,(proven_gap_upper,centre_radius_upper,h))
    if gap<0 or rho<=0 or h<=0:
        raise ValueError('nonnegative proved gap and positive radius/smear required')
    if rho*rho<F(n-1,n):
        raise ValueError('radius upper bound incompatible with prescribed mean square')
    r=sqrt_bounds(2*h*h/n,scale)[1]
    R=max(sqrt_bounds(F(2),scale)[1],rho+r)
    penalty=n*smear_penalty_upper(h)
    energy_gap=gap+penalty
    electrostatic=R*sqrt_bounds(energy_gap,scale)[1]/n
    smoothing=2*r/3
    conversion=1-sqrt_bounds(F(n-1,n),scale)[0] if original_z_scale else F(0)
    energy_upper=smoothing+electrostatic+conversion
    # Independent coupling also gives W1<=sqrt(second moments sum), because
    # the reference disk has zero mean. Its second moment is exactly one.
    moment_upper=sqrt_bounds(F(2) if original_z_scale else F(2*n-1,n),scale)[1]
    upper=min(energy_upper,moment_upper)
    return {'n':n,'h':h,'proven_gap_upper':gap,'centre_radius_upper':rho,
            'full_support_radius_upper':R,'smear_radius_upper':r,
            'smear_penalty_upper':penalty,'smoothing_W1_upper':smoothing,
            'electrostatic_W1_upper':electrostatic,'scale_conversion_upper':conversion,
            'energy_transport_W1_upper':energy_upper,'independent_coupling_W1_upper':moment_upper,
            'W1_upper':upper,'universal_quantization_W1_lower':quantization_lower(n,scale),
            'point_scale':'original z' if original_z_scale else 'Coulomb x',
            'scope':'Necessary full-plane W1 bound under separately proved energy gap and support hypotheses; no optimizer or lattice is certified.'}


def effective_global_gap(n,scale=10**12):
    """Previously proved attainable-polar/disk deficit upper, valid n>=6."""
    require_n(n)
    if n<6:
        raise ValueError('this existing effective polar remainder needs n>=6')
    return F(n,100)+F(77,6)*sqrt_bounds(F(n),scale)[1]+F(35,4)


def global_transport_bound(n,proven_gap_upper=None,h_values=(F(1),F(97,100)),
                           scale=10**12):
    """For GLOBAL maximizers, use their proved fixed radius cap.

    A supplied smaller gap must be separately proved, e.g. from an exact
    attainable lower configuration. With no input, use effective_global_gap.
    """
    require_n(n)
    gap=effective_global_gap(n,scale) if proven_gap_upper is None else F(proven_gap_upper)
    hs=tuple(map(F,h_values))
    if not hs or any(h<=0 for h in hs):
        raise ValueError('nonempty positive smearing choices required')
    rho=F(n-1,n)*sqrt_bounds(F(2),scale)[1]
    rows=[transport_bound(n,gap,rho,h,scale=scale) for h in hs]
    best=min(rows,key=lambda row:row['W1_upper'])
    return {'n':n,'W1_upper':best['W1_upper'],'chosen_h':best['h'],
            'universal_quantization_W1_lower':best['universal_quantization_W1_lower'],
            'proven_gap_upper':gap,'candidates':rows,
            'scope':'Necessary for every global maximizer, using the established subset radius cap and the specified proved energy budget.'}


def asymptotic_transport_certificate(h=F(97,100),gamma_gap_upper=F(1,800),scale=10**12):
    """Conditional on the repo's external linear-energy/gamma-gap theorem.

    Gives a limsup coefficient, not a finite-n onset for that coefficient.
    """
    h,gap=F(h),F(gamma_gap_upper)
    if h<=0 or gap<0:
        raise ValueError('positive smear and nonnegative limiting gap bound required')
    psi=smear_penalty_upper(h)
    root2lo,root2hi=sqrt_bounds(F(2),scale)
    upper=F(2,3)*h*root2hi+sqrt_bounds(2*(gap+psi),scale)[1]
    return {'h':h,'limiting_gap_upper':gap,'psi_upper':psi,
            'limsup_sqrt_n_W1_upper':upper,
            'universal_liminf_sqrt_n_W1_lower':F(2,3)*root2lo,
            'headline_upper':F(981,1000),
            'headline_certified':upper<F(981,1000),
            'scope':'Asymptotic only; uses the pre-existing linear energy theorem and triangular trial bound. No effective starting n or optimal lattice is inferred.'}


def cross_axes_background_interval(scale=10**12):
    """E_mu0 min(|Re w|,|Im w|)=8*(sqrt(2)-1)/(3*pi)."""
    lo,hi=sqrt_bounds(F(2),scale)
    plo,phi=pi_interval()
    return 8*(lo-1)/(3*phi),8*(hi-1)/(3*plo)


def lipschitz_class_check(empirical_interval,background_interval,lipschitz_upper,W1_upper):
    """Reject a whole class if one proved Lipschitz statistic is incompatible."""
    elo,ehi=map(F,empirical_interval)
    blo,bhi=map(F,background_interval)
    L,w=F(lipschitz_upper),F(W1_upper)
    if elo>ehi or blo>bhi or L<0 or w<0:
        raise ValueError('ordered statistic intervals and nonnegative bounds required')
    separation=max(F(0),elo-bhi,blo-ehi)
    return {'statistic_separation_lower':separation,'allowed_separation_upper':L*w,
            'excluded':separation>L*w,
            'scope':'Conditional class exclusion; statistic intervals and Lipschitz constant must be proved independently.'}


def box_population(n,box,W1_upper,eta):
    """Safe point-count interval for an axis-aligned ORIGINAL-z rectangle.

    The input W1 upper must compare empirical z to the uniform sqrt(2) disk.
    An expanded rectangle bounds the Euclidean neighborhood from above.
    An eroded rectangle inside the background disk gives a mass lower.
    """
    require_n(n)
    left,right,bottom,top=map(F,box)
    w,eta=F(W1_upper),F(eta)
    if left>=right or bottom>=top or w<0 or eta<=0:
        raise ValueError('nonempty rectangle, nonnegative W1 and positive margin required')
    plo,phi=pi_interval()
    outer_area=(right-left+2*eta)*(top-bottom+2*eta)
    outer_upper=min(F(1),outer_area/(2*plo))
    il,ir,ib,it=left+eta,right-eta,bottom+eta,top-eta
    inner_lower=F(0)
    if il<ir and ib<it and max(x*x+y*y for x in (il,ir) for y in (ib,it))<=2:
        inner_lower=(ir-il)*(it-ib)/(2*phi)
    lower=max(F(0),inner_lower-w/eta)
    upper=min(F(1),outer_upper+w/eta)
    lo=n*lower
    hi=n*upper
    count_lower=-((-lo.numerator)//lo.denominator)
    count_upper=hi.numerator//hi.denominator
    return {'n':n,'box':[left,right,bottom,top],'margin':eta,
            'population_lower':max(0,count_lower),'population_upper':min(n,count_upper),
            'eroded_background_mass_lower':inner_lower,
            'expanded_background_mass_upper':outer_upper,
            'scope':'Necessary count range from full-plane W1. Thin or small boxes may yield [0,n]; no full grid transport plan is claimed.'}


if __name__=='__main__':
    import argparse
    import json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n',type=int)
    args=parser.parse_args()
    report=global_transport_bound(args.n)
    for key in ('W1_upper','universal_quantization_W1_lower'):
        report[key]=str(outward_decimal((report[key],report[key]))[1 if key=='W1_upper' else 0])
    del report['candidates']
    print(json.dumps(report,default=str,indent=2))
