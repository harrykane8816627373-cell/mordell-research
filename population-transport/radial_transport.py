"""Truncated inverse-radius transport joined to same-circle capacities.

The analytic proof is RADIAL_TRANSPORT.md. All acceptance uses Fraction.
Energy transport is in Coulomb x scale; the independent CDF comparison
uses original z scale and says so explicitly in its output.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'radial-angular'))
sys.path.insert(0,str(ROOT/'configuration-bounds'))
sys.path.insert(0,str(ROOT/'collective-stability'))
from radial_layers import smearing_penalty_upper
from radial_distribution import root_upper,cdf_error_upper,polar_gap_upper
from ring_population import pi_interval


def require_n(n):
    if type(n) is not int or n<2:
        raise ValueError('integer n>=2 required')


def sqrt_lower(q,scale=10**12):
    q=F(q)
    if q<0 or type(scale) is not int or scale<1:
        raise ValueError('nonnegative rational and positive root scale required')
    return F(isqrt(q.numerator*scale*scale//q.denominator),scale)


def ceil_fraction(q):
    q=F(q)
    return -((-q.numerator)//q.denominator)


def energy_mean_lower(n,proven_gap_upper,a,h=1,nonradial_delta_lower=0,scale=10**12):
    """Lower bound on n^-1 sum min(1/a,1/|x_i|), Q(x)=n-1.

    The radial Dirichlet budget is delta+n*psi(h)-nonradial. A nonzero
    nonradial lower input needs independent proof for this SAME h and an
    orthogonal angular energy; a second radial lower cannot be subtracted.
    """
    require_n(n)
    gap,a,h,angular=map(F,(proven_gap_upper,a,h,nonradial_delta_lower))
    if gap<0 or a<=0 or a*a>2 or h<=0 or angular<0:
        raise ValueError('nonnegative proved gaps, 0<a<=sqrt(2), and h>0 required')
    penalty=n*smearing_penalty_upper(h)
    budget=gap+penalty-angular
    if budget<0:
        raise ValueError('orthogonal angular lower already exceeds the energy budget')
    target=sqrt_lower(2,scale)-a/2
    energy_error=root_upper(budget/(2*n*n*a*a),2,scale)
    smoothing_error=h*h/(4*n*a**3)
    raw=target-energy_error-smoothing_error
    return {'n':n,'scale':'Coulomb x, sum |x_i|^2=n-1','a':a,'h':h,
            'equilibrium_mean_lower':target,'radial_energy_budget_upper':budget,
            'Dirichlet_error_upper':energy_error,'disk_mean_error_upper':smoothing_error,
            'mean_lower':max(F(0),raw),'unclipped_mean_lower':raw,
            'scope':'Energy transport only. Capacity conclusions additionally require global-maximizer circle constraints.'}


def affine_layers(n,proven_gap_upper,a,h=1,nonradial_delta_lower=0,scale=10**12):
    """Necessary layers from m-1<=sqrt(8*n)*R_x, including singletons."""
    data=energy_mean_lower(n,proven_gap_upper,a,h,nonradial_delta_lower,scale)
    denominator=root_upper(8*n,2,scale)+1/F(a)
    real=n*data['mean_lower']/denominator
    return {'n':n,'a':F(a),'h':F(h),'capacity':'spectral affine circle capacity',
            'mean_certificate':data,'per_layer_weight_upper':denominator,
            'real_layer_lower':real,'minimum_layers':max(1,ceil_fraction(real)),
            'scope':'Necessary for global maxima; uses the independently proved affine same-circle capacity.'}


def pi_layers(n,proven_gap_upper,a,h=1,nonradial_delta_lower=0,scale=10**12):
    """Independent separation-only baseline, with singleton cutoff checked."""
    data=energy_mean_lower(n,proven_gap_upper,a,h,nonradial_delta_lower,scale)
    pilo,pihi=pi_interval()
    if F(a)*pilo*sqrt_lower(n,scale)<1:
        raise ValueError('cutoff not certified above 1/(pi*sqrt(n)); singleton capacity unavailable')
    denominator=pihi*root_upper(n,2,scale)
    real=n*data['mean_lower']/denominator
    return {'n':n,'a':F(a),'h':F(h),'capacity':'separation-only pi baseline',
            'mean_certificate':data,'per_layer_weight_upper':denominator,
            'real_layer_lower':real,'minimum_layers':max(1,ceil_fraction(real)),
            'scope':'Necessary for global maxima, using strict pair separation; no spectral capacity is used.'}


def cdf_mean_from_error(epsilon,a,scale=10**12):
    """Original-z target disk and certified support R<=sqrt(2) required."""
    epsilon,a=F(epsilon),F(a)
    if not 0<=epsilon<=1 or a<=0 or a*a>2:
        raise ValueError('CDF error in [0,1] and 0<a<=sqrt(2) required')
    if a*a>=2*epsilon:
        raw=sqrt_lower(2,scale)-a/2-epsilon/a+epsilon*sqrt_lower(F(1,2),scale)
    else:
        raw=sqrt_lower(2,scale)+epsilon*sqrt_lower(F(1,2),scale)-root_upper(2*epsilon,2,scale)
    return max(sqrt_lower(F(1,2),scale),raw)


def cdf_mean_lower(n,proven_gap_upper,a,scale=10**12):
    """Independent original-z CDF transport, using global R_z<sqrt(2).

    Clips the lower CDF at zero before integration, improving epsilon/a.
    This interface is slower asymptotically than direct energy transport.
    """
    require_n(n)
    gap,a=F(proven_gap_upper),F(a)
    if gap<0 or a<=0 or a*a>2:
        raise ValueError('nonnegative gap and 0<a<=sqrt(2) required')
    epsilon=cdf_error_upper(n,gap,scale)
    return {'n':n,'scale':'original z, sum |z_i|^2=n','a':a,'cdf_error_upper':epsilon,
            'mean_lower':cdf_mean_from_error(epsilon,a,scale),
            'scope':'Requires the global-maximizer original-radius cap R_z<sqrt(2); no mean constraint is used in clipping.'}


def cdf_affine_layers(n,proven_gap_upper,a,scale=10**12):
    data=cdf_mean_lower(n,proven_gap_upper,a,scale)
    denominator=root_upper(8*(n-1),2,scale)+1/F(a)
    real=n*data['mean_lower']/denominator
    return {'n':n,'a':F(a),'mean_certificate':data,
            'per_layer_weight_upper':denominator,'real_layer_lower':real,
            'minimum_layers':max(1,ceil_fraction(real)),
            'scope':'Independent CDF route, in original z scale, joined to affine global circle capacity.'}


def finite_family_best(n,proven_gap_upper,multipliers=(F(1,2),F(3,4),F(1),F(5,4),F(3,2),F(2)),
                       h_values=(F(1),F(19,20)),scale=10**12):
    """Best certificate among a prescribed finite rational cutoff/smear set.

    This maximizes only the listed certificate family, not configurations.
    """
    require_n(n)
    base=root_upper(F(1,n),4,scale)
    hs=tuple(map(F,h_values)); multipliers=tuple(map(F,multipliers))
    if not hs or not multipliers or any(h<=0 for h in hs) or any(m<=0 for m in multipliers):
        raise ValueError('nonempty positive rational multiplier families required')
    rows=[]
    for multiplier in multipliers:
        a=base*multiplier
        if a*a>2:
            continue
        for h in hs:
            rows.append(affine_layers(n,proven_gap_upper,a,h,scale=scale))
    if not rows:
        raise ValueError('no allowed cutoff in the supplied finite family')
    best=max(rows,key=lambda row:row['real_layer_lower'])
    return {'n':n,'candidate_count':len(rows),'best':best,
            'scope':'Largest exact lower certificate among the stated finite cutoff/smear choices only.'}


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('n',type=int)
    args=p.parse_args()
    gap=polar_gap_upper(args.n)
    print(json.dumps(finite_family_best(args.n,gap),default=str,indent=2))
