"""Coupled radial-bin quantization, with exact convex dual certificates.

Proof: RADIAL_QUANTIZATION.md. The bins use u=|x|/sqrt(2), where Q(x)=n-1.
Layer/population inputs are necessary-condition hypotheses, not constructions.
"""
from bisect import bisect_right
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'radial-angular'))
from radial_layers import smearing_penalty_upper, sqrt_upper


def require_n(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n >= 2 required')


def bin_cost(left, width, h, y):
    """Convex decreasing f(y)=a*(w-2*h*y)_+^3/(3*y^2), for y>0."""
    a, w, h, y = map(F, (left,width,h,y))
    if a < 0 or w <= 0 or h <= 0 or y <= 0:
        raise ValueError('a>=0, w>0, h>0, y>0 required')
    return a*max(F(0),w-2*h*y)**3/(3*y*y)


def bin_derivative(left, width, h, y):
    a, w, h, y = map(F, (left,width,h,y))
    bin_cost(a,w,h,y)
    if w <= 2*h*y:
        return F(0)
    return a*(-2*w**3/y**3+6*h*w*w/y**2-8*h**3)/3


def cubic_grid_lower(target, scale=10**8):
    """Largest k/scale with z^2*(3+2*z)<=target, using integers only."""
    target = F(target)
    if target < 0 or type(scale) is not int or scale < 1:
        raise ValueError('nonnegative rational and positive integer scale required')
    def valid(k):
        return (3*k*k*scale+2*k**3)*target.denominator <= target.numerator*scale**3
    lo, hi = 0, scale
    while valid(hi):
        hi *= 2
    while lo+1 < hi:
        mid = (lo+hi)//2
        if valid(mid):
            lo = mid
        else:
            hi = mid
    return F(lo,scale)


def dual_constants(bins=128, h=F(197,200), slope=F(7,100), root_scale=10**8):
    """Construct one globally valid dual line B-lambda*C for bin budgets.

    All nonzero bins use a rational tangent with slope no smaller than
    -lambda. No optimizer or unverified stationary equation is involved.
    """
    h, lam = F(h), F(slope)
    if type(bins) is not int or bins < 2 or h <= 0 or lam <= 0:
        raise ValueError('integer bins>=2, positive h and positive slope required')
    width, intercept, anchors = F(1,bins), F(0), []
    for j in range(1,bins):
        left = F(j,bins)
        target = 3*lam/(left*(2*h)**3)
        z = cubic_grid_lower(target,root_scale)
        tangent_slope = left*(2*h)**3*z*z*(3+2*z)/3
        b = left*width*(2*h)**2*z*z
        if tangent_slope > lam:
            raise ArithmeticError('tangent slope outside dual budget')
        intercept += b
        anchors.append({'bin':j, 'z_lower':z})
    return {'bins':bins, 'active_bins':bins-1, 'h':h, 'slope':lam,
            'intercept':intercept, 'root_scale':root_scale, 'anchors':anchors,
            'penalty_upper':smearing_penalty_upper(h),
            'scope':'Exact dual line for all bin layer allocations; the first bin is safely omitted.'}


def fixed_layer_certificate():
    result = dual_constants()
    short_penalty = F(1,4400)
    if result['penalty_upper'] > short_penalty:
        raise ArithmeticError('short smearing penalty not certified')
    result['penalty_upper'] = short_penalty
    result['asymptotic_gap_upper'] = F(1,800)
    result['asymptotic_coefficient_lower'] = (
        result['intercept']-short_penalty-F(1,800))/result['slope']
    result['headline_coefficient'] = F(13,30)
    if result['asymptotic_coefficient_lower'] <= F(13,30):
        raise ArithmeticError('13/30 asymptotic coefficient not certified')
    return result


def uniform_delta_lower(n, layers, constants=None, scale=10**12):
    """Necessary energy cost for every shape with at most `layers` radii."""
    require_n(n)
    if type(layers) is not int or not 1 <= layers <= n:
        raise ValueError('integer layer count in [1,n] required')
    data = fixed_layer_certificate() if constants is None else constants
    inverse_root = sqrt_upper(F(1,n),scale)
    rate = data['intercept']-data['penalty_upper']-data['slope']*(layers+data['active_bins'])*inverse_root
    return {'n':n, 'layers_upper':layers, 'delta_lower':n*max(F(0),rate),
            'delta_per_n_lower':rate, 'bins':data['bins'], 'h':data['h'],
            'scope':'At most the stated number of distinct radii; no angular or stationarity assumption.'}


def necessary_layers(n, proven_gap_upper, constants=None, scale=10**12):
    """Integer necessary layer count from a separately PROVED delta upper."""
    require_n(n)
    gap = F(proven_gap_upper)
    if gap < 0 or type(scale) is not int or scale < 1:
        raise ValueError('nonnegative proved energy gap and positive root scale required')
    data = fixed_layer_certificate() if constants is None else constants
    coefficient = (data['intercept']-data['penalty_upper']-gap/n)/data['slope']
    root_lower = F(isqrt(n*scale*scale),scale)
    real_lower = root_lower*max(F(0),coefficient)-data['active_bins']
    integer_lower = max(1,-((-real_lower.numerator)//real_lower.denominator))
    return {'n':n, 'proven_gap_upper':gap, 'minimum_layers':integer_lower,
            'real_layer_lower':real_lower,
            'scope':'Necessary only. Caller must independently prove the supplied gap upper bound.'}


def allocation_delta_lower(n, local_layer_upper, h=F(197,200), scale=10**12):
    """Conditional bound from per-bin upper counts of distinct center radii.

    Bin j is j/J <= |x|/sqrt(2) < (j+1)/J, with the outer endpoint 1
    included in the last bin. Radii above 1 may occur and are ignored.
    Point-population upper counts are also safe upper counts for layers.
    """
    require_n(n)
    counts = list(local_layer_upper)
    if len(counts) < 2 or any(type(k) is not int or k < 0 for k in counts) or sum(counts) > n:
        raise ValueError('at least two nonnegative integer bin counts with sum <= n required')
    h = F(h)
    if h <= 0:
        raise ValueError('positive smear factor required')
    bins, total, rows = len(counts), F(0), []
    inverse_root = sqrt_upper(F(1,n),scale)
    for j in range(1,bins):
        y = (counts[j]+1)*inverse_root
        value = bin_cost(F(j,bins),F(1,bins),h,y)
        total += value
        rows.append({'bin':j, 'layer_upper':counts[j], 'cost_per_n_lower':value})
    penalty = smearing_penalty_upper(h)
    return {'n':n, 'bins':bins, 'h':h, 'local_counts':counts,
            'radial_cost_per_n_lower':total, 'penalty_per_n_upper':penalty,
            'delta_lower':n*max(F(0),total-penalty), 'bin_costs':rows,
            'scope':'Conditional on all per-bin layer/population upper counts; no existence claimed.'}


def exact_radial_census(squared_radii, bins=128):
    """Exact census using only rational q=|x|^2 comparisons, Q(x)=n-1."""
    qs = list(map(F,squared_radii))
    n = len(qs)
    require_n(n)
    if type(bins) is not int or bins < 2 or any(q < 0 for q in qs) or sum(qs) != n-1:
        raise ValueError('nonnegative Coulomb squared radii summing to n-1 and bins>=2 required')
    edges = [2*F(j,bins)**2 for j in range(bins+1)]
    population, layers, outside = [0]*bins, [set() for _ in range(bins)], 0
    for q in qs:
        if q > 2:
            outside += 1
            continue
        j = min(bins-1,bisect_right(edges,q)-1)
        population[j] += 1
        layers[j].add(q)
    return {'n':n, 'bins':bins, 'population':population,
            'local_layers':[len(v) for v in layers], 'outside_population':outside,
            'total_distinct_radii':len(set(qs)),
            'scope':'Exact radial census only; no angles or optimality are inferred.'}


if __name__ == '__main__':
    import json
    data = fixed_layer_certificate()
    print(json.dumps({k:v for k,v in data.items() if k != 'anchors'},default=str,indent=2))
