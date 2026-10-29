"""Exact independent arithmetic checks for coupled radial quantization."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

from radial_quantization import (bin_cost,bin_derivative,cubic_grid_lower,
    dual_constants,fixed_layer_certificate,uniform_delta_lower,necessary_layers,
    allocation_delta_lower,exact_radial_census)


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def explicit_polynomial_cost(a,w,h,y):
    if w <= 2*h*y:
        return F(0)
    return a*(w**3/y**2-6*h*w*w/y+12*h*h*w-8*h**3*y)/3


def check_core_geometry(bins,epsilon,centers):
    """Exact interval union checks, independent of the production dual code."""
    checks = 0
    for j in range(1,bins):
        a,b = F(j,bins),F(j+1,bins)
        left,right = a+epsilon,b-epsilon
        if left >= right:
            continue
        local = [u for u in centers if a <= u < b or (j == bins-1 and u == b)]
        for u in centers:
            if u not in local:
                require(u+epsilon <= left or u-epsilon >= right,
                        'centers outside bin cannot enter its open core')
                checks += 1
        intervals = sorted((max(left,u-epsilon),min(right,u+epsilon)) for u in local
                           if max(left,u-epsilon) < min(right,u+epsilon))
        merged = []
        for x,y in intervals:
            if merged and x <= merged[-1][1]:
                merged[-1] = (merged[-1][0],max(merged[-1][1],y))
            else:
                merged.append((x,y))
        coverage = sum((y-x for x,y in merged),F(0))
        require(coverage <= 2*epsilon*len(local),'local support coverage budget')
        require(len(merged)+1 <= len(local)+1,'local empty-gap count budget')
        checks += 2
    return checks


def run():
    root_checks = 0
    for numerator in range(31):
        for denominator in range(1,8):
            target = F(numerator,denominator)
            z = cubic_grid_lower(target,1000)
            require(z*z*(3+2*z) <= target,'cubic root lower grid endpoint')
            v = z+F(1,1000)
            require(v*v*(3+2*v) > target,'maximal cubic root lower grid endpoint')
            root_checks += 1
    tangent_checks = 0
    for a,w,h,z in product((F(1,8),F(1,2),F(7,8)),(F(1,8),F(1,3)),
                          (F(197,200),F(1),F(6,5)),(F(0),F(1,10),F(1),F(2))):
        y0 = w/(2*h*(1+z))
        slope = a*(2*h)**3*z*z*(3+2*z)/3
        intercept = a*w*(2*h)**2*z*z
        require(bin_derivative(a,w,h,y0) == -slope,'tangent derivative identity')
        require(bin_cost(a,w,h,y0)+slope*y0 == intercept,'tangent intercept identity')
        for y in (y0/3,y0/2,y0,2*y0,w/(2*h),w/h):
            value = bin_cost(a,w,h,y)
            require(value == explicit_polynomial_cost(a,w,h,y),'independent bin expansion')
            require(value >= intercept-slope*y,'global convex supporting line')
            tangent_checks += 1
    geometry_checks = 0
    centers = sorted(set([F(j,16) for j in range(20)]+[F(1,3),F(3,4),F(1)]))
    for bins in (2,4,8):
        for epsilon in (F(1,100),F(197,20000),F(1,20)):
            geometry_checks += check_core_geometry(bins,epsilon,centers)
    allocation_checks = 0
    small = dual_constants(bins=4,h=F(1),slope=F(1,10),root_scale=10**6)
    for counts in product(range(4),repeat=4):
        n = 10000
        exact = allocation_delta_lower(n,counts,h=1)
        dual = uniform_delta_lower(n,max(1,sum(counts)),small)
        # An all-zero census would force all points outside the target bins;
        # the algebraic comparison remains valid, but no existence is claimed.
        require(exact['delta_lower'] >= dual['delta_lower'],'coupled dual below every local allocation')
        allocation_checks += 1
    constants = fixed_layer_certificate()
    B,lam,h = constants['intercept'],constants['slope'],constants['h']
    independently_summed = F(0)
    for row in constants['anchors']:
        j,z = row['bin'],row['z_lower']
        a,w = F(j,128),F(1,128)
        target = 3*lam/(a*(2*h)**3)
        require(z*z*(3+2*z) <= target,'fixed certificate valid tangent slope')
        require((z+F(1,10**8))**2*(3+2*(z+F(1,10**8))) > target,
                'fixed certificate maximal root grid')
        independently_summed += a*w*(2*h)**2*z*z
    require(B == independently_summed,'fixed certificate intercept sum')
    require(B == F(52147603526844565517847,1638400000000000000000000),'fixed B bytes-independent value')
    y = F(3,397)
    log_upper = 2*sum((y**(2*j+1)/(2*j+1) for j in range(4)),F(0))
    log_upper += 2*y**9/(9*(1-y*y))
    require((h*h-1)/2+log_upper < F(1,4400),'independent short smear penalty')
    margin = B-lam*F(13,30)-F(1,4400)-F(1,800)
    require(margin > 0,'13/30 asymptotic contradiction margin')
    require(margin == F(960516385870662088951,54067200000000000000000000),
            'exact positive asymptotic margin')
    # Per-bin point populations safely upper-bound local distinct radii.
    qs = [F(9999,10000)]*10000
    census = exact_radial_census(qs,bins=8)
    layer_value = allocation_delta_lower(10000,census['local_layers'])
    population_value = allocation_delta_lower(10000,census['population'])
    require(layer_value['delta_lower'] >= population_value['delta_lower'],
            'actual layers strengthen population upper counts')
    require(population_value['delta_lower'] > 0,'empty-bin population costs positive')
    require(census['total_distinct_radii'] == 1 and sum(census['population']) == 10000,
            'exact normalized radial census')
    # Rational square boundaries require no numerical roots.
    boundary = exact_radial_census([F(0),F(1,2),F(3,2)],bins=2)
    require(boundary['population'] == [1,2],'half-open radial bin boundaries')
    invalid = [lambda:bin_cost(-1,1,1,1),lambda:cubic_grid_lower(-1),
               lambda:dual_constants(1),lambda:uniform_delta_lower(8,0),
               lambda:necessary_layers(8,-1),lambda:allocation_delta_lower(8,[1,-1]),
               lambda:exact_radial_census([1,1]),lambda:allocation_delta_lower(8,[8,1])]
    for call in invalid:
        try:
            call()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid input accepted')
    return {'passed':True,'arithmetic':'integer/Fraction only; no floating optimizer',
        'cubic_root_checks':root_checks,'convex_tangent_checks':tangent_checks,
        'bin_core_geometry_checks':geometry_checks,'local_allocation_checks':allocation_checks,
        'rejected_invalid_cases':len(invalid),'certificate':constants,
        'positive_asymptotic_margin':margin,
        'conditional_population_example':{'n':10000,'bins':8,
            'local_population':census['population'],'local_layers':census['local_layers'],
            'population_delta_lower':population_value['delta_lower'],
            'layer_delta_lower':layer_value['delta_lower']},
        'scope':'Finite arithmetic checks. RADIAL_QUANTIZATION.md proves all-n bin inequalities and the asymptotic layer theorem; no optimum configuration is classified.'}


if __name__ == '__main__':
    report = run()
    destination = Path(__file__).with_name('radial_quantization_verification.json')
    destination.write_text(json.dumps(report,default=str,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'certificate'},default=str,indent=2))
