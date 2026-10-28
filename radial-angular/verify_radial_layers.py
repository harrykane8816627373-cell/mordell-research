"""Independent exact checks for radial_layers.py and the stated constants."""
from fractions import Fraction as F
from pathlib import Path
import json

from radial_layers import (sqrt_upper, log_bounds, smearing_penalty_upper,
    layer_delta_lower, radial_cover_delta_lower, minimum_layers, free_gap_polynomial_lower,
    free_gap_energy_bounds, population_gap_energy_bounds,
    upper_support_le, lower_support_ge, radial_gap_certificate)


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def independent_log(q, terms=80):
    """Direct atanh series, independently without power-of-two reduction."""
    q = F(q)
    if q < 1:
        a, b = independent_log(1/q, terms)
        return -b, -a
    t = (q-1)/(q+1)
    lo = 2*sum((t**(2*j+1)/(2*j+1) for j in range(terms)), F(0))
    return lo, lo+2*t**(2*terms+1)/((2*terms+1)*(1-t*t))


def run():
    root_checks = 0
    for a in range(40):
        for b in range(1, 12):
            value = sqrt_upper(F(a, b), 10000)
            require(value*value >= F(a, b), 'outward square root')
            if value:
                require((value-F(1,10000))**2 < F(a,b), 'minimal root grid endpoint')
            root_checks += 1
    logarithm_checks = 0
    for q in (F(1), F(99,100), F(100,99), F(3,2), F(2), F(7,4), F(4), F(1,8)):
        a, b = log_bounds(q)
        c, d = independent_log(q)
        require(a <= d and c <= b, 'independent log intervals overlap')
        logarithm_checks += 1
    # Exact derivative identity for log(v)-3(v^2-1)/(v^2+4v+1).
    derivative_checks = 0
    for v in (F(1,7), F(1), F(6,5), F(2), F(10)):
        denominator = v*v+4*v+1
        numerator = denominator**2-12*v*(v*v+v+1)
        require(numerator == (v-1)**4, 'gap logarithm derivative numerator')
        derivative_checks += 1
    gap_checks = 0
    for a, b in ((F(0), F(1)), (F(1,100), F(1)), (F(1,5), F(6,5)),
                 (F(1,2), F(3,4)), (F(1), F(7,5))):
        low, high = free_gap_energy_bounds(a,b)
        require(low <= high, 'free gap interval ordering')
        require(low >= free_gap_polynomial_lower(a,b), 'polynomial lower bound retained')
        for n, population in ((8,0),(8,1),(8,4),(8,8)):
            if a == 0 and population:
                continue
            lo, hi = population_gap_energy_bounds(n,population,a,b)
            require(lo <= hi and hi >= low, 'fixed population versus freely minimized gap')
            if a:
                l,h = independent_log(b/a)
                c = F(population,n)
                constant = -c*(b*b-a*a)/2+(b**4-a**4)/16
                require(lo <= c*c*h+constant and hi >= c*c*l+constant,
                        'independent fixed population log integration')
            gap_checks += 1
    support_checks = 0
    for center in (F(0), F(1,3), F(1), F(3,2)):
        for r in (F(1,10), F(1,2), F(1)):
            for endpoint in (F(0), F(1,4), F(1), F(2)):
                require(upper_support_le(center**2,r**2,endpoint**2) ==
                        (center+r <= endpoint), 'upper support radical elimination')
                require(lower_support_ge(center**2,r**2,endpoint**2) ==
                        (max(center-r,F(0)) >= endpoint), 'lower support radical elimination')
                support_checks += 2
    layer_monotonicity_checks = 0
    for n in (2,8,30,100,1000):
        for h in (F(1),F(99,100),F(4,5),F(6,5)):
            values = [layer_delta_lower(n,k,h)['delta_lower'] for k in range(1,min(n,40)+1)]
            require(all(a >= b for a,b in zip(values,values[1:])), 'layer bound decreases with K')
            layer_monotonicity_checks += len(values)-1
    for n in (100,1000,1000000):
        base = layer_delta_lower(n,2)
        require(radial_cover_delta_lower(n,2,0)['delta_lower'] == base['delta_lower'],
                'zero-width bands equal exact layers')
        require(radial_cover_delta_lower(n,2,F(1,10))['delta_lower'] <= base['delta_lower'],
                'positive band width weakens the necessary energy cost')
    # Finite n: K >= 2 and K <= (3/22)sqrt(n) force the contradictory constants.
    coefficient = F(3,22)
    finite_lower = 9*(1-2*coefficient)**4/(64*(F(3,2)*coefficient)**2)
    finite_upper = F(1,100)+F(77,6)*coefficient/2+F(35,4)*coefficient**2/4
    require(finite_lower == F(1024,1089), 'finite layer lower constant')
    require(finite_upper == F(179211,193600), 'finite polar upper constant')
    require(finite_lower > finite_upper, 'all-n 3/22 contradiction')
    # No equal-radius instability input needed for the weaker universal 1/10 version.
    require(F(36,25) > F(1,100)+F(77,60)+F(35,400), 'all-n 1/10 contradiction')
    # Optimized smear h=99/100 raises the asymptotic layer guarantee above .4.
    h, coefficient = F(99,100), F(101,250)
    asymptotic_lower = 9*(1-2*h*coefficient)**4/(64*coefficient**2)-smearing_penalty_upper(h)
    require(asymptotic_lower > F(1,800), 'asymptotic 101/250 contradiction')
    require(smearing_penalty_upper(F(1)) == 0, 'default smear has no penalty')
    # Entire thick radial class, exact integral, no numerical circle areas.
    thick = radial_gap_certificate(100,[(F(49,50),F(1),100)],[(F(0),F(4,5))])
    require(thick['delta_lower'] == 256, 'thick-ring exact radius-gap certificate')
    require(thick['delta_lower'] > F(1657,12), 'thick ring excludes all global maximizers')
    population = radial_gap_certificate(100,[(F(1,4),F(1,4),50),
        (F(173,100),F(173,100),50)],[(F(7,10),F(11,10))])
    require(population['gaps'][0]['cumulative_population'] == 50, 'actual gap population')
    require(population['delta_lower'] > 0, 'population-weighted gap cost positive')
    invalid_cases = [
        lambda: layer_delta_lower(1,1), lambda: layer_delta_lower(8,0),
        lambda: layer_delta_lower(8,2,F(0)), lambda: free_gap_energy_bounds(1,2),
        lambda: population_gap_energy_bounds(8,1,0,F(1,2)),
        lambda: radial_gap_certificate(8,[(0,0,1),(1,1,7)],[(F(1,3),F(2,3))]),
        lambda: radial_gap_certificate(8,[(0,0,1),(1,1,7)],[(0,F(1,4)),(F(1,5),F(1,3))]),
        lambda: minimum_layers(8,-1)]
    for call in invalid_cases:
        try:
            call()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid input accepted')
    examples = []
    for n in (100,1000,10000,1000000,100000000):
        gap = F(n,100)+F(77,6)*sqrt_upper(n)+F(35,4)
        row = minimum_layers(n,gap)
        k = row['minimum_layers']
        require(1 <= k <= n, 'necessary layer population range')
        if k > 1:
            require(row['last_excluded_certificate']['delta_lower'] > gap,
                    'binary search last excluded level')
        examples.append(row)
    return {'passed': True, 'arithmetic': 'integer/Fraction and outward logarithm intervals',
        'outward_root_checks': root_checks, 'logarithm_checks': logarithm_checks,
        'derivative_identity_checks': derivative_checks, 'gap_integral_checks': gap_checks,
        'support_comparison_checks': support_checks,
        'layer_monotonicity_checks': layer_monotonicity_checks,
        'rejected_invalid_cases': len(invalid_cases),
        'finite_layer_constant': '3/22', 'finite_lower_comparison': finite_lower,
        'finite_upper_comparison': finite_upper, 'asymptotic_layer_constant': '101/250',
        'asymptotic_contradiction_lower': asymptotic_lower,
        'thick_ring_delta_lower': thick['delta_lower'],
        'population_gap_certificate': population, 'effective_examples': examples,
        'scope': 'Finite arithmetic checks; the all-n and asymptotic claims are proved in RADIAL_LAYERS.md.'}


if __name__ == '__main__':
    report = run()
    destination = Path(__file__).with_name('radial_layers_verification.json')
    destination.write_text(json.dumps(report,default=str,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,indent=2))
