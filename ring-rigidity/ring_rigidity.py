"""Necessary full-stationarity and Hessian conditions, never sufficiency tests."""
from collections import Counter
from fractions import Fraction as F
from math import comb


def require(ok, message):
    if not ok:
        raise ValueError(message)


def population_case(a, b, c):
    require(all(type(x) is int and x >= 2 for x in (a, b, c)), 'populations >=2 required')
    a, b, c = sorted((a, b, c))
    if a == c:
        return 'equal_population_not_excluded'
    if b == c:
        return 'two_largest_equal'
    if a + b != c:
        return 'unique_largest_unbalanced'
    if a != b:
        return 'unique_largest_balanced_unequal'
    return 'unique_largest_balanced_equal'


def fourier_coefficients(m):
    require(type(m) is int and m >= 2, 'rotation order >=2 required')
    # k=m denotes the constant Fourier mode, not a new mode k=0.
    return [F((k-2)*(k-m), 2) for k in range(1, m+1)]


def rotational_bound(m, orbits, center=0):
    require(type(m) is int and m >= 2, 'rotation order >=2 required')
    require(type(orbits) is int and orbits >= 1, 'positive orbit count required')
    require(type(center) is int and center in (0, 1), 'center in {0,1} required')
    n = m*orbits + center
    floor = m*m//4
    rhs = 2*(n-1+center)
    ratio = F(floor, rhs)
    return {'n':n, 'rotation_order':m, 'nonzero_orbits':orbits, 'center':center,
            'floor_m_squared_over_four':floor, 'necessary_upper_bound':rhs,
            'operator_norm_over_kappa_lower':str(ratio),
            'full_hessian_instability_proved':ratio > 1,
            'hessian_largest_eigenvalue_over_kappa_lower':str(2*(ratio-1)) if ratio > 1 else None,
            'necessary_condition_only':True}


def hermite_polynomial(n):
    """Monic physicists' Hermite polynomial; exact rational coefficients."""
    require(type(n) is int and n >= 0, 'nonnegative degree required')
    if n == 0:
        return [F(1)]
    previous, current = [F(1)], [F(0), F(1)]
    for k in range(1, n):
        nxt = [F(0)] + current
        for j, coefficient in enumerate(previous):
            nxt[j] -= F(k, 2)*coefficient
        previous, current = current, nxt
    return current


def certificate():
    census = Counter()
    for a in range(2, 101):
        for b in range(a, 101):
            for c in range(b, 101-a-b):
                census[population_case(a,b,c)] += 1
    profiles = [(2,2,2), (2,2,3), (2,2,4), (2,3,5), (3,4,4),
                (3,4,5), (3,3,6), (4,4,4), (8,8,16), (9,10,10)]
    boundaries = []
    for orbits in range(1, 13):
        for center in (0,1):
            threshold = next(m for m in range(2, 8*orbits+4)
                             if rotational_bound(m,orbits,center)['full_hessian_instability_proved'])
            boundaries.append({'orbits':orbits, 'center':center, 'first_excluded_order':threshold,
                               'at_threshold':rotational_bound(threshold,orbits,center),
                               'immediately_below':rotational_bound(threshold-1,orbits,center)})
    return {
        'date':'2026-11-05',
        'scope':'unit-strength unrestricted pointwise stationarity; no global optimum claimed',
        'three_distinct_regular_rings_require_equal_populations':True,
        'population_census':{'minimum_population':2, 'maximum_total_noncentral_population':100,
                             'counts':dict(sorted(census.items())), 'profiles_checked':sum(census.values()),
                             'finite_census_replaces_general_proof':False},
        'examples':[{'populations':list(p), 'case':population_case(*p)} for p in profiles],
        'rotational_order_necessary_condition':'floor(m^2/4) <= 2*(n-1+h), n=K*m+h, h in {0,1}',
        'three_ring_thresholds':{'without_center':24, 'with_center':25},
        'orbit_thresholds':boundaries,
        'stationary_but_unstable_examples':[
            {'n':n,'center':n%2,'populations':[2,2,2],
             'polynomial_in_increasing_degree_order':list(map(str,hermite_polynomial(n))),
             'squared_radius_isolation_intervals':intervals,
             'force_multiplier':'1',
             'full_hessian_direction':'v_i = i*(z_i^2-(n-1)/2)',
             'direction_norm_squared':str(F(n*(n-1)*(n-2),4)),
             'positive_second_variation':str(F(n*(n-1)*(n-2),2))}
            for n,intervals in [(6,[['0','1/4'],['1/4','2'],['2','8']]),
                                (7,[['0','1'],['1','3'],['3','8']])]
        ]}


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(),indent=2))
