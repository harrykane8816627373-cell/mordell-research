"""New explicit lower bounds -> smaller global deficit -> necessary layers.

Only n=30,100 are claimed as combined examples. Logarithms are enclosed
after factoring the radial proxy and the small angular corrections; no
large rational distance product is expanded inside a logarithm.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from pathlib import Path
import json
import sys

from population_search import population_search
from radial_layers import minimum_layers, sqrt_upper, layer_delta_lower

PRIOR = Path(__file__).resolve().parents[1]/'joint-bounds'
if str(PRIOR) not in sys.path:
    sys.path.insert(0, str(PRIOR))
from subset_moment import linear_logs, log_interval, outward_decimal, disk_log_interval


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def log1p_interval(q, terms=4):
    """Short outward atanh expansion for log(1+q), -1<q<1.

    This avoids taking high powers of a large rational very close to one.
    The omitted absolute tail is bounded by a geometric series.
    """
    q = F(q)
    if not -1 < q < 1 or type(terms) is not int or terms < 1:
        raise ValueError('need -1<q<1 and positive integer series length')
    y = q/(2+q)
    partial = 2*sum((y**(2*j+1)/(2*j+1) for j in range(terms)), F(0))
    tail = 2*abs(y)**(2*terms+1)/((2*terms+1)*(1-y*y))
    return (partial, partial+tail) if y >= 0 else (partial-tail, partial)


def candidate_log_interval(result, terms=28, angular_terms=4):
    """Factored logarithm of the ACTUAL selected regular-ring product."""
    n = result['n']
    rings, signs = result['rings'], result['signs']
    N = n*(n-1)//2
    items = [(-N, F(n-1))]
    for m, b in rings:
        items.extend(((m, F(m)), (F(m*b, 2), F(b))))
    lo, hi = linear_logs(items, terms=terms)
    factor_check = result['surrogate']
    for i, j in combinations(range(len(rings)), 2):
        mi, a = rings[i]
        mj, b = rings[j]
        g = gcd(mi, mj)
        L = mi*mj//g
        relative_sign = signs[j]**(mi//g)*signs[i]**(mj//g)
        correction = -relative_sign*F(a, b)**(L//2)
        lower, upper = log1p_interval(correction, angular_terms)
        lo += 2*g*lower
        hi += 2*g*upper
        factor_check *= (1+correction)**(2*g)
    require(factor_check == result['value'], 'factored angular/radial product equals exact actual witness')
    return outward_decimal((lo, hi), 12)


def decimal_string(q, places=12):
    scale = 10**places
    value = F(q)*scale
    require(value.denominator == 1, 'decimal mesh required')
    v = value.numerator
    sign = '-' if v < 0 else ''
    v = abs(v)
    return '%s%d.%0*d' % (sign, v//scale, places, v % scale)


def run():
    independent_checks = 0
    for q in (F(-3,4), F(-1,3), F(0), F(1,5), F(3,4)):
        lo, hi = log1p_interval(q)
        a, b = log_interval(1+q)
        require(lo <= a <= b <= hi, 'independent ordinary-log enclosure inside short log1p enclosure')
        independent_checks += 1
    small = population_search(14)
    lo, hi = candidate_log_interval(small)
    a, b = log_interval(small['value'], terms=10)
    require(lo <= b and a <= hi, 'factored candidate log overlaps independent direct rational logarithm')
    rows = []
    for n in (30, 100):
        candidate = population_search(n)
        candidate_log = candidate_log_interval(candidate)
        disk_log = disk_log_interval(n)
        gap_interval = outward_decimal((disk_log[0]-candidate_log[1], disk_log[1]-candidate_log[0]), 12)
        require(0 <= gap_interval[0] <= gap_interval[1], 'nonnegative rigorous global gap enclosure')
        gap_upper = gap_interval[1]
        simple_three_layer_lower = None
        if n == 100:
            simple_three_layer_lower = layer_delta_lower(100, 3, h=1)['delta_lower']
            require(simple_three_layer_lower == F(9, 4), 'simple exact K<=3 lower bound')
            require(gap_upper < simple_three_layer_lower, 'all global 100-point optimizers have at least four radii')
        old_budget = F(n,100)+F(77,6)*sqrt_upper(n)+F(35,4)
        require(gap_upper < old_budget, 'new attainable lower strictly tightens global deficit budget')
        old_layers = minimum_layers(n, old_budget)
        new_layers = minimum_layers(n, gap_upper)
        require(new_layers['minimum_layers'] >= old_layers['minimum_layers'], 'smaller global gap cannot weaken layer count')
        if new_layers['minimum_layers'] > 1:
            require(new_layers['last_excluded_certificate']['delta_lower'] > gap_upper,
                    'last excluded layer count has exact contradiction')
        rows.append({'n': n, 'candidate_center': candidate['center'],
            'candidate_rings': candidate['rings'], 'candidate_signs': candidate['signs'],
            'log_candidate_interval': candidate_log,
            'log_U_over_candidate_interval': gap_interval,
            'new_global_optimizer_gap_upper_decimal': decimal_string(gap_upper),
            'old_effective_polar_gap_upper': old_budget,
            'old_necessary_layers': old_layers['minimum_layers'],
            'new_necessary_layers': new_layers['minimum_layers'],
            'simple_h1_three_layer_delta_lower': simple_three_layer_lower,
            'last_excluded_layer_certificate': new_layers['last_excluded_certificate']})
    return {'passed': True, 'independent_log1p_interval_checks': independent_checks,
        'independent_factored_vs_direct_log_checks': 1, 'computed_point_counts': [30, 100],
        'scope': 'Actual new attainable candidates give proven upper bounds on log(U_n/M_n), which are fed into the same universal layer theorem. These are necessary layer counts for global maximizers, not their actual layer counts or a global-optimality claim.',
        'rows': rows}


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('joint_budget_verification.json').write_text(
        json.dumps(result, default=str, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, default=str, ensure_ascii=False, indent=2))
