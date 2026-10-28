"""Exact radial empty-gap energies with integer population constraints.

Analytic statement: INTEGER_GAPS.md.  Population assumptions are necessary
conditions on the supplied class, not an assertion that the class exists.
Independent per-gap minimization is safe but does not enumerate all global
population sequences or all configurations.
"""
from fractions import Fraction as F
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'joint-bounds'))
from subset_moment import linear_logs, log_interval, outward_decimal


def _require_integer(value, lower, name):
    if type(value) is not int or value < lower:
        raise ValueError(name + ' must be an integer >= ' + str(lower))


def _validate_gap(a, b):
    a, b = F(a), F(b)
    if not 0 <= a < b <= 2:
        raise ValueError('need a squared-radius gap 0<=a<b<=2')
    return a, b


def gap_energy_interval(a, b, mass, terms=28):
    """Enclose 1/2 integral_a^b (mass-s/2)^2 ds/s.

    The caller must prove the smeared radial CDF is CONSTANT on (a,b).
    The endpoints are squared radii in the Coulomb x scale.
    """
    a, b = _validate_gap(a, b)
    mass = F(mass)
    if not 0 <= mass <= 1:
        raise ValueError('a probability mass in [0,1] is required')
    if a == 0:
        if mass:
            raise ValueError('nonzero constant CDF down to radius zero has infinite energy')
        return b*b/16, b*b/16
    constant = (b*b-a*a)/16-mass*(b-a)/2
    return linear_logs(((mass*mass/2, b/a),), constant, terms)


def free_gap_energy_interval(a, b, terms=28):
    """Free-real-constant minimum, for comparison only (not integer aware)."""
    a, b = _validate_gap(a, b)
    if a == 0:
        return b*b/16, b*b/16
    lo, hi = log_interval(b/a, terms)
    constant, coefficient = (b*b-a*a)/16, (b-a)**2/8
    return constant-coefficient/lo, constant-coefficient/hi


def integer_gap_interval(n, a, b, count_lower=0, count_upper=None,
                         terms=28, candidate_budget=32):
    """Minimize the exact gap energy over stated integer cumulative counts.

    Only O(1) candidates are normally needed.  An outward log interval
    brackets the continuous minimizer; convexity excludes all other counts.
    A precision/budget failure rejects the calculation rather than omitting
    potentially optimal counts.  Counts are normalized by n.
    """
    _require_integer(n, 2, 'n')
    _require_integer(count_lower, 0, 'count_lower')
    _require_integer(candidate_budget, 1, 'candidate_budget')
    if count_upper is None:
        count_upper = n
    _require_integer(count_upper, 0, 'count_upper')
    if not count_lower <= count_upper <= n:
        raise ValueError('need 0<=count_lower<=count_upper<=n')
    a, b = _validate_gap(a, b)
    if a == 0:
        if count_lower:
            raise ValueError('no finite-energy CDF with positive count at radius zero')
        return {'energy_interval': gap_energy_interval(a, b, 0),
                'candidate_counts': [0], 'count_range': [count_lower, count_upper]}
    loglo, loghi = log_interval(b/a, terms)
    optimum_lo = n*(b-a)/(2*loghi)
    optimum_hi = n*(b-a)/(2*loglo)
    lo_floor = optimum_lo.numerator // optimum_lo.denominator
    hi_ceil = -((-optimum_hi.numerator) // optimum_hi.denominator)
    lower = min(count_upper, max(count_lower, lo_floor))
    upper = min(count_upper, max(count_lower, hi_ceil))
    if upper-lower+1 > candidate_budget:
        raise ValueError('log precision insufficient for the integer-candidate budget')
    candidates = list(range(lower, upper+1))
    intervals = [gap_energy_interval(a, b, F(j, n), terms) for j in candidates]
    interval = min(lo for lo, hi in intervals), min(hi for lo, hi in intervals)
    return {'energy_interval': interval, 'candidate_counts': candidates,
            'count_range': [count_lower, count_upper],
            'continuous_count_interval': (optimum_lo, optimum_hi)}


def upper_support_le(q, r2, boundary):
    """(sqrt(q)+sqrt(r2))^2 <= boundary, without floating roots."""
    q, r2, boundary = map(F, (q, r2, boundary))
    if min(q, r2, boundary) < 0:
        raise ValueError('nonnegative squared quantities required')
    c = boundary-q-r2
    return c >= 0 and c*c >= 4*q*r2


def lower_support_ge(q, r2, boundary):
    """max(sqrt(q)-sqrt(r2),0)^2 >= boundary, without floating roots."""
    q, r2, boundary = map(F, (q, r2, boundary))
    if min(q, r2, boundary) < 0:
        raise ValueError('nonnegative squared quantities required')
    if boundary == 0:
        return True
    c = q-r2-boundary
    return c >= 0 and c*c >= 4*r2*boundary


def band_gap_population(n, bands, a, b):
    """Prove an empty smeared gap and its integer population from q-bands.

    bands=(q_lower,q_upper,population), in x scale with sum|x_i|²=n-1.
    A band which might touch the interior makes this certificate unavailable.
    """
    _require_integer(n, 2, 'n')
    a, b = _validate_gap(a, b)
    bands = [(F(lo), F(hi), population) for lo, hi, population in bands]
    if any(not 0 <= lo <= hi or type(population) is not int or population < 1
           for lo, hi, population in bands):
        raise ValueError('valid bands and positive integer populations required')
    if sum(population for lo, hi, population in bands) != n:
        raise ValueError('band populations must sum to n')
    if not sum(lo*j for lo, hi, j in bands) <= n-1 <= sum(hi*j for lo, hi, j in bands):
        raise ValueError('bands incompatible with the prescribed second moment')
    r2, count = F(2, n), 0
    for lo, hi, population in bands:
        if upper_support_le(hi, r2, a):
            count += population
        elif not lower_support_ge(lo, r2, b):
            raise ValueError('a band may place smeared mass inside the proposed gap')
    return count


def band_gap_certificate(n, bands, gaps, terms=28):
    """Add exact energies of disjoint empty gaps for a whole radial class.

    The returned delta lower bound can be added to an independently proved
    angular projection penalty, because the radial/angular energies are
    orthogonal.  Two different radial lower bounds generally take a maximum.
    """
    gaps = [tuple(map(F, gap)) for gap in gaps]
    if not gaps:
        raise ValueError('at least one gap is required')
    if any(a >= b or a < 0 or b > 2 for a, b in gaps):
        raise ValueError('valid squared-radius gaps in [0,2] required')
    if any(gaps[i][1] > gaps[i+1][0] for i in range(len(gaps)-1)):
        raise ValueError('ordered gaps with disjoint interiors required')
    lo_total = hi_total = F(0)
    previous_count, records = 0, []
    for a, b in gaps:
        count = band_gap_population(n, bands, a, b)
        if count < previous_count:
            raise ArithmeticError('CDF population decreased across ordered gaps')
        previous_count = count
        lo, hi = gap_energy_interval(a, b, F(count, n), terms)
        lo_total += lo
        hi_total += hi
        records.append({'squared_radius_gap': [str(a), str(b)],
                        'cdf_count': count, 'cdf_mass': str(F(count, n)),
                        'energy_interval': [str(x) for x in outward_decimal((lo, hi))]})
    return {'n': n, 'radial_energy_interval': (lo_total, hi_total),
            'delta_lower': n*n*lo_total, 'gaps': records,
            'scope': 'Conditional class certificate; exact populations on disjoint empty smeared radial gaps. No global domain covering is asserted.'}


if __name__ == '__main__':
    import json
    example = band_gap_certificate(100,
        [(F(1,25), F(1,20), 10), (F(109,100), F(11,10), 90)],
        [(F(1,5), F(4,5))])
    example['radial_energy_interval'] = [str(x) for x in outward_decimal(example['radial_energy_interval'])]
    example['delta_lower'] = str(outward_decimal((example['delta_lower'], example['delta_lower']))[0])
    print(json.dumps(example, indent=2))
