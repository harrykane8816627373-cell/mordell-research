"""Sharp two-moment logarithmic envelope for the existing subset recurrence.

Only standard-library exact arithmetic is used.  The analytic argument and
the distinction between a stronger recurrence and the global upper envelope
are in SUBSET_MOMENT.md.  No subset or configuration enumeration is implied.
"""
from fractions import Fraction as F
from functools import lru_cache


def _integer(value, lower, name):
    if type(value) is not int or value < lower:
        raise ValueError(name + ' must be an integer >= ' + str(lower))


@lru_cache(maxsize=1024)
def log_interval(value, terms=28):
    """Return rational lower/upper bounds for log(value), value > 0."""
    value = F(value)
    _integer(terms, 1, 'terms')
    if value <= 0:
        raise ValueError('positive logarithm input required')
    if value < 1:
        lo, hi = log_interval(1 / value, terms)
        return -hi, -lo
    if value == 1:
        return F(0), F(0)
    exponent = max(0, value.numerator.bit_length()
                   - value.denominator.bit_length())
    reduced = value / (1 << exponent)
    if reduced < 1:
        exponent -= 1
        reduced *= 2
    if reduced >= 2:
        exponent += 1
        reduced /= 2

    def core(x):
        y = (x - 1) / (x + 1)
        square, power, low = y * y, y, F(0)
        for j in range(terms):
            low += 2 * power / (2 * j + 1)
            power *= square
        tail = 2 * power / ((2 * terms + 1) * (1 - square))
        return low, low + tail

    low, high = core(reduced)
    low_two, high_two = core(F(2))
    return low + exponent * low_two, high + exponent * high_two


def linear_logs(items, constant=0, terms=28):
    """Enclose constant + sum(coefficient * log(argument))."""
    low = high = F(constant)
    for coefficient, argument in items:
        coefficient = F(coefficient)
        lo, hi = log_interval(F(argument), terms)
        if coefficient >= 0:
            low += coefficient * lo
            high += coefficient * hi
        else:
            low += coefficient * hi
            high += coefficient * lo
    return low, high


def outward_decimal(interval, places=12):
    """Round a rational interval outward to a decimal rational mesh."""
    _integer(places, 0, 'places')
    low, high = interval
    scale = 10 ** places
    scaled_low, scaled_high = low * scale, high * scale
    return (F(scaled_low.numerator // scaled_low.denominator, scale),
            F(-((-scaled_high.numerator) // scaled_high.denominator), scale))


def two_moment_parameters(mu, upper, variance_lower):
    """Return the extremal support a and probability p at a.

    Requires 0 < v < mu*(upper-mu).  The zero-variance and degenerate
    support cases should instead use the ordinary log(mu) bound.
    """
    mu, upper, v = map(F, (mu, upper, variance_lower))
    if not (0 < mu < upper and 0 < v < mu * (upper - mu)):
        raise ValueError('need 0 < mu < L and 0 < v < mu*(L-mu)')
    a = mu - v / (upper - mu)
    p = (upper - mu) ** 2 / ((upper - mu) ** 2 + v)
    return a, p


def two_moment_log_interval(mu, upper, variance_lower, terms=28):
    """Enclose p*log(a)+(1-p)*log(upper), an upper bound on E log X.

    The caller must prove X>0, EX=mu, X<=upper, and Var X>=variance_lower.
    Parameter validity alone does not establish those hypotheses.
    """
    a, p = two_moment_parameters(mu, upper, variance_lower)
    return linear_logs(((p, a), (1 - p, F(upper))), terms=terms)


def subset_epsilon(n, m):
    """The prior explicit strict-subset slack; zero for a singleton."""
    _integer(n, 3, 'n')
    _integer(m, 1, 'm')
    if m >= n:
        raise ValueError('proper subset required')
    if m == 1:
        return F(0)
    B = F(3 * n - m - 1, m - 1)
    degree = m * (n - m)
    geometric_sum = (B ** degree - 1) / (B - 1)
    return 1 / (2 * (1 + B) * geometric_sum)


def subset_parameters(n, m, explicit_slack=True):
    """Existing universal mean, positive variance lower, and subset support.

    Valid at global maxima, n>=4 and 2<=m<=n-2.  The optional old rational
    subset slack is used in both the point-radius and complement cuts.
    """
    _integer(n, 4, 'n')
    _integer(m, 2, 'm')
    if m > n - 2:
        raise ValueError('positive universal variance needs m<=n-2')
    h = n - m
    point_slack = subset_epsilon(n, n - 1) if explicit_slack else F(0)
    complement_slack = subset_epsilon(n, h) if explicit_slack else F(0)
    mu = F(n * (m - 1), n - 1)
    upper = min(F(n), F(2 * m, n) * (n - 1 - point_slack),
                F(n) - F(h * (h - 1) + 2 * complement_slack, n - 1))
    variance = F(n * (m - 1) * h * (h - 1),
                 m * (n - 1) ** 2 * (n - 2))
    a, p = two_moment_parameters(mu, upper, variance)
    return {'mu': mu, 'L': upper, 'v': variance, 'a': a, 'p': p}


def recurrence_intervals(n, m, log_Mm_terms, log_Mm_constant=0,
                         explicit_slack=True, terms=28):
    """New and prior upper-bound expressions, given a PROVED log M_m upper.

    log_Mm_terms is an iterable of rational (coefficient, argument) pairs;
    together with log_Mm_constant it must bound log M_m from above.
    Returns enclosures of the two formula values, not an enclosure of M_n.
    """
    data = subset_parameters(n, m, explicit_slack)
    mu, upper, v, a, p = (data[key] for key in ('mu', 'L', 'v', 'a', 'p'))
    N, Nm = F(n * (n - 1), 2), F(m * (m - 1), 2)
    cardinality_ratio = N / Nm
    common = [(cardinality_ratio * F(c), F(x)) for c, x in log_Mm_terms]
    constant = cardinality_ratio * F(log_Mm_constant)
    new_terms = common + [(N * p, a), (N * (1 - p), upper), (-N, F(m))]
    T = upper / mu
    discount = N * v / (mu * mu * (T - 1) ** 2)
    old_terms = common + [(N, mu / m), (discount, T)]
    old_constant = constant - discount * (T - 1)
    data['new'] = linear_logs(new_terms, constant, terms)
    data['prior'] = linear_logs(old_terms, old_constant, terms)
    return data


def disk_log_interval(n, terms=28):
    """Existing disk bound log U_n; no claim of a new disk constant."""
    _integer(n, 2, 'n')
    N = F(n * (n - 1), 2)
    return linear_logs(((F(n, 2), F(n)), (N, F(2 * n, n - 1))),
                       -N / 2, terms)


def known_small_log_terms(m):
    """Proved sharp M_m from the repository; m=6,7 retain external input."""
    if type(m) is not int or not 2 <= m <= 7:
        raise ValueError('known sharp inputs restricted to 2<=m<=7')
    if m <= 5:
        return [(F(m), F(m))]
    return [(F(m * (m - 1), 2), F(m)),
            (-F((m - 1) * (m - 2), 2), F(m - 1))]


def bounded_comparison(max_n=30):
    """Exact comparisons only for 4<=n<=max_n and known 2<=m<=7 inputs."""
    _integer(max_n, 4, 'max_n')
    records = []
    for n in range(4, max_n + 1):
        disk = disk_log_interval(n)
        for m in range(2, min(7, n - 2) + 1):
            data = recurrence_intervals(n, m, known_small_log_terms(m))
            new, old = data.pop('new'), data.pop('prior')
            improvement = (old[0] - new[1], old[1] - new[0])
            disk_gap = (new[0] - disk[1], new[1] - disk[0])
            records.append({'n': n, 'm': m,
                            'parameters': {key: str(val) for key, val in data.items()},
                            'prior_minus_new': [str(x) for x in outward_decimal(improvement)],
                            'new_minus_disk': [str(x) for x in outward_decimal(disk_gap)],
                            'strictly_improves_prior': improvement[0] > 0,
                            'new_still_above_disk': disk_gap[0] > 0})
    return records


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-n', type=int, default=30)
    options = parser.parse_args()
    print(json.dumps({'scope': 'Formula comparison, not subset enumeration.',
                      'records': bounded_comparison(options.max_n)}, indent=2))
