"""Exact replay for the sharp two-moment subset-recursion supplement.

Checks are explicit exceptions, so python -O does not disable acceptance.
Finite cases check the implementation; SUBSET_MOMENT.md proves the theorem.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from subset_moment import (bounded_comparison, linear_logs, log_interval,
                           outward_decimal, subset_epsilon, subset_parameters,
                           two_moment_log_interval, two_moment_parameters)


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def independent_log(q, terms=56):
    """Separate direct atanh sum after repeated binary reductions."""
    q = F(q)
    if q < 1:
        lo, hi = independent_log(1 / q, terms)
        return -hi, -lo
    count = 0
    while q >= 2:
        q /= 2
        count += 1
    def direct(value):
        y = (value - 1) / (value + 1)
        lo = 2 * sum((y ** (2*j + 1) / (2*j + 1)
                      for j in range(terms)), F(0))
        tail = 2*y**(2*terms + 1) / ((2*terms + 1)*(1-y*y))
        return lo, lo + tail
    lo, hi = direct(q)
    lo2, hi2 = direct(F(2))
    return lo + count*lo2, hi + count*hi2


def independent_average(atoms):
    lo = hi = F(0)
    for x, probability in atoms:
        left, right = independent_log(x)
        lo += probability * left
        hi += probability * right
    return lo, hi


def moment_test(atoms, variance_factor=F(1)):
    require(sum(p for x, p in atoms) == 1, 'test probability mass')
    mu = sum(p*x for x, p in atoms)
    variance = sum(p*(x-mu)**2 for x, p in atoms)
    upper = max(x for x, p in atoms)
    v = variance * variance_factor
    a, p = two_moment_parameters(mu, upper, v)
    require(p*a+(1-p)*upper == mu, 'extremal law mean')
    require(p*(a-mu)**2+(1-p)*(upper-mu)**2 == v,
            'extremal law variance')
    actual = independent_average(atoms)
    envelope = two_moment_log_interval(mu, upper, v)
    at_equality = v == variance and all(x in (a, upper) for x, prob in atoms)
    if at_equality:
        require(actual[0] <= envelope[1] and envelope[0] <= actual[1],
                'sharp two-point equality enclosure')
    else:
        require(actual[1] < envelope[0], 'strict moment-log bound on test law')
    return at_equality


def centered_subset_law(points, m):
    n = len(points)
    xbar = sum(F(x) for x, y in points) / n
    ybar = sum(F(y) for x, y in points) / n
    centered = [(F(x)-xbar, F(y)-ybar) for x, y in points]
    total = sum(x*x+y*y for x, y in centered)
    counts = {}
    subsets = list(combinations(centered, m))
    for subset in subsets:
        sum_x = sum(x for x, y in subset)
        sum_y = sum(y for x, y in subset)
        q = (sum(x*x+y*y for x, y in subset)
             - (sum_x*sum_x+sum_y*sum_y)/m) * n / total
        counts[q] = counts.get(q, 0) + 1
    atoms = [(q, F(count, len(subsets))) for q, count in sorted(counts.items())]
    require(sum(q*p for q, p in atoms) == F(n*(m-1), n-1),
            'direct subset mean, with original normalization')
    return atoms


def main(write_report=False):
    log_checks = 0
    for q in (F(1), F(1, 1024), F(2), F(7, 5), F(5, 7),
              F(3, 2)**15, F(2)**100, F(123456789, 987654321)):
        left, right = log_interval(q)
        a, b = independent_log(q)
        require(left <= a <= b <= right, 'independent log enclosure')
        low, high = outward_decimal((left, right), 14)
        require(low <= left <= right <= high, 'outward interval rounding')
        log_checks += 1
    for q in (F(3, 7), F(2), F(19, 11)):
        lo, hi = linear_logs(((1, q), (1, 1/q)))
        require(lo <= 0 <= hi, 'log reciprocal cancellation')
        log_checks += 1

    laws = [
        [(F(1), F(1, 2)), (F(3), F(1, 2))],
        [(F(1, 5), F(1, 7)), (F(2), F(2, 7)), (F(4), F(4, 7))],
        [(F(1), F(1, 4)), (F(2), F(1, 2)), (F(5), F(1, 4))],
    ]
    families = [
        [(1, 0), (0, 1), (-1, 0), (0, -1)],
        [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)],
        [(x, y) for x in (-1, 0, 1) for y in (-1, 1)],
        [(x, y) for x in (-1, 0, 1) for y in (-1, 0, 1)
         if (x, y) != (0, 0)],
    ]
    enumerated_subsets = 0
    for points in families:
        for m in range(2, len(points)-1):
            law = centered_subset_law(points, m)
            enumerated_subsets += len(list(combinations(points, m)))
            if len(law) > 1:
                laws.append(law)
    equality_checks = 0
    for law in laws:
        equality_checks += int(moment_test(law))
        require(not moment_test(law, F(1, 2)), 'strict lower-variance test')

    # Check the sign of the Hermite quadratic coefficient, using interval logs.
    for a, upper in ((F(1), F(3)), (F(1, 7), F(5, 2)), (F(5), F(6))):
        loglo, loghi = log_interval(upper/a)
        require(loghi < (upper-a)/a, 'negative quadratic coefficient')
        for k in range(1, 13):
            x = F(k, 13)*upper
            if x == a:
                continue
            delta = x-a
            coefficient = delta*delta/(upper-a)**2
            # q(x)-log(x), keeping all logarithm signs outward.
            lo, hi = linear_logs(((1, a), (coefficient, upper/a), (-1, x)),
                                 delta/a-coefficient*(upper-a)/a)
            require(lo > 0, 'Hermite majorant at rational test point')

    bad_parameters = [(0, 1, F(1, 10)), (1, 1, F(1, 10)),
                      (1, 2, 0), (1, 2, 1), (1, 2, 2)]
    for args in bad_parameters:
        try:
            two_moment_parameters(*args)
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid moment parameters accepted')
    require(subset_epsilon(8, 7) == F(2187, 9217846),
            'prior subset slack copied with exact normalization')
    parameter_checks = 0
    for n in range(4, 101):
        for m in (2, n-2):
            strong = subset_parameters(n, m)
            weak = subset_parameters(n, m, explicit_slack=False)
            require(strong['L'] <= weak['L'], 'subset slack strengthens support')
            require(0 < strong['a'] < strong['mu'] < strong['L'],
                    'positive extremal support and probabilities')
            parameter_checks += 1

    comparison = bounded_comparison(30)
    require(len(comparison) == 147, 'bounded cardinality comparison count')
    require(all(row['strictly_improves_prior'] for row in comparison),
            '147 exact strict improvements of prior recurrence formula')
    require(all(row['new_still_above_disk'] for row in comparison),
            'do not claim disk-envelope improvement from this recurrence')
    report = {
        'status': 'PASS',
        'scope': 'Finite exact implementation checks; the general proof is SUBSET_MOMENT.md.',
        'log_interval_checks': log_checks,
        'moment_laws_checked': len(laws),
        'variance_choices_per_law': 2,
        'two_point_equalities_checked': equality_checks,
        'directly_enumerated_test_subsets': enumerated_subsets,
        'support_parameter_checks': parameter_checks,
        'formula_comparisons': len(comparison),
        'new_formula_beats_prior_in_all_comparisons': True,
        'new_formula_still_above_disk_in_all_comparisons': True,
        'known_M6_M7_input': 'Repository sharp values retain the external physical-stationary classification assumption.',
        'records': comparison,
    }
    if write_report:
        Path(__file__).with_name('subset_moment_comparison.json').write_text(
            json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: value for key, value in report.items() if key != 'records'},
                     ensure_ascii=False, indent=2))
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    main(args.write_report)
