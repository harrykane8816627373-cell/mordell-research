"""Exhaustive integer-grid bounds for the planar squared Vandermonde.

No floating-point optimizer is used. A grid value is a lower bound;
the analytic rounding theorem supplies its global upper certificate.
The complete search is combinatorial and rapidly becomes impractical.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, factorial, isqrt, prod


def check_n(n):
    if type(n) is not int or n < 2:
        raise ValueError('n must be an integer >= 2')


def ceil_sqrt_ratio(a, b):
    if a < 0 or b <= 0:
        raise ValueError('nonnegative numerator, positive denominator required')
    r = isqrt(a // b)
    return r + int(r * r * b < a)


def grid_radius(n, q):
    check_n(n)
    if type(q) is not int or q < 1:
        raise ValueError('q must be a positive integer')
    return ceil_sqrt_ratio(2 * (n - 1) * q * q, n)


def rational_value(points):
    """Value after centering and scaling integer points to Q=n."""
    n = len(points)
    check_n(n)
    if len(set(points)) != n:
        return Fraction(0)
    distances = [(a - c) ** 2 + (b - d) ** 2
                 for (a, b), (c, d) in combinations(points, 2)]
    D = sum(distances)
    return Fraction(n * n, D) ** (n * (n - 1) // 2) * prod(distances)


def subset_floor(n, m):
    """Published global-optimizer subset floor, including its rational slack."""
    if m == n:
        return Fraction(n)
    if m == 1:
        return Fraction(0)
    if not 2 <= m < n:
        raise ValueError('invalid subset size')
    B = Fraction(3 * n - m - 1, m - 1)
    epsilon = 1 / (2 * (1 + B) * sum((B ** j for j in range(m * (n - m))), Fraction(0)))
    return (m * (m - 1) + 2 * epsilon) / (n - 1)


def passes_rounded_subset_cuts(points, q):
    """Safe cuts for an unnormalized nearest-grid image of a global maximum.

    These weakened inequalities incorporate rounding error. They do not
    assert that arbitrary grid configurations obey optimizer conditions.
    """
    n = len(points)
    check_n(n)
    grid_radius(n, q)
    for m in range(2, n + 1):
        threshold = m * q * (q - 1) * subset_floor(n, m) - Fraction(m * m * (q - 1), 2)
        for block in combinations(points, m):
            D = sum((a - c) ** 2 + (b - d) ** 2
                    for (a, b), (c, d) in combinations(block, 2))
            if D < threshold:
                return False
    return True


def exhaustive_grid(n, q, max_configurations=1000000):
    """Enumerate every n-element subset of the specified integer square.

    A budget refusal raises before any search; it never returns a partial
    sample as if it were an exhaustive certificate. Labels need not be
    enumerated since the objective is permutation invariant.
    """
    A = grid_radius(n, q)
    points = [(a, b) for a in range(-A, A + 1) for b in range(-A, A + 1)]
    if len(points) < n:
        raise ValueError('integer box has fewer than n distinct points')
    count = comb(len(points), n)
    if max_configurations is not None and count > max_configurations:
        raise ValueError('complete grid has %s configurations, above budget %s'
                         % (count, max_configurations))
    N = n * (n - 1) // 2
    # Independent integer comparisons avoid repeated Fraction reductions.
    distance = {(i, j): (points[i][0] - points[j][0]) ** 2
                         + (points[i][1] - points[j][1]) ** 2
                for i in range(len(points)) for j in range(i + 1, len(points))}
    best_num, best_den, best_indices, visited = 0, 1, None, 0
    for indices in combinations(range(len(points)), n):
        ds = [distance[i, j] for i, j in combinations(indices, 2)]
        numerator, denominator = prod(ds), sum(ds) ** N
        if numerator * best_den > best_num * denominator:
            best_num, best_den, best_indices = numerator, denominator, indices
        visited += 1
    if visited != count or best_indices is None:
        raise ArithmeticError('incomplete exhaustive grid')
    witness = tuple(points[i] for i in best_indices)
    value = Fraction(n ** (2 * N) * best_num, best_den)
    if value != rational_value(witness):
        raise ArithmeticError('witness normalization mismatch')
    return {'n': n, 'q': q, 'A': A, 'grid_points': len(points),
            'configurations': visited, 'witness': [list(p) for p in witness],
            'value': str(value), 'complete': True}


def rounding_exponent(n, q):
    """Return a rational eta with G <= M <= G*exp(eta), n>=3.

    A rational upper bound tau on sqrt((n-1)/2)/q is used throughout.
    """
    check_n(n)
    grid_radius(n, q)  # validate q even when this function is called alone
    if n == 2:
        return Fraction(0), Fraction(0)
    if 2 * q * q <= n - 1:
        raise ValueError('rounding theorem requires 2*q*q > n-1')
    scale = 65536
    while True:
        tau = Fraction(ceil_sqrt_ratio((n - 1) * scale * scale, 2 * q * q), scale)
        if tau < 1:
            break
        scale *= 2
    simple_tau = Fraction(ceil_sqrt_ratio(n - 1, 2), q)
    if simple_tau < 1:
        tau = min(tau, simple_tau)
    eta = Fraction(n * (n - 1), 2 * q * q) * (1 + Fraction(n, 6) * tau / (1 - tau))
    return eta, tau


def exp_upper(x, order=24):
    """Exact rational Taylor upper bound for exp(x), x>=0."""
    x = Fraction(x)
    if x < 0:
        raise ValueError('x must be nonnegative')
    if type(order) is not int or order < 0:
        raise ValueError('nonnegative integer order required')
    while x >= order + 2:
        order += 1
    term = Fraction(1)
    total = term
    for j in range(1, order + 1):
        term *= x / j
        total += term
    omitted = term * x / (order + 1)
    return total + omitted / (1 - x / (order + 2))


def relative_grid_plan(n, epsilon):
    """Sufficient q for M <= (1+epsilon)*G, before running the grid.

    Planning does not execute or certify the combinatorial search itself.
    """
    check_n(n)
    epsilon = Fraction(epsilon)
    if epsilon <= 0:
        raise ValueError('epsilon must be positive')
    if n == 2:
        q = 1
    else:
        a = ceil_sqrt_ratio(n - 1, 2)
        target = Fraction(7 * n * (n - 1), 12) * (1 + epsilon) / epsilon
        q = max((n + 1) * a, ceil_sqrt_ratio(target.numerator, target.denominator))
        eta, _ = rounding_exponent(n, q)
        if eta > epsilon / (1 + epsilon):
            raise ArithmeticError('integer stopping plan failed')
    A = grid_radius(n, q)
    return {'n': n, 'epsilon': str(epsilon), 'q': q, 'A': A,
            'grid_points': (2 * A + 1) ** 2,
            'configurations': str(comb((2 * A + 1) ** 2, n)),
            'search_executed': False}


def certified_record(n, q, max_configurations=1000000):
    eta, tau = rounding_exponent(n, q)
    row = exhaustive_grid(n, q, max_configurations)
    value = Fraction(row['value'])
    factor = exp_upper(eta)
    row.update(eta=str(eta), tau=str(tau),
               upper_factor=str(factor), upper_bound=str(value * factor))
    return row


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int)
    parser.add_argument('--q', type=int, default=2)
    parser.add_argument('--max-configurations', type=int, default=1000000)
    parser.add_argument('--plan', help='positive rational relative error; do not run search')
    args = parser.parse_args()
    result = relative_grid_plan(args.n, args.plan) if args.plan else certified_record(args.n, args.q, args.max_configurations)
    print(json.dumps(result, ensure_ascii=False, indent=2))
