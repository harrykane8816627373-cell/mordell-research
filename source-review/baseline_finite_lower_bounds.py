"""Exact finite-n lower bounds for the normalized squared Vandermonde.
All acceptance arithmetic is integer/Fraction; no floating maximization.
Coefficient moments are existing project results. Partition DP is an
implementation, not a replacement for the all-n proof in README.md.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import permutations
from math import factorial, comb


def require_n(n):
    if not isinstance(n, int) or n < 2:
        raise ValueError('n must be an integer at least 2')


def perm_sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))) % 2 else 1


class CoefficientMoments:
    """Finite exact coefficient norm, with permutation orbits combined.
    Cost grows rapidly with n and order. No claim of polynomial complexity.
    V=det(z_i**j) for j=0,...,n-1; its orientation has no effect on norms.
    """
    def __init__(self, n):
        require_n(n)
        self.n = n
        self.N = n * (n - 1) // 2
        self._permutations = None

    @property
    def signed_permutations(self):
        if self._permutations is None:
            self._permutations = tuple((p, perm_sign(p)) for p in permutations(range(self.n)))
        return self._permutations

    @lru_cache(maxsize=None)
    def _square_coefficient(self, alpha):
        """Choose two exponent permutations with prescribed component sums."""
        n = self.n
        @lru_cache(maxsize=None)
        def solve(i, used_a, used_b):
            if i == n:
                return 1
            total = 0
            for a in range(n):
                b = alpha[i] - a
                if used_a >> a & 1 or not 0 <= b < n or used_b >> b & 1:
                    continue
                inversions = (used_a >> (a + 1)).bit_count() + (used_b >> (b + 1)).bit_count()
                total += (-1 if inversions % 2 else 1) * solve(i + 1, used_a | (1 << a), used_b | (1 << b))
            return total
        return solve(0, 0, 0)

    def coefficient(self, order, alpha):
        n = self.n
        if len(alpha) != n or sum(alpha) != order * self.N or min(alpha) < 0 or max(alpha) > order * (n - 1):
            return 0
        if order % 2 and len(set(alpha)) < n:
            return 0
        parity = sum(alpha[i] < alpha[j] for i in range(n) for j in range(i + 1, n)) % 2
        canonical = tuple(sorted(alpha, reverse=True))
        sign = -1 if order % 2 and parity else 1
        return sign * self._canonical_coefficient(order, canonical)

    @lru_cache(maxsize=None)
    def _canonical_coefficient(self, order, alpha):
        n = self.n
        if order == 0:
            return int(not any(alpha))
        if order == 1:
            return (-1) ** self.N if alpha == tuple(range(n - 1, -1, -1)) else 0
        if order == 2:
            return self._square_coefficient(alpha)
        total = 0
        for p, sign in self.signed_permutations:
            remainder = tuple(a - b for a, b in zip(alpha, p))
            total += sign * self.coefficient(order - 1, remainder)
        return total

    def canonical_exponents(self, order):
        """Dominated partitions; odd powers have distinct exponents."""
        n, total = self.n, order * self.N
        strict = order % 2
        def visit(prefix, remaining, maximum):
            i = len(prefix)
            slots = n - i
            if not slots:
                if remaining == 0:
                    yield tuple(prefix)
                return
            maximum = min(maximum, remaining)
            minimum_tail = (slots - 1) * (slots - 2) // 2 if strict else 0
            for value in range(maximum, -1, -1):
                tail = remaining - value
                next_max = value - strict
                if tail < minimum_tail or tail > (slots - 1) * max(next_max, 0):
                    continue
                prefix_sum = total - remaining + value
                dominance = order * (i + 1) * (2 * n - i - 2) // 2
                if prefix_sum > dominance:
                    continue
                yield from visit(prefix + [value], tail, next_max)
        yield from visit([], total, order * (n - 1))

    @lru_cache(maxsize=None)
    def coefficient_norm(self, order):
        if order < 0:
            raise ValueError('nonnegative moment order required')
        if order == 0:
            return 1
        if order == 1:
            return factorial(self.n) * __import__('math').prod(factorial(j) for j in range(self.n))
        total = 0
        for alpha in self.canonical_exponents(order):
            c = self._canonical_coefficient(order, alpha)
            if c == 0:
                continue
            counts = {a: alpha.count(a) for a in set(alpha)}
            orbit = factorial(self.n)
            for multiplicity in counts.values():
                orbit //= factorial(multiplicity)
            total += orbit * c * c * __import__('math').prod(factorial(a) for a in alpha)
        return total

    def moment(self, order):
        return Fraction(self.n ** (order * self.N) * factorial(self.n - 2) * self.coefficient_norm(order), factorial(order * self.N + self.n - 2))

    def ratio(self, order):
        return self.moment(order + 1) / self.moment(order)

    def rayleigh_lower(self):
        """Explicit rational improvement over ratio mu3/mu2 for n>=3."""
        if self.n == 2:
            return Fraction(4), Fraction(0)
        m1, m2, m3 = (self.moment(j) for j in (1, 2, 3))
        A, B = m1 * m3 - m2 * m2, m3 - m1 * m2
        if A <= 0 or B <= 0:
            raise ArithmeticError('nonconstant moment positivity failed')
        a = A / B
        denominator = m2 - 2 * a * m1 + a * a
        result = (m3 - 2 * a * m2 + a * a * m1) / denominator
        if denominator <= 0 or result <= m3 / m2:
            raise ArithmeticError('Rayleigh improvement check failed')
        return result, a


def gaussian_lower(n):
    require_n(n)
    N = n * (n - 1) // 2
    product = __import__('math').prod(factorial(j) for j in range(n))
    return Fraction(factorial(n) * factorial(n - 2) * n ** N * product, factorial(N + n - 2))


def central_polygon_lower(n):
    require_n(n)
    return Fraction(n ** (n * (n - 1) // 2), (n - 1) ** ((n - 1) * (n - 2) // 2))


def hermite_lower(n):
    require_n(n)
    return Fraction(__import__('math').prod(j ** j for j in range(1, n + 1)), (n - 1) ** (n * (n - 1) // 2))


def subset_slack(n, m):
    require_n(n)
    if not 2 <= m < n:
        raise ValueError('strict subset result requires 2 <= m < n')
    D = m * (n - m)
    B = Fraction(3 * n - m - 1, m - 1)
    epsilon = 1 / (2 * (1 + B) * sum((B ** j for j in range(D)), Fraction(0)))
    return B, epsilon


def stopping_check(n, k, epsilon):
    """Exact sufficient test for B_k <= M_n <= (1+eps)*B_k."""
    require_n(n)
    epsilon = Fraction(epsilon)
    if k < 0 or epsilon <= 0:
        raise ValueError('k >= 0 and epsilon > 0 required')
    ell = k + 1
    D = comb(ell * n * (n - 1) // 2 + n - 2, n - 2)
    a, b = epsilon.numerator, epsilon.denominator
    return D * b ** ell <= (a + b) ** ell


def finite_bound_record(n, moment_order=0):
    """A rational lower bound for every integer n>=2; optional exact moments."""
    from triangular_lattice import lattice_certificate, triangular_norm, gaussian_norm, improved_centered_polygon
    candidates = {
        'regular_polygon': Fraction(n ** n),
        'centered_polygon': central_polygon_lower(n),
        'gaussian_sphere_average': gaussian_lower(n),
        'square_integer_disk': lattice_certificate(n, gaussian_norm)['value'],
        'triangular_integer_disk': lattice_certificate(n, triangular_norm)['value'],
    }
    old_improvement = improved_centered_polygon(n)
    if old_improvement is not None:
        candidates['existing_explicit_center_improvement'] = old_improvement
    moments = []
    if moment_order:
        engine = CoefficientMoments(n)
        moments = [engine.moment(k) for k in range(moment_order + 1)]
        for k in range(moment_order):
            candidates['moment_ratio_' + str(k)] = moments[k + 1] / moments[k]
        if moment_order >= 3:
            candidates['rational_rayleigh'] = engine.rayleigh_lower()[0]
    best = max(candidates.values())
    return {
        'n': n,
        'moment_order_computed': moment_order,
        'candidates': {key: str(value) for key, value in candidates.items()},
        'moments': [str(value) for value in moments],
        'lower_bound': str(best),
        'winners': [key for key, value in candidates.items() if value == best],
        'global_sharpness_claimed': False,
    }


def main():
    import argparse,json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int, nargs='+')
    parser.add_argument('--moment-order', type=int, choices=(0, 1, 2, 3), default=0,
                        help='0 uses fast constructive/average bounds; 3 adds Rayleigh. Exact high moments grow rapidly in cost.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rows = [finite_bound_record(n, args.moment_order) for n in args.n]
    value = json.dumps({'format': 'mordell-finite-lower-v1', 'rows': rows}, indent=2)
    if args.output:
        args.output.write_text(value + '\n', encoding='utf-8')
    else:
        print(value)


if __name__ == '__main__':
    main()
