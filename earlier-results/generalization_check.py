#!/usr/bin/env python3
"""Exact supporting checks for the general variational estimates.

This program checks algebraic identities and the rational constants in the
six-point local-gap argument.  It is not a substitute for the analytic
arguments in the manuscript and is not a certificate of a global six-point
discriminant bound.  Only the Python standard library is required.
"""

from fractions import Fraction as F
from pathlib import Path
import json


class Poly:
    """Small exact multivariate polynomial, with named indeterminates."""
    def __init__(self, terms=None):
        self.terms = {m: F(c) for m, c in (terms or {}).items() if c}

    @staticmethod
    def var(name):
        return Poly({((name, 1),): F(1)})

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Poly) else Poly({(): F(value)})

    def __add__(self, other):
        other = self.coerce(other)
        result = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            result[monomial] = result.get(monomial, F(0)) + coefficient
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        result = {}
        for m, c in self.terms.items():
            for v, d in other.terms.items():
                powers = dict(m)
                for name, exponent in v:
                    powers[name] = powers.get(name, 0) + exponent
                monomial = tuple(sorted(powers.items()))
                result[monomial] = result.get(monomial, F(0)) + c*d
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Only nonnegative integer powers are supported")
        result = self.coerce(1)
        for _ in range(exponent):
            result = result*self
        return result


def main():
    rational_checks = {}
    identity_checks = []

    def less(name, lhs, rhs):
        lhs, rhs = F(lhs), F(rhs)
        if not lhs < rhs:
            raise ValueError("Failed strict inequality: " + name)
        rational_checks[name] = {
            "lhs": str(lhs), "rhs": str(rhs), "margin": str(rhs-lhs)
        }

    def identity(name, lhs, rhs):
        if (Poly.coerce(lhs)-Poly.coerce(rhs)).terms:
            raise ValueError("Failed polynomial identity: " + name)
        identity_checks.append(name)

    # Six points, 0 < s <= 1/200. Squaring these positive constants
    # verifies the only square-root comparisons needed in this argument.
    less("radius_product_deviation", F(1, 600), F(1, 24)**2)
    separation = (F(12, 25)-F(3, 1150))/F(25, 6)
    less("half_angle_sine_squared", F(1, 10), separation)
    less("incidence_factor_squared", F(16, 3), F(7, 3)**2)
    less("angular_error", F(7, 368), F(1, 50))
    less("sqrt_two_lower_squared", F(7, 5)**2, F(2))
    less("sqrt_three_over_two_lower_squared", F(6, 7)**2, F(3, 4))
    less("consecutive_pair_distance", F(9, 10), F(23, 6)*F(17, 35)**2)
    less("two_step_pair_distance", F(8, 3), F(23, 6)*F(59, 70)**2)
    less("opposite_pair_distance", F(100, 27), F(23, 6)*F(69, 70)**2)
    edge, two_step, opposite = F(10, 9), F(3, 8), F(27, 100)
    eigenvalues = [edge+3*two_step+2*opposite,
                  3*edge+3*two_step, 4*edge+2*opposite]
    less("six_point_laplacian", max(eigenvalues), F(5))

    # The polynomial identity underlying the general B z^k recurrence.
    # Summing it over j != i gives the recurrence without a CAS.
    x, y = Poly.var("x"), Poly.var("y")
    for k in range(1, 13):
        tail = sum((k-1-r)*x**(k-2-r)*y**r for r in range(k-1))
        identity("pairwise_B_recurrence_k="+str(k),
                 x**k-y**k,
                 k*x**(k-1)*(x-y)-(x-y)**2*tail)

    # This pair identity gives the mixed moment recurrence after pairing
    # the stationary equations multiplied by z_i^(k+1).
    for k in range(11):
        identity("mixed_moment_pair_identity_k="+str(k),
                 x**(k+1)-y**(k+1),
                 (x-y)*sum(x**(k-r)*y**r for r in range(k+1)))

    # Gram coefficients are checked as identities in n and independent
    # formal moment symbols.  The factors 4 n clear all denominators.
    n, s = Poly.var("n"), Poly.var("s")
    T, U, V = (Poly.var(name) for name in ("T", "U", "V"))
    q, w2, m, p2, cross = (
        Poly.var(name) for name in ("q", "w2", "m", "p2", "cross"))
    a, b = n-1, n-2
    identity("mixed_moment_coefficient_for_B_z4", 2*(2*n-3)-3*b, n)
    identity("Gram_D22_times_4n",
             a*a*(n*(n+s)-q)-4*n*a*a*s,
             a*a*(n*n-3*n*s-q))
    B4_norm = 4*a*a*V-12*a*b*U+9*b*b*T-n*q
    identity("Gram_D44_times_4n",
             n*a*a*V-a*a*w2-4*n*B4_norm,
             -15*n*a*a*V+48*n*a*b*U-36*n*b*b*T+4*n*n*q-a*a*w2)
    B2_B4_inner = 2*a*a*m-(2*a+3*b)*(2*n-3)*p2+3*a*b*p2
    identity("Gram_D24_times_4n",
             n*a*a*m-a*a*cross-4*n*B2_B4_inner,
             -7*n*a*a*m+4*n*(7*n*n-22*n+18)*p2-a*a*cross)

    # Algebraic consequences of projecting away the constant vector:
    # ||P e_i||^2 = (n-1)/n, lambda=(n-1)/2.  These coefficient identities
    # check the sharper entry bound and the squared column-norm bound.
    identity("projected_entry_bound_times_2n", a*(n-1), a*a)
    identity("projected_column_bound_times_4n", a*a*(n-1), a**3)

    # A further projected Cauchy inequality is useful for restricting q:
    # q <= n*s*(n+s)/(s+n*(n-2)^2/(n-1)^2), for s>0.
    # Verify its rearrangement and its five-point specialization.
    identity("projected_moment_q_rearrangement",
             a*a*n*s*(n+s)-a*a*s*q-n*b*b*q,
             a*a*n*s*(n+s)-(a*a*s+n*b*b)*q)
    identity("five_point_q_bound_numerator", 16*5*s*(5+s), 80*s*(5+s))
    identity("five_point_q_bound_denominator", 16*s+5*9, 16*s+45)

    # Exact centered regular-pentagon-plus-origin counterexample for n=6.
    # Its normalized ring radius squared is 6/5.  The ring discriminant
    # contributes 5^5 and total scaling contributes (6/5)^15.
    six_value = F(5**5)*F(6, 5)**15
    six_regular_value = F(6**6)
    less("six_point_counterexample", six_regular_value, six_value)
    six_ratio = six_value/six_regular_value
    if six_ratio != F(6**9, 5**10):
        raise ValueError("Incorrect six-point counterexample ratio")

    result = {
        "status": "passed",
        "scope": "Supporting algebra and constants; not a replacement for the analytic proof",
        "six_point_radial_gap_threshold": "1/200",
        "rational_checks": rational_checks,
        "formal_polynomial_identity_checks": identity_checks,
        "six_point_comparison_laplacian_eigenvalues": [str(v) for v in eigenvalues],
        "six_point_counterexample_ratio": str(six_ratio)
    }
    output = Path(__file__).with_name("generalization_result.json")
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
