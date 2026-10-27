"""Exact independent subset-variance checks; standard library only."""

from fractions import Fraction as F
from itertools import combinations
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def squared(v):
    return v[0] * v[0] + v[1] * v[1]


def centered(points):
    n = len(points)
    ax = sum(F(p[0]) for p in points) / n
    ay = sum(F(p[1]) for p in points) / n
    return [(F(p[0]) - ax, F(p[1]) - ay) for p in points]


def variance_formula(n, m, r4, p2):
    h = n - m
    base = F(n * (m - 1) * h * (h - 1), m * (n - 1)**2 * (n - 2))
    cr = F((m - 1) * h * ((n - 1) * (m - 1) - 2),
           m * (n - 1) * (n - 2) * (n - 3))
    cp = F((m - 1) * h * (h - 1),
           m * n * (n - 1) * (n - 2) * (n - 3))
    return base + cr * (r4 - n) + cp * p2, base


def check_configuration(points, totals):
    """Directly center and normalize squared lengths, without irrational roots."""
    z = centered(points)
    n = len(z)
    q = sum(squared(p) for p in z)
    require(q > 0, "constant configuration")
    scale = F(n, 1) / q
    radii = [scale * squared(p) for p in z]
    r4 = sum(r*r for r in radii)
    real_p2 = sum(p[0]**2 - p[1]**2 for p in z) * scale
    imag_p2 = sum(2*p[0]*p[1] for p in z) * scale
    p2 = real_p2**2 + imag_p2**2
    pair_fourth = F(0)
    for i, j in combinations(range(n), 2):
        d2 = scale * squared((z[i][0] - z[j][0], z[i][1] - z[j][1]))
        pair_fourth += d2*d2
    require(pair_fourth == n*r4 + 2*n*n + p2, "fourth-distance identity")
    totals["fourth_distance_identities"] += 1
    for m in range(2, n):
        subset_values = []
        for ids in combinations(range(n), m):
            bx = sum(z[i][0] for i in ids) / m
            by = sum(z[i][1] for i in ids) / m
            direct_q = scale * sum(squared((z[i][0] - bx, z[i][1] - by)) for i in ids)
            subset_values.append(direct_q)
        mu = sum(subset_values) / len(subset_values)
        var = sum((v-mu)**2 for v in subset_values) / len(subset_values)
        require(mu == F(n*(m-1), n-1), "subset mean identity")
        predicted, base = variance_formula(n, m, r4, p2)
        require(var == predicted, "subset variance identity")
        require(var >= base, "universal variance lower bound")
        if m == n-1:
            require(var == F(n, (n-1)**2)*(r4-n), "deletion endpoint")
            totals["deletion_endpoints"] += 1
        if r4 == n and p2 == 0:
            require(var == base, "isotropic equal-radius equality")
            totals["equality_cases"] += 1
        totals["subset_variance_identities"] += 1
        totals["subsets_enumerated"] += len(subset_values)


def run():
    totals = dict(fourth_distance_identities=0, subset_variance_identities=0,
                  deletion_endpoints=0, equality_cases=0, subsets_enumerated=0,
                  coefficient_sign_checks=0)
    for n in range(4, 10):
        cases = [
            [(i, i*i) for i in range(n)],
            [(i, 0) for i in range(n)],
            [((i*i+3*i+1) % 17 - 8, (i*i*i+5*i+2) % 19 - 9) for i in range(n)],
        ]
        for points in cases:
            check_configuration(points, totals)
    square = [(F(1), F(0)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(-1))]
    rotated = [(F(3, 5)*x-F(4, 5)*y, F(4, 5)*x+F(3, 5)*y) for x, y in square]
    check_configuration(square, totals)
    check_configuration(square+rotated, totals)
    for n in range(4, 101):
        for m in range(2, n):
            h = n-m
            require((n-1)*(m-1)-2 > 0, "positive fourth-moment coefficient")
            require(h*(h-1) >= 0, "nonnegative p2 coefficient")
            totals["coefficient_sign_checks"] += 1
    return {"passed": True, "arithmetic": "fractions.Fraction", **totals}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
