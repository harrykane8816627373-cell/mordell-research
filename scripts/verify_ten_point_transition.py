"""Exact n=10 feasible witness from Gueron--Shafrir, as recorded by Ji--Xu--Yu.

No floating point, dependencies, stationarity or optimality assumption is used.
Run from any directory, including with -B -S -O.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def verify():
    points = [(-11183, 0), (11183, 0)]
    for a, b in ((32765, 14669), (12323, 31572)):
        points.extend((sx*a, sy*b) for sx, sy in product((-1, 1), repeat=2))
    require(len(set(points)) == 10, 'ten distinct integer points')
    require(sum(x for x, y in points) == sum(y for x, y in points) == 0,
            'centered witness')
    mass = sum(x*x+y*y for x, y in points)
    require(mass == 9999608174, 'exact witness mass')
    distances = [(x-u)**2+(y-v)**2 for (x, y), (u, v) in combinations(points, 2)]
    require(len(distances) == 45 and min(distances) > 0, '45 positive distances')
    delta = 1
    for distance in distances:
        delta *= distance
    # z=sqrt(10/mass)*W has Q=10 and Delta(z)=(10/mass)^45*delta.
    # C_10=10^45/9^36. These positive factors cancel exactly.
    lhs, rhs = 9**36*delta, mass**45
    require(lhs > rhs, 'strict improvement over center plus regular nonagon')
    ratio = Fraction(lhs, rhs)
    scale = 10**8
    lower = ratio.numerator*scale//ratio.denominator
    require(Fraction(lower, scale) <= ratio < Fraction(lower+1, scale),
            'certified rational ratio interval')
    return {
        'passed': True,
        'arithmetic': 'Python integers and Fraction; checks survive -O',
        'n': 10,
        'points': points,
        'mass': mass,
        'pair_count': len(distances),
        'strict_integer_comparison': '9^36 * product(squared distances) > mass^45',
        'delta_over_C10_interval': [str(Fraction(lower, scale)), str(Fraction(lower+1, scale))],
        'integer_margin_sha256': sha256(str(lhs-rhs).encode('ascii')).hexdigest(),
        'conclusion': 'M_10 > C_10; no claim of global optimality of the witness',
        'source': 'Gueron--Shafrir (1999), coordinates recorded in arXiv:2611.03088v1, Proposition 5.1',
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
