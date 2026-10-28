"""Portable all-n rational polar value and audited sign/centering refinement.

The original polar family is from Fengcan Jiang's supplied revised source.
No dependency beyond the Python standard library is required.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd, prod


def polar_rings(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n >= 2 required')
    k = 0
    while 3*(k+1)*(k+2) <= n:
        k += 1
    n0 = 3*k*(k+1)
    r = n-n0
    rings = [(6*j, 6*j*j) for j in range(1, k+1)]
    if r:
        rings.append((r, n+n0))
    if sum(m for m, b in rings) != n or sum(m*b for m, b in rings) != n*n:
        raise ArithmeticError('population or prescribed moment failed')
    return n0, r, tuple(rings)


def ring_value(n, signs=None, center_singleton=False):
    n0, r, rings = polar_rings(n)
    signs = (1,)*len(rings) if signs is None else tuple(signs)
    if len(signs) != len(rings) or any(s not in (-1, 1) for s in signs):
        raise ValueError('one +/-1 sign required per ring')
    N = n*(n-1)//2
    numerator = prod(m**m*b**(m*(m-1)//2) for m, b in rings)
    for i, j in combinations(range(len(rings)), 2):
        m, a = rings[i]
        mm, b = rings[j]
        g = gcd(m, mm)
        L = m*mm//g
        if L % 2:
            raise ArithmeticError('unexpected odd cross exponent')
        cross = abs(signs[j]**(m//g)*b**(L//2)
                    - signs[i]**(mm//g)*a**(L//2))
        if cross < b**(L//2)-a**(L//2):
            raise ArithmeticError('phase changed a difference to something smaller')
        numerator *= cross**(2*g)
    value = Fraction(numerator, n**N)
    if center_singleton and r == 1:
        value *= Fraction(n**3, n**3-n-n0)**N
    return value


def refined_polar(n):
    n0, r, rings = polar_rings(n)
    minus = (-1,)*len(rings)
    flip = (1,)+minus[1:]
    candidates = [(ring_value(n, signs, True), signs) for signs in (minus, flip)]
    value, signs = max(candidates, key=lambda p:p[0])
    return {'n':n, 'rings':rings, 'signs':signs, 'singleton_centered':r == 1,
            'value':value, 'original_value':ring_value(n)}


def twelve_two_ring_maximum():
    """Maximum over two concentric regular hexagons, all radii/phases.

    Only this family is solved: u+v=2, uv=20/27 and opposite sixth powers.
    The unrestricted twelve-point global maximum is not identified.
    """
    return Fraction(2**102*5**15, 3**57)


def verify(max_n=100):
    count = 0
    for n in range(2, max_n+1):
        record = refined_polar(n)
        original, improved = record['original_value'], record['value']
        if improved < original or (n >= 7 and improved <= original):
            raise ArithmeticError('uniform rational improvement failed')
        if n <= 5 and original != n**n:
            raise ArithmeticError('small polygon value failed')
        count += 1
    if ring_value(18) != Fraction(2**162*455**12, 3**111):
        raise ArithmeticError('original eighteen-point source value failed')
    if refined_polar(12)['value'] != Fraction(3**27*7**12, 2**30):
        raise ArithmeticError('twelve-point refinement failed')
    t = Fraction(20, 27)
    direct = 6**12*t**15*(8-6*t)**12
    if direct != twelve_two_ring_maximum() or direct <= refined_polar(12)['value']:
        raise ArithmeticError('optimized twelve-point family value failed')
    if 120-162*t != 0 or not 0 < t < 1:
        raise ArithmeticError('twelve-point critical parameter failed')
    return {'passed':True,'point_counts_checked':count,'n_min':2,'n_max':max_n,
            'arithmetic':'integers and Fraction only',
            'scope':'Finite checks exercise the formula; POLAR_REFINEMENT.md gives the all-n proof.'}


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int, nargs='?')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        result = verify()
    elif args.n is not None:
        result = refined_polar(args.n)
        result['value'] = str(result['value'])
        result['original_value'] = str(result['original_value'])
        if args.n == 12:
            result['two_equal_six_ring_family_maximum'] = str(twelve_two_ring_maximum())
    else:
        parser.error('specify n or --verify')
    print(json.dumps(result, ensure_ascii=False, indent=2))
