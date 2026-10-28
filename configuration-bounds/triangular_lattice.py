#!/usr/bin/env python3
"""Exact finite-n lower bounds from triangular and Gaussian integer lattices.

No third-party dependencies. All deciding comparisons use integers/Fraction;
Decimal values in the human-readable table are display values only.

For triangular coordinates omega=exp(i*pi/3), q(a,b)=a*a+a*b+b*b.
The point rule is global lexicographic sorting by (q,a,b), including ties.
The finite box [-n,n]^2 contains the first n lattice points: outside the
box q>=(n+1)^2/2, while n axis points inside have q<=ceil((n-1)/2)^2.
The inequality between these two thresholds is strict. The same box works
for the square lattice q(a,b)=a*a+b*b.

If D=n*sum(q(v))-q(sum(v)), then D=sum_{i<j}q(v_i-v_j)>0.
Centering and multiplying the complex points by n/sqrt(D) gives moment n.
Thus (n*n/D)^N * product_{i<j}q(v_i-v_j) is a rational lower bound M_n.
This is a constructive lower bound, not a classification or sharpness claim.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from math import factorial
from pathlib import Path
from typing import Callable

Point = tuple[int, int]
Norm = Callable[[Point], int]


def triangular_norm(v: Point) -> int:
    a, b = v
    return a*a + a*b + b*b


def gaussian_norm(v: Point) -> int:
    a, b = v
    return a*a + b*b


def first_points(n: int, norm: Norm = triangular_norm) -> list[Point]:
    if n < 2:
        raise ValueError('n must be at least 2')
    points = [(a, b) for a in range(-n, n+1) for b in range(-n, n+1)]
    points.sort(key=lambda v: (norm(v), v[0], v[1]))
    chosen = points[:n]
    axis_threshold = ((n-1+1)//2)**2  # ceil((n-1)/2)^2
    if not (2*axis_threshold < (n+1)**2):
        raise ArithmeticError('finite enumeration proof threshold failed')
    if norm(chosen[-1]) > axis_threshold:
        raise ArithmeticError('first-n threshold exceeds axis witnesses')
    if len(set(chosen)) != n:
        raise ArithmeticError('repeated point in lattice rule')
    return chosen


def lattice_certificate(n: int, norm: Norm = triangular_norm) -> dict:
    points = first_points(n, norm)
    sum_a = sum(v[0] for v in points)
    sum_b = sum(v[1] for v in points)
    norm_sum = sum(norm(v) for v in points)
    inertia_integer = n*norm_sum - norm((sum_a, sum_b))
    pair_sum = 0
    pair_product = 1
    for i, u in enumerate(points):
        for v in points[i+1:]:
            distance_squared = norm((u[0]-v[0], u[1]-v[1]))
            if distance_squared <= 0:
                raise ArithmeticError('nonpositive pair distance')
            pair_sum += distance_squared
            pair_product *= distance_squared
    if inertia_integer != pair_sum or inertia_integer <= 0:
        raise ArithmeticError('inertia/pair identity failed')
    # A second centering check using integer numerators n*v-sum(v).
    centered_numerator_norm_sum = sum(
        norm((n*a-sum_a, n*b-sum_b)) for a, b in points
    )
    if centered_numerator_norm_sum != n*inertia_integer:
        raise ArithmeticError('centered moment identity failed')
    pair_count = n*(n-1)//2
    value = Fraction(n*n, inertia_integer)**pair_count * pair_product
    return {
        'points': points,
        'sum_coordinates': (sum_a, sum_b),
        'sum_norms': norm_sum,
        'D': inertia_integer,
        'pair_squared_distance_product': pair_product,
        'centered_numerator_norm_sum': centered_numerator_norm_sum,
        'value': value,
    }


def gaussian_average(n: int) -> Fraction:
    pair_count = n*(n-1)//2
    factorial_product = 1
    for k in range(n):
        factorial_product *= factorial(k)
    return Fraction(
        factorial(n)*factorial(n-2)*n**pair_count*factorial_product,
        factorial(pair_count+n-2),
    )


def polygon(n: int) -> Fraction:
    return Fraction(n**n)


def centered_polygon(n: int) -> Fraction:
    pair_count = n*(n-1)//2
    return Fraction(n**pair_count, (n-1)**((n-1)*(n-2)//2))


def hermite_collinear(n: int) -> Fraction:
    product = 1
    for k in range(1, n+1):
        product *= k**k
    return Fraction(product, (n-1)**(n*(n-1)//2))


def improved_centered_polygon(n: int) -> Fraction | None:
    if n == 10:
        u, v = Fraction(1, 8), Fraction(39, 32)
        return 4*u*8**8*v**28*(u**4+v**4)**4
    if n >= 11:
        m = n-1
        pair_count = n*(n-1)//2
        q = (m-2)**2//4
        b = 1 if m % 2 == 0 else 2
        A = Fraction(b*m*(q-m-3), 2)
        H = m+2
        t = min(Fraction(1, 2*H), A/(4*pair_count*H**3))
        x = Fraction(2, 3)*A*t*t
        if A <= 0 or x <= 0:
            raise ArithmeticError('old perturbation parameter not positive')
        return centered_polygon(n)*(1+x)
    return None


def encoded_fraction(value: Fraction) -> dict[str, str]:
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator)}


def decimal_string(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits+12
        x = Decimal(value.numerator)/Decimal(value.denominator)
        return format(x, f'.{digits-1}g')


def json_certificate(certificate: dict) -> dict:
    return {
        key: encoded_fraction(value) if isinstance(value, Fraction) else value
        for key, value in certificate.items()
    }


def audit(max_n: int = 30) -> tuple[dict, str]:
    rows = []
    table = [
        '# Exact triangular-lattice finite lower-bound audit',
        '',
        'All comparisons are exact rational comparisons. Ratios below are decimal displays only.',
        'The all-n bound follows from the centering and homogeneity proof in the script; this finite table does not prove a universal comparison between lattice families.',
        '',
        '| n | triangular / square-grid | triangular / prior envelope | envelope winner(s) |',
        '|---:|---:|---:|:---|',
    ]
    strict_improvements = []
    for n in range(2, max_n+1):
        triangular = lattice_certificate(n, triangular_norm)
        square = lattice_certificate(n, gaussian_norm)
        baseline = {
            'gaussian_average': gaussian_average(n),
            'regular_polygon': polygon(n),
            'centered_polygon': centered_polygon(n),
            'gaussian_integer_disk': square['value'],
        }
        improved = improved_centered_polygon(n)
        if improved is not None:
            baseline['old_explicit_center_improvement'] = improved
        collinear = hermite_collinear(n)
        if collinear > polygon(n) or (n >= 3 and collinear >= polygon(n)):
            raise ArithmeticError('Hermite redundancy check failed')
        prior_max = max(baseline.values())
        combined = dict(baseline)
        combined['triangular_integer_disk'] = triangular['value']
        combined_max = max(combined.values())
        winners = [name for name, value in combined.items() if value == combined_max]
        stronger = triangular['value'] > prior_max
        if stronger:
            strict_improvements.append(n)
        if n == 2 and triangular['value'] != 4:
            raise ArithmeticError('n=2 exact value failed')
        if n == 3 and (triangular['points'] != [(0, 0), (-1, 0), (-1, 1)] or triangular['value'] != 27):
            raise ArithmeticError('n=3 tie handling or equilateral value failed')
        if n == 7 and triangular['value'] != centered_polygon(n):
            raise ArithmeticError('n=7 full first hexagonal shell failed')
        if n == 5 and square['value'] != centered_polygon(n):
            raise ArithmeticError('square first shell comparison failed')
        if n == 10 and improved <= Fraction(113, 100)*centered_polygon(n):
            raise ArithmeticError('old n=10 improvement comparison failed')
        rows.append({
            'n': n,
            'triangular': json_certificate(triangular),
            'gaussian_integer_disk': json_certificate(square),
            'prior_baseline': {name: encoded_fraction(value) for name, value in baseline.items()},
            'hermite_collinear': encoded_fraction(collinear),
            'triangular_strictly_improves_prior_envelope': stronger,
            'combined_envelope': encoded_fraction(combined_max),
            'combined_winners': winners,
        })
        table.append(
            f"| {n} | {decimal_string(triangular['value']/square['value'], 12)} "
            f"| {decimal_string(triangular['value']/prior_max, 12)} | {', '.join(winners)} |"
        )
    table.extend([
        '',
        'Within this finite range, the triangular rule strictly improves the prior envelope at n = '+', '.join(map(str, strict_improvements))+'.',
        '',
        'Special exact checks: T_2=4; T_3=27 with points (0,0),(-1,0),(-1,1); T_7=7^21/6^15. The square first shell gives its centered-square value at n=5. No claim of sharpness is made for n>=8.',
        '',
        'The prior envelope includes the Gaussian average, unit regular polygon, centered regular polygon, Gaussian-integer disk rule, and the already-proved explicit center improvement for n>=10. Hermite is recorded in JSON but is dominated by the regular polygon for every n>=3.',
        '',
        'The old n>=11 log improvement is converted to the rational lower bound C_n(1+x_m); it is a direct corollary of existing work, not a new perturbation theorem.',
        '',
        'No assertion about the exact linear asymptotic coefficient or global crystallization is used.',
        '',
        'Reproduce: python -B -S -O triangular_lattice.py --max-n 30',
    ])
    payload = {
        'format': 'mordell-triangular-lattice-exact-v1',
        'normalization': 'M_n=max_{sum|z_i|^2=n} product_{i<j}|z_i-z_j|^2',
        'omega': 'exp(i*pi/3)',
        'triangular_norm': 'a^2+a*b+b^2',
        'sorting': ['norm', 'a', 'b'],
        'finite_enumeration_box': '[-n,n]^2; contains global first n including ties',
        'proof': 'D=n*sum(q(v_i))-q(sum(v_i))=sum_{i<j}q(v_i-v_j)>0; centered scaled points n*(v_i-mean(v))/sqrt(D) have total squared modulus n',
        'max_n': max_n,
        'strict_improvements_in_recorded_range': strict_improvements,
        'rows': rows,
    }
    return payload, '\n'.join(table)+'\n'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-n', type=int, default=30)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    if args.max_n < 7:
        parser.error('--max-n must be at least 7 to include all special exact checks')
    payload, report = audit(args.max_n)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'triangular_lattice_values.json').write_text(
        json.dumps(payload, ensure_ascii=False, indent=2)+'\n', encoding='utf-8'
    )
    (args.output_dir/'triangular_lattice_audit.md').write_text(report, encoding='utf-8')
    print('Exact integer/rational lattice checks passed for n=2..'+str(args.max_n)+'.')
    print('Strict improvements over the prior envelope: '+', '.join(map(str, payload['strict_improvements_in_recorded_range'])))


if __name__ == '__main__':
    main()
