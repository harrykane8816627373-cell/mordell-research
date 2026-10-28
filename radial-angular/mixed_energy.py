"""Exact radial-weighted angular Gram certificates, beyond holomorphic moments.

Proof and scope: MIXED_ENERGY.md. Fraction-only acceptance, no optimizer.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'joint-bounds'))
from angular_bounds import sqrt_upper, log_bounds, cmul


def require_index(k, degree):
    if type(k) is not int or k < 1 or type(degree) is not int or degree < 0:
        raise ValueError('integer angular k>=1 and radial degree>=0 required')


def gram_matrix(k, degree, support_squared):
    require_index(k, degree)
    s = F(support_squared)
    if s <= 0:
        raise ValueError('positive squared support required')
    return [[F((k+j)*(k+l), k+j+l)*s**(k+j+l)
             for l in range(degree+1)] for j in range(degree+1)]


def solve(matrix, rhs):
    """Exact symmetric positive-definite solve, rejecting nonpositive pivots."""
    n = len(matrix)
    if n == 0 or len(rhs) != n or any(len(row) != n for row in matrix):
        raise ValueError('nonempty square matrix and matching vector required')
    a = [list(map(F, row))+[F(b)] for row, b in zip(matrix, rhs)]
    for j in range(n):
        if a[j][j] <= 0:
            raise ValueError('positive Gram pivot required')
        pivot = a[j][j]
        for i in range(j+1, n):
            ratio = a[i][j]/pivot
            for l in range(j, n+1):
                a[i][l] -= ratio*a[j][l]
    x = [F(0)]*n
    for j in range(n-1, -1, -1):
        x[j] = (a[j][n]-sum((a[j][l]*x[l] for l in range(j+1, n)), F(0)))/a[j][j]
    return x


def projection_value(k, support_squared, real_moments, imaginary_moments):
    """A*H^-1 A for ACTUAL smeared mixed moments A, not divided by n."""
    if not real_moments or len(real_moments) != len(imaginary_moments):
        raise ValueError('matching nonempty real and imaginary moment lists required')
    h = gram_matrix(k, len(real_moments)-1, support_squared)
    re, im = list(map(F, real_moments)), list(map(F, imaginary_moments))
    ar, ai = solve(h, re), solve(h, im)
    return sum((a*b+c*d for a, b, c, d in zip(re, ar, im, ai)), F(0))


def dual_interval_bound(k, support_squared, coefficients, moment_rectangles):
    """Entire-class angular bound with a fixed complex rational test vector.

    Each moment rectangle is (Re lower, Re upper, Im lower, Im upper) for
    the actual A_kj. Independent proofs of these intervals and support are
    caller hypotheses. Bounds for different k may be summed; same-k tests
    must be put into one Gram certificate or combined by taking a maximum.
    """
    coefficients = [(F(a), F(b)) for a, b in coefficients]
    if not coefficients or len(coefficients) != len(moment_rectangles):
        raise ValueError('matching nonempty coefficients and rectangles required')
    h = gram_matrix(k, len(coefficients)-1, support_squared)
    pairing = F(0)
    for (a, b), rect in zip(coefficients, moment_rectangles):
        lo, hi, ilo, ihi = map(F, rect)
        if lo > hi or ilo > ihi:
            raise ValueError('ordered moment intervals required')
        pairing += min(a*lo, a*hi)+min(b*ilo, b*ihi)
    norm = sum((h[j][l]*(coefficients[j][0]*coefficients[l][0]
                        + coefficients[j][1]*coefficients[l][1])
                for j in range(len(h)) for l in range(len(h))), F(0))
    return max(F(0), 2*pairing-norm)


def disk_weight(k, j, q, smear_squared):
    """Disk mean of w^k |w|^(2j) equals center^k times this radial weight."""
    require_index(k, j)
    q, r2 = F(q), F(smear_squared)
    if q < 0 or r2 <= 0:
        raise ValueError('nonnegative center norm and positive smear required')
    return sum((F(comb(k+j, l)*comb(j, l), l+1)*q**(j-l)*r2**l
                for l in range(j+1)), F(0))


def shape_certificate(points, max_mode=3, radial_degree=1, h=1,
                      support_squared=None, radial_delta_lower=0):
    """Normalize centered rational points to Q(x)=n-1 and project all modes.

    radial_delta_lower is a separately proved lower bound for n^2 D(bar nu)
    using this SAME smear h*sqrt(2/n). Zero is always valid. An externally
    supplied support must contain every smeared point disk.
    """
    require_index(max_mode, radial_degree)
    pts = [(F(a), F(b)) for a, b in points]
    n = len(pts)
    if n < 2 or len(set(pts)) != n:
        raise ValueError('at least two distinct rational points required')
    h, radial = F(h), F(radial_delta_lower)
    if h <= 0 or radial < 0:
        raise ValueError('positive smear factor and nonnegative radial lower bound required')
    mx = sum((a for a, b in pts), F(0))/n
    my = sum((b for a, b in pts), F(0))/n
    v = [(a-mx, b-my) for a, b in pts]
    qraw = [a*a+b*b for a, b in v]
    factor = F(n-1)/sum(qraw)
    qs = [factor*q for q in qraw]
    r2 = F(2, n)*h*h
    qmax = max(qs)
    if support_squared is None:
        s = qmax+r2+2*sqrt_upper(qmax*r2)
    else:
        s = F(support_squared)
        c = s-qmax-r2
        if c < 0 or c*c < 4*qmax*r2:
            raise ValueError('supplied support does not contain smeared points')
    rows, total, old = [], F(0), F(0)
    powers = [(F(1), F(0)) for _ in v]
    for k in range(1, max_mode+1):
        powers = [cmul(a, b) for a, b in zip(powers, v)]
        raw_re, raw_im = [], []
        for j in range(radial_degree+1):
            weights = [disk_weight(k, j, q, r2) for q in qs]
            raw_re.append(sum((p[0]*w for p, w in zip(powers, weights)), F(0)))
            raw_im.append(sum((p[1]*w for p, w in zip(powers, weights)), F(0)))
        # All moments share factor^(k/2); its square is rational.
        value = factor**k*projection_value(k, s, raw_re, raw_im)
        base = factor**k*(raw_re[0]**2+raw_im[0]**2)/(k*s**k)
        if value < base:
            raise ArithmeticError('nested projection must dominate the original moment')
        total += value
        old += base
        rows.append({'k': k, 'projection': value, 'holomorphic_only': base,
                     'raw_real_moments': raw_re, 'raw_imaginary_moments': raw_im})
    # L(r)-L(r_opt)=-n psi(h). Upper-bound psi when subtracting it.
    penalty_upper = n*((h*h-1)/2-log_bounds(h)[0])
    return {'n': n, 'radial_degree': radial_degree, 'max_mode': max_mode,
            'normalization_factor_squared': factor, 'smear_factor': h,
            'support_squared': s, 'smear_squared': r2,
            'angular_delta_lower': total, 'holomorphic_delta_lower': old,
            'radial_delta_lower': radial, 'smear_penalty_upper': penalty_upper,
            'delta_lower': max(F(0), total+radial-penalty_upper), 'modes': rows,
            'scope': 'Necessary energy certificate; caller must prove any supplied nonzero radial lower bound. No optimality or full-domain cover.'}


if __name__ == '__main__':
    import json
    example = [(3, 0), (-3, 0), (4, 0), (-4, 0), (0, 5), (0, -5)]
    print(json.dumps(shape_certificate(example, support_squared=3), default=str, indent=2))
