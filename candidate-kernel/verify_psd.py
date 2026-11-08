"""Exact PSD engineering checks; no continuous logarithmic majorant claimed.

The frozen n=10 blocks are repaired by adding 2^-22 I, with a uniform
perturbation cap from exact univariate Bernstein coefficients. The frozen
n=12 blocks are strictly positive definite without any repair. These facts
do not remove either saved counterexample to the two-point majorant.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import json

HERE = Path(__file__).resolve().parent
U = F(15403, 10000)
SHIFT = F(1, 2**22)


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def decode(q):
    return [decode(v) for v in q] if isinstance(q, list) else F.from_float(float.fromhex(q))


def positive_ldl(Q):
    """Prove strict positivity by an exact reconstructed LDL^T factorization."""
    n = len(Q)
    need(n > 0 and all(len(row) == n for row in Q), 'nonempty square block')
    need(all(Q[i][j] == Q[j][i] for i in range(n) for j in range(n)), 'exact symmetry')
    L = [[F(i == j) for j in range(n)] for i in range(n)]
    D = []
    for i in range(n):
        pivot = Q[i][i] - sum(L[i][j]**2 * D[j] for j in range(i))
        need(pivot > 0, 'positive exact LDL pivot')
        D.append(pivot)
        for k in range(i+1, n):
            L[k][i] = (Q[k][i] - sum(L[k][j]*L[i][j]*D[j] for j in range(i))) / pivot
    need(all(Q[i][j] == sum(L[i][k]*D[k]*L[j][k] for k in range(n))
             for i in range(n) for j in range(n)), 'exact LDL reconstruction')
    return D


def block_report(Qs, shift):
    rows = []
    for k, Q in enumerate(Qs):
        repaired = [[a + shift*(i == j) for j, a in enumerate(row)] for i, row in enumerate(Q)]
        D = positive_ldl(repaired)
        unit = 2**100
        scaled = min(D)*unit
        lower = F(scaled.numerator//scaled.denominator - 1, unit)
        need(lower > 0, 'displayed positive pivot lower bound')
        rows.append({'k': k, 'dimension': len(Q), 'all_exact_LDL_pivots_positive': True,
                     'exact_LDL_reconstruction_checked': True,
                     'all_pivots_strictly_above_dyadic_lower': str(lower)})
    return rows


def bernstein_absolute(coefficients):
    """Exact polynomial absolute bound on [0,U], by Bernstein convex hull.

    If p(u)=sum c_i*u^i has degree at most d, its Bernstein coefficient of
    index j on [0,U] is sum_{i<=j} c_i*U^i*C(j,i)/C(d,i). The Bernstein
    basis is nonnegative and sums to one throughout this entire interval.
    """
    d = len(coefficients)-1
    b = [sum((coefficients[i]*U**i*F(comb(j, i), comb(d, i)) for i in range(j+1)), F(0))
         for j in range(d+1)]
    return max(map(abs, b))


def multiply(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def ten_point_perturbation_cap(data):
    nodes = decode(data['radial_nodes'])
    means = decode(data['radial_moment_centers'])
    # For a physical pair, |H_k| <= (uv)^(k/2) <= U^k. Let A_k bound
    # the squared Euclidean norm of the corresponding radial features.
    # Then |Delta K| <= SHIFT*sum(U^k*A_k) by Cauchy--Schwarz.
    radial = [F(1)]
    for a in nodes:
        radial = multiply(radial, [-a, F(1)])
    radial = [10*x for x in multiply(radial, radial)]
    a0 = sum(max(abs(mean), abs(U**j-mean))**2 for j, mean in enumerate(means, 1))
    a0 += bernstein_absolute(radial)**2
    A = [a0]
    for transform in decode(data['transforms']):
        s = len(transform)
        cols = len(transform[0])
        norm_columns = F(0)
        for j in range(cols):
            # E acts on ((u-1)^i); expand each projected feature in u.
            coefficients = [sum(transform[i][j]*comb(i, r)*(-1)**(i-r) for i in range(r, s))
                            for r in range(s)]
            norm_columns += bernstein_absolute(coefficients)**2
        gram = [[sum(transform[i][j]*transform[l][j] for j in range(cols))
                 for l in range(s)] for i in range(s)]
        # On [0,U] subset [0,2], each |u-1|^i <= 1. The absolute row sum
        # bounds the spectral norm of this real symmetric Gram matrix.
        need(U < 2, 'unit bound on shifted radial powers')
        gram_bound = s*max(sum(map(abs, row)) for row in gram)
        A.append(min(norm_columns, gram_bound))
    kernel_cap = SHIFT*sum(U**k*a for k, a in enumerate(A))
    # h=p(u)+p(v)-K and Delta p=-Delta K(u,u,u)/(2*9).
    h_cap = F(10, 9)*kernel_cap
    need(h_cap < F(231, 10**6), 'small uniform repair cost')
    need(45*h_cap < F(104, 10000), 'total repair budget')
    return {'domain': '0<=u,v<=15403/10000; t^2<=uv',
            'radial_feature_squared_norm_caps': [str(a) for a in A],
            'absolute_kernel_change_upper': str(kernel_cap),
            'absolute_h_change_upper': str(h_cap),
            'absolute_h_change_strictly_less_than': '231/1000000',
            '45_times_h_change_upper': str(45*h_cap),
            '45_times_h_change_strictly_less_than': '104/10000'}


def verify():
    paths = {n: HERE/f'data/n{n}_exact_interpretation.json' for n in (10, 12)}
    raw = {n: path.read_bytes() for n, path in paths.items()}
    data = {n: json.loads(blob) for n, blob in raw.items()}
    ten, twelve = data[10], data[12]
    need((ten['n'], ten['mass'], ten['feature_order'], ten['harmonics_through']) == (10, 9, 6, 9),
         'fixed ten-point fitted kernel')
    need((twelve['n'], twelve['mass'], twelve['feature_order'], twelve['harmonics_through']) == (12, 11, 5, 11),
         'fixed twelve-point fitted kernel')
    Q10 = decode(ten['Q_blocks'])
    Q12 = decode(twelve['Q_blocks'])
    need(len(Q10) == 10 and len(Q12) == 12, 'all fixed harmonic blocks present')
    repaired = block_report(Q10, SHIFT)
    original = block_report(Q12, F(0))
    return {'passed': True,
            'arithmetic': 'exact Fraction; exact LDL reconstruction; univariate Bernstein convex hull',
            'source_sha256': {paths[n].relative_to(HERE).as_posix(): hashlib.sha256(raw[n]).hexdigest()
                              for n in (10, 12)},
            'n10': {'diagonal_shift_per_Q_block': '2^-22',
                    'all_repaired_Q_blocks_strictly_positive_definite': True,
                    'blocks': repaired, 'uniform_perturbation': ten_point_perturbation_cap(ten),
                    'continuous_log_majorant_accepted': False},
            'n12': {'diagonal_shift_per_Q_block': '0',
                    'all_original_Q_blocks_strictly_positive_definite': True,
                    'blocks': original, 'continuous_log_majorant_accepted': False},
            'continuous_log_majorant_bound_proved': False,
            'global_upper_bound_proved': False}


def main():
    print(json.dumps(verify(), indent=2))


if __name__ == '__main__':
    main()
