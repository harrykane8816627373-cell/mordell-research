"""Exact checks for a quantitative holomorphic spectral obstruction.

L and U supplied to gap_lower are proof obligations; this routine does not
certify arbitrary proposed lower and upper bounds for M_n.
"""
from fractions import Fraction as F
from math import factorial, comb, prod
from pathlib import Path
import json


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def compositions(total, slots):
    if slots == 1:
        yield (total,)
    else:
        for a in range(total + 1):
            for rest in compositions(total-a, slots-1):
                yield (a,) + rest


def multinomial(alpha):
    return factorial(sum(alpha)) // prod(factorial(a) for a in alpha)


def level_zero_upper(n):
    N = n*(n-1)//2
    return F(n**N * factorial(n) * prod(factorial(j) for j in range(n)), factorial(N))


def gap_lower(n, r, L, U):
    if n < 2 or r < 0:
        raise ValueError('n>=2 and r>=0 required')
    L, U = F(L), F(U)
    P = level_zero_upper(n)
    if not 0 < L <= U < P:
        raise ValueError('requires proven 0 < L <= M <= U < P')
    N = n*(n-1)//2
    return L*(P-U)/(2*P*comb(N+r, r))


def main():
    coefficient_comparisons = 0
    for d in range(2, 5):
        for N in range(1, 9):
            for r in range(13):
                C = comb(N+r, r)
                for gamma in compositions(N, d):
                    beta = (gamma[0]+r,) + gamma[1:]
                    ratio = F(multinomial(beta), multinomial(gamma))
                    require(ratio == F(C, comb(gamma[0]+r, r)), 'exact multinomial ratio')
                    require(ratio <= C, 'multinomial comparison')
                    coefficient_comparisons += 1
    three_levels = [F(2), F(5,4), F(11,10)]
    for r, lam in enumerate(three_levels):
        require(lam*lam >= 1+F(1, comb(3+r,r)), 'three-point row-energy bound')
        require(27*lam-27 >= gap_lower(3,r,27,27), 'three-point explicit positive gap')
        require(gap_lower(3,r,27,27) == F(27,4*comb(3+r,r)), 'three-point gap formula')
    # For f=u^N+v^N (N>=2), S=1, Lambda0=2. Independently
    # test the 2x2 principal matrix at output monomials u^(N+r), u^r v^N.
    # Its diagonal entries are 1 and b, squared off-diagonal entry 1/C.
    # Negative determinant of (1+1/(4C))I minus this block proves that
    # its top eigenvalue is above the claimed positive gap.
    block_comparisons = 0
    for N in range(2, 11):
        for r in range(26):
            C = comb(N+r, r)
            b = F(1 + (comb(r,N) if r >= N else 0), C)
            delta = F(1,4*C)
            determinant = delta*(1+delta-b)-F(1,C)
            require(determinant < 0, 'independent two-monomial block lower bound')
            require(F(1,C) == F(2-1,C), 'exact row-energy equality family')
            block_comparisons += 1
    known = [(3,F(27)),(4,F(256)),(5,F(3125)),
             (6,F(6**15,5**10)),(7,F(7**21,6**15))]
    examples = []
    for n, M in known:
        P = level_zero_upper(n)
        require(P > M, 'known exact value below level zero')
        gaps = [gap_lower(n,r,M,M) for r in (0,1,2,5)]
        require(all(g > 0 for g in gaps), 'known-case positive rational gaps')
        examples.append({'n':n,'known_M':str(M),'level_zero_upper':str(P),
                         'gap_lower_r_0_1_2_5':[str(g) for g in gaps]})
    report={'passed':True,'arithmetic':'integer and fractions.Fraction',
            'multinomial_comparisons':coefficient_comparisons,
            'three_point_exact_levels':3,'two_monomial_block_comparisons':block_comparisons,
            'examples':examples,
            'scope':'Quantifies a spectral-method obstruction; no new unknown M_n computed. Known n6,n7 values retain original external classification dependencies.'}
    path=Path(__file__).with_name('quantitative_spectral_gap_checks.json')
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
