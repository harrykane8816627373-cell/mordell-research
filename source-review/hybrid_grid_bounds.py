"""Combine an existing row-norm grid estimate with later strict pair bounds.

Only integer/Fraction arithmetic is used. These are sufficient analytic
error factors for COMPLETE specified grids, not sampled configurations.
"""
from fractions import Fraction as F
from math import isqrt


def require_nq(n, q):
    if type(n) is not int or n < 2 or type(q) is not int or q < 1:
        raise ValueError('n>=2 and q>=1 must be integers')


def ceil_root_ratio(a, b, degree):
    if a < 0 or b <= 0 or degree not in (2, 4):
        raise ValueError('invalid root arguments')
    floor = isqrt(a // b)
    if degree == 4:
        floor = isqrt(floor)
    return floor + int(floor ** degree * b < a)


def root_upper(a, degree=2, scale=65536):
    if a < 0 or type(scale) is not int or scale < 1:
        raise ValueError('invalid root enclosure')
    return F(ceil_root_ratio(a * scale ** degree, 1, degree), scale)


def hybrid_exponents(n, q):
    """Return positive rational upper exponents for the global grid error.

    The mathematical row bound is an existing result of the revised
    source; its denominator is improved with the later strict pair gap.
    Every term in `candidates` is a valid global logarithmic error upper
    bound under the declared hypotheses. Their minimum remains valid.
    """
    require_nq(n, q)
    if n == 2:
        return {'n': n, 'q': q, 'tau': F(0), 'leading': F(0),
                'candidates': {'two_point_exact': F(0)}, 'best': F(0)}
    if 2*q*q <= n-1:
        raise ValueError('need 2*q*q > n-1 for collision-free rounding')
    first_admissible_q = isqrt((n-1)//2) + 1
    scale = 65536
    while True:
        # Enclose the constant BEFORE division by q; precision remains
        # meaningful even when q is very large.
        c = F(ceil_root_ratio((n-1)*scale*scale, 2, 2), scale)
        tau = c/q
        # For n>2, the true new pair-gap constant is strictly below the
        # old constant (n-1)/sqrt(n). Retain this useful strict comparison.
        # Make the chosen precision depend only on n, not on q. This
        # also makes the resulting rational error bound monotone in q.
        if c < first_admissible_q and n*c*c < (n-1)**2:
            break
        scale *= 2
    leading = F(n*(n-1), 2*q*q)
    pair = leading * (1 + F(n, 6)*tau/(1-tau))
    root4_n = root_upper(n, 4, scale)
    root2_nm1 = root_upper(n-1, 2, scale)
    row = leading + root4_n*(n-1)**2*root2_nm1/(3*(1-tau)*q**3)
    candidates = {'strict_pair_sum': pair, 'old_row_norm_new_pair_gap': row}
    if q >= n:
        # Existing revised-source rational majorant. It is retained as a
        # fallback, independently of the outward root precision above.
        candidates['existing_uniform_row_majorant'] = leading + F(35*n**3, 72*q**3)
    return {'n': n, 'q': q, 'tau': tau, 'leading': leading,
            'root4_n_upper':root4_n, 'root2_nm1_upper':root2_nm1,
            'candidates':candidates, 'best':min(candidates.values())}


def relative_plan(n, epsilon):
    """Smallest q certified by this explicit monotone rational bound.

    This plans a COMPLETE grid; it does not execute the grid itself.
    """
    require_nq(n, 1)
    epsilon = F(epsilon)
    if epsilon <= 0:
        raise ValueError('positive relative error required')
    if n == 2:
        return {'n':2, 'q':1, 'epsilon':epsilon, 'eta':F(0), 'search_executed':False}
    target = epsilon/(1+epsilon)
    lower = isqrt((n-1)//2) + 1
    # The existing 35/72 uniform bound guarantees this initial high q
    # is sufficient when epsilon/(1+epsilon)<1.
    upper = max(n, ceil_root_ratio(n*n*(1+epsilon).numerator*epsilon.denominator,
                                  (1+epsilon).denominator*epsilon.numerator, 2))
    while hybrid_exponents(n, upper)['best'] > target:
        upper *= 2
    while lower < upper:
        mid = (lower+upper)//2
        if hybrid_exponents(n, mid)['best'] <= target:
            upper = mid
        else:
            lower = mid+1
    eta = hybrid_exponents(n, lower)['best']
    if not eta <= target < 1:
        raise ArithmeticError('relative plan failed')
    return {'n':n, 'q':lower, 'epsilon':epsilon, 'eta':eta,
            'relative_factor_upper':1/(1-eta), 'search_executed':False}


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('n', type=int)
    p.add_argument('--q', type=int)
    p.add_argument('--epsilon', default='1/100')
    args = p.parse_args()
    result = hybrid_exponents(args.n, args.q) if args.q is not None else relative_plan(args.n, args.epsilon)
    print(json.dumps(result, default=str, ensure_ascii=False, indent=2))
