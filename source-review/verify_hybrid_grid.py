"""Finite exact checks for analytic hybrid constants, not grid enumeration."""
from fractions import Fraction as F
import json
from hybrid_grid_bounds import (ceil_root_ratio, root_upper, hybrid_exponents,
                                relative_plan)


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def audit():
    root_cases = 0
    for a in range(101):
        for b in range(1, 17):
            for degree in (2,4):
                r = ceil_root_ratio(a,b,degree)
                require(r**degree*b >= a and (r==0 or (r-1)**degree*b < a), 'integer root cap')
                root_cases += 1
    exponent_cases = 0
    for n in range(3, 101):
        for q in (n, 2*n, 10*n, 1000000*n):
            row = hybrid_exponents(n,q)
            t = row['tau']
            require(2*q*q*t*t >= n-1 and 0<t<1, 'new pair enclosure')
            require(n*q*q*t*t < (n-1)**2, 'strictly below old b')
            require(row['root4_n_upper']**4 >= n, 'fourth root enclosure')
            require(row['root2_nm1_upper']**2 >= n-1, 'square root enclosure')
            require(row['best'] == min(row['candidates'].values()), 'minimum of valid bounds')
            require(row['best'] <= row['candidates']['existing_uniform_row_majorant'], 'old fallback')
            exponent_cases += 1
    plans = []
    for n in (2,3,5,7,8,10,30,100):
        for eps in (F(1,10), F(1,100)):
            row = relative_plan(n, eps)
            require(row['eta'] <= eps/(1+eps), 'relative exponent')
            if n>=3:
                require(row['relative_factor_upper'] <= 1+eps, 'rational stopping factor')
                if 2*(row['q']-1)**2 > n-1:
                    require(hybrid_exponents(n,row['q']-1)['best'] > eps/(1+eps), 'minimal sufficient q')
            plans.append({'n':n,'epsilon':str(eps),'q':row['q']})
    for q in (1,2,10,1000000):
        require(hybrid_exponents(2,q)['best']==0, 'two-point exact')
    # Recheck the rational root facts used by the existing 35/72 bound.
    require(8**4 < 2*7**4 and 25<27, 'existing rational majorants')
    return {'passed':True,'integer_root_cases':root_cases,
            'hybrid_exponent_cases':exponent_cases,'plans':plans,
            'scope':'Exact finite checks of constants and schedules; no new complete grid enumeration.'}


if __name__ == '__main__':
    print(json.dumps(audit(),indent=2))
