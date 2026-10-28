"""Exact constants and elementary identities for RADIAL_LAYERS.md.

The all-n electrostatic argument is analytic; these finite checks do not
replace its Dirichlet representation or rotational averaging proof.
"""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import json


def require(ok, name):
    if not ok:
        raise ArithmeticError(name)


def log_bounds(x, terms=64):
    x = F(x)
    if x <= 0:
        raise ValueError('positive logarithm argument required')
    if x < 1:
        a, b = log_bounds(1/x, terms)
        return -b, -a
    t = (x-1)/(x+1)
    low = 2*sum((t**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    tail = 2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
    return low, low+tail


def exp_upper(x, terms=40):
    x = F(x)
    require(x >= 0 and x < terms+2, 'exponential geometric tail')
    value = sum((x**j/factorial(j) for j in range(terms+1)), F(0))
    return value + x**(terms+1)/factorial(terms+1)/(1-x/(terms+2))


def radial_gap_lower(n, K, root_scale=1000000):
    """A rational lower bound on the finite radial-energy gap expression.

    This necessary bound alone does not certify the existence of a
    maximizing configuration with K radii.
    """
    if type(n) is not int or n < 2 or type(K) is not int or K < 1:
        raise ValueError('integer n>=2, K>=1 required')
    if type(root_scale) is not int or root_scale < 1:
        raise ValueError('positive integer root scale required')
    root_lower = F(isqrt(n*root_scale*root_scale), root_scale)
    length_lower = max(F(0), 2-F(8*K)/root_lower)
    return F(n*n, 192*(K+1)**2)*length_lower**3


def main():
    # Gamma(1/3)=3 int_0^infty exp(-u^3)du. The odd Taylor
    # polynomial is a pointwise lower bound on [0,3].
    gamma_integral_lower = 3*sum((F((-1)**j*3**(3*j+1), (3*j+1)*factorial(j))
                                  for j in range(152)), F(0))
    gamma_lower = F(26789385, 10000000)
    require(gamma_integral_lower > gamma_lower, 'Gamma(1/3) lower certificate')
    def atan_partial(x, last):
        return sum((F((-1)**j)*x**(2*j+1)/(2*j+1) for j in range(last+1)), F(0))
    pi_upper = F(355,113)
    machin_upper = 16*atan_partial(F(1,5),10)-4*atan_partial(F(1,239),9)
    require(machin_upper < pi_upper, 'Machin pi upper certificate')
    ratio_upper = 8*exp_upper(F(3,2))*pi_upper**3/(3*gamma_lower**6)
    require(ratio_upper < F(401,400), 'gamma_D - gamma_triangle < 1/800')
    # 401/400 < exp(1/400), so log(ratio_upper)/2 < 1/800.
    log12_lower = 2*log_bounds(2)[0]+log_bounds(3)[0]
    require(log12_lower > F(62,25), 'gamma_D - gamma_polar < 1/100')
    require(F(1,3) <= F(49,144), 'sqrt(n/3) <= 7 sqrt(n)/12')
    finite_upper = F(1,100)+F(77,96)+F(35,1024)
    require(finite_upper < F(9,8), 'all-n finite K > sqrt(n)/16 contradiction')
    enhanced_upper = F(1,100)+F(77,144)+F(35,2304)
    require(enhanced_upper < F(64,81), 'enhanced finite K > sqrt(n)/12 contradiction')
    c = F(11,50)
    f = (2-8*c)**3/(192*c*c)
    require(f == F(9,6050) and f > F(1,800), 'asymptotic rational constant')
    # Exact integral over one radial gap, including arbitrary CDF values.
    gap_checks = 0
    for a in (F(0), F(1,7), F(3,5), F(1)):
        for length in (F(1,10), F(1,4), F(1,2)):
            b = a+length
            if b > 2:
                continue
            optimum = (a+b)/4
            for offset in (F(-1,3), F(0), F(2,5)):
                value = optimum+offset
                integral = value*value*(b-a)-value*(b*b-a*a)/2+(b**3-a**3)/12
                require(integral == length**3/48+length*offset*offset, 'gap integration identity')
                require(integral/4 >= length**3/192, 'gap energy lower bound')
                gap_checks += 1
    width_checks = 0
    for a in range(21):
        R = F(a,20)
        for b in range(1,21):
            r = F(b,20)
            width = (R+r)**2-max(R-r, F(0))**2
            require(width <= 4*r, 'unit-radius support width')
            width_checks += 1
    for n in range(2,101):
        for K in range(1,11):
            require(radial_gap_lower(n,K) >= 0, 'rational radial gap nonnegative')
    report = {'passed':True,'arithmetic':'integers and Fraction; no floating decisions',
              'gamma_lower':str(gamma_lower),'pi_upper':str(pi_upper),
              'gamma_gap_strictly_below':'1/800',
              'finite_K_factor':'1/16', 'enhanced_finite_K_factor':'1/12',
              'asymptotic_K_factor':'11/50',
              'gap_integral_cases':gap_checks,'support_width_cases':width_checks,
              'scope':'Constants and finite identity checks; RADIAL_LAYERS.md proves the all-n energy argument.'}
    (Path(__file__).parent/'radial_layers_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
