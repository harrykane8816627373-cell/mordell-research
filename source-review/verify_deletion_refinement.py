"""Exact coefficient and logarithm-enclosure checks for deletion refinement."""
from fractions import Fraction as F
import json


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def log_bounds(x, terms=40):
    x = F(x)
    if x <= 0:
        raise ValueError('positive logarithm argument required')
    y = (x-1)/(x+1)
    partial = sum((2*y**(2*j+1)/(2*j+1) for j in range(terms)), F(0))
    tail = 2*abs(y)**(2*terms+1)/((2*terms+1)*(1-y*y))
    return partial-tail, partial+tail


def run():
    rows = []
    for n in range(3, 101):
        mu = F(n*(n-2), n-1)
        T = F(n,1)/mu
        require(T == F(n-1,n-2), 'deleted mean ratio')
        N = F(n*(n-1),2)
        variance_scale = F(n,(n-1)**2)
        prefactor = N*variance_scale/(mu*mu)
        require(prefactor == F(n-1,2*(n-2)**2), 'variance penalty rescaling')
        require(prefactor/(T-1)**2 == F(n-1,2), 'h(T) logarithm factor')
        lo, hi = log_bounds(T)
        c_lower = F(n-1,2)*(F(1,n-2)-hi)
        c_upper = F(n-1,2)*(F(1,n-2)-lo)
        require(c_lower > F(1,4*(n-2)) > F(1,4*(n-1)), 'strict coefficient comparison')
        rows.append({'n':n,'c_lower':str(c_lower),'c_upper':str(c_upper)})
    # Rational two-point mean-one distributions used in the sharpness proof.
    sharpness_cases = 0
    for T in (F(3,2), F(2), F(3), F(5)):
        for p in (F(1,100),F(1,1000)):
            x = (1-p*T)/(1-p)
            require(x > 0 and p*T+(1-p)*x == 1, 'positive sharpness distribution')
            direct_var = p*(T-1)**2+(1-p)*(x-1)**2
            require(direct_var == p*(T-1)**2/(1-p), 'sharpness variance expression')
            sharpness_cases += 1
    return {'passed':True,'coefficient_cases':len(rows),
            'sharpness_distribution_cases':sharpness_cases,
            'arithmetic':'Fractions and explicit logarithm-tail intervals',
            'scope':'Finite coefficient checks; DELETION_REFINEMENT.md proves the general theorem.',
            'sample_enclosures':[rows[0],rows[3],rows[-1]]}


if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
