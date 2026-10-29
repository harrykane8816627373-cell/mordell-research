"""Layer population constraints from the author's strict subset separation.

No regular angular spacing is assumed. Proof: RING_POPULATION.md.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt


@lru_cache(maxsize=1)
def pi_interval():
    """Machin identity with outward alternating rational arctangent sums."""
    def atan(q, terms=16):
        total=sum(((-1)**j*q**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
        next_term=(-1)**terms*q**(2*terms+1)/F(2*terms+1)
        return min(total,total+next_term),max(total,total+next_term)
    a,b=atan(F(1,5)),atan(F(1,239))
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]


def population_cost(populations):
    populations=tuple(populations)
    if not populations or any(type(m) is not int or m<1 for m in populations):
        raise ValueError('nonempty positive integer populations required')
    return sum(m**3 for m in populations if m>=2)


def minimum_cost(n,layers,max_cases=1000000):
    """Complete minimum over ALL positive integer populations of n in K layers.

    Condition on the number of singleton layers, then balance other parts.
    Budget refusal occurs before enumerating any singleton case.
    """
    if type(n) is not int or n<2 or type(layers) is not int or not 1<=layers<=n:
        raise ValueError('integer n>=2 and 1<=layers<=n required')
    if type(max_cases) is not int or max_cases<1:
        raise ValueError('positive integer case budget required')
    if layers==n:
        return {'cost':0,'singletons':n,'balanced_small':0,'small_count':0,
                'large_count':0,'cases':1}
    first=max(0,2*layers-n)
    planned=layers-first
    if planned>max_cases:
        raise ValueError('complete singleton allocation exceeds case budget')
    best=None
    for singles in range(first,layers):
        groups=layers-singles
        points=n-singles
        small,large_count=divmod(points,groups)
        cost=(groups-large_count)*small**3+large_count*(small+1)**3
        row={'cost':cost,'singletons':singles,'balanced_small':small,
             'small_count':groups-large_count,'large_count':large_count,'cases':planned}
        if best is None or cost<best['cost']:
            best=row
    return best


def minimum_layers(n,max_cases=1000000):
    """Necessary layers for a GLOBAL optimizer, using exact integer costs."""
    if type(n) is not int or n<2:
        raise ValueError('integer n>=2 required')
    cap=pi_interval()[1]**2*n*(n-1)
    upper=min(n,isqrt(n)+1)
    while minimum_cost(n,upper,max_cases)['cost']>=cap:
        upper=min(n,2*upper)
    lower=0
    while lower+1<upper:
        middle=(lower+upper)//2
        if minimum_cost(n,middle,max_cases)['cost']>=cap:
            lower=middle
        else:
            upper=middle
    last=minimum_cost(n,lower,max_cases) if lower else None
    return {'n':n,'minimum_layers':upper,'last_excluded_layers':lower,
            'last_excluded_cost':last,'cost_cap_upper':cap,
            'scope':'Necessary global-maximizer layer count from strict subset separation and angular circle packing; no sufficiency or exact optimizer claimed.'}


def population_certificate(populations):
    populations=tuple(populations)
    cost=population_cost(populations)
    n=sum(populations)
    if n<2:
        raise ValueError('at least two total points required')
    floor=F(cost)/((n-1)*pi_interval()[1]**2)
    return {'n':n,'cost':cost,'moment_lower':floor,
            'excluded_from_global_maxima':floor>=n}
