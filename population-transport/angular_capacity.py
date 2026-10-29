"""Exact angular population capacities from the author's stability operator.

No equal angular spacing is assumed. See ANGULAR_CAPACITY.md. These are
necessary conditions; passing them never certifies local/global maximality.
"""
from fractions import Fraction as F
from math import isqrt


def require_n(n):
    if type(n) is not int or n<2:
        raise ValueError('integer n>=2 required')


def population_cost(m):
    if type(m) is not int or m<1:
        raise ValueError('positive integer population required')
    return m*(2*(m//2))**2


def layer_capacity(n,population,radius_squared_upper):
    require_n(n)
    m=population
    if type(m) is not int or not 1<=m<=n:
        raise ValueError('integer population in [1,n] required')
    q=F(radius_squared_upper)
    if q<0:
        raise ValueError('nonnegative squared-radius upper required')
    even=2*(m//2)
    gap=even*even-8*(n-1)*q
    strict=(m>=3 and m%2==1)
    return {'n':n,'population':m,'radius_squared_upper':q,
            'minimum_radius_squared':F(even*even,8*(n-1)),
            'radius_lower_is_strict':strict,
            'excluded_from_local_maxima':gap>0 or (gap==0 and strict),
            'scope':'A single exact-radius layer; all angles allowed.'}


def population_certificate(populations):
    p=tuple(populations)
    if not p:
        raise ValueError('nonempty positive populations required')
    cost=sum(population_cost(m) for m in p)
    n=sum(p)
    require_n(n)
    cap=8*n*(n-1)
    strict=any(m>=3 and m%2 for m in p)
    return {'n':n,'layers':len(p),'cost':cost,'cost_cap':cap,
            'strict_cost_cap':strict,
            'excluded_from_local_maxima':cost>cap or (cost==cap and strict)}


def balanced_odds(total,count):
    if count==0:
        if total:
            raise ValueError('zero groups require zero population')
        return []
    if count<0 or total<count or (total-count)%2:
        raise ValueError('positive odd populations must have the required parity')
    small,extra=divmod((total-count)//2,count)
    m=1+2*small
    return [(m,count-extra)]+([(m+2,extra)] if extra else [])


def minimum_cost(n,layers,max_cases=1000000):
    """Exact minimum over ALL positive populations, with parity compression.

    At most one even population occurs in a minimizer. The remaining odd
    populations are balanced. One-even cases are fully enumerated, with
    refusal BEFORE incomplete enumeration if the budget is insufficient.
    """
    require_n(n)
    if type(layers) is not int or not 1<=layers<=n:
        raise ValueError('integer layers in [1,n] required')
    if type(max_cases) is not int or max_cases<1:
        raise ValueError('positive integer case budget required')
    def candidate(even):
        groups=balanced_odds(n-even,layers-int(even>0))
        if even:
            groups=[(even,1)]+groups
        groups=[(m,c) for m,c in groups if c]
        return {'cost':sum(c*population_cost(m) for m,c in groups),
                'population_groups':groups,
                'strict':any(m>=3 and m%2 and c for m,c in groups)}
    if (n-layers)%2==0:
        row=candidate(0)
        return {**row,'all_minimizers_strict':row['strict'],'cases':1}
    if layers==1:
        row=candidate(n)
        return {**row,'all_minimizers_strict':False,'cases':1}
    planned=(n-layers+1)//2
    if planned>max_cases:
        raise ValueError('complete unique-even population search exceeds case budget')
    best=None
    for even in range(2,n-layers+2,2):
        row=candidate(even)
        if best is None or row['cost']<best['cost']:
            best={**row,'all_minimizers_strict':row['strict']}
        elif row['cost']==best['cost']:
            best['all_minimizers_strict'] &= row['strict']
    return {**best,'cases':planned}


def minimum_layers(n,max_cases=1000000):
    require_n(n)
    cap=8*n*(n-1)
    def excluded(k):
        row=minimum_cost(n,k,max_cases)
        return row['cost']>cap or (row['cost']==cap and row['all_minimizers_strict'])
    lower,upper=0,n
    while lower+1<upper:
        mid=(lower+upper)//2
        if excluded(mid):lower=mid
        else:upper=mid
    return {'n':n,'minimum_layers':upper,'cost_cap':cap,
            'last_excluded_layers':lower,
            'last_excluded_cost':minimum_cost(n,lower,max_cases) if lower else None,
            'scope':'Necessary for all local maxima, using exact parity-dependent layer populations; not sufficient.'}


def elementary_layer_lower(n):
    """Rational endpoint below n/(1+2sqrt(2(n-1)))."""
    require_n(n)
    scale=10**12
    value=2*(n-1)*scale*scale
    r=isqrt(value)
    if r*r<value:r+=1
    return F(n)/(1+2*F(r,scale))


def chebyshev_value_derivative(degree,x):
    if type(degree) is not int or degree<0:
        raise ValueError('nonnegative integer polynomial degree required')
    x=F(x)
    if x<1:
        raise ValueError('this certificate requires x>=1')
    if degree==0:return F(1),F(0)
    a,da,b,db=F(1),F(0),x,F(1)
    for k in range(1,degree):
        a,da,b,db=b,db,2*x*b-a,2*b+2*x*db-da
    return b,db


def annulus_capacity(n,population,radius_lower,radius_upper,use_global_distance=False):
    """Entire thin radial band, arbitrary radii and angular placements.

    Uses local basic separation by default. Optional stronger separation
    restricts the conclusion to GLOBAL maxima. The width condition is
    essential to both positivity and the Jensen convexity proof.
    """
    require_n(n)
    m=population
    if type(m) is not int or not 2<=m<=n:
        raise ValueError('integer band population in [2,n] required')
    if type(use_global_distance) is not bool:
        raise ValueError('global-distance flag must be boolean')
    a,b=F(radius_lower),F(radius_upper)
    if not 0<=a<=b or b<=0:
        raise ValueError('ordered nonnegative radii with positive outer radius required')
    separation=F(4,n-1) if use_global_distance else F(2*n,(n-1)**2)
    w=b-a
    if 2*w*w>=separation:
        raise ValueError('band width must satisfy 2*width^2 < separation^2')
    even=2*(m//2)
    eta=1-2*w*w/separation
    x=1+w*w/(2*b*b)
    t,dt=chebyshev_value_derivative(even//2,x)
    reciprocal_sum=dt/(2*b*b*(t+1))
    lhs=eta*reciprocal_sum
    l=F(n-1,2)
    return {'n':n,'population':m,'even_subset_population':even,
            'radius_lower':a,'radius_upper':b,'squared_separation_lower':separation,
            'phase_real_factor_lower':eta,'odd_skip_reciprocal_sum_lower':reciprocal_sum,
            'operator_lower':lhs,'operator_upper':l,'stability_slack':l-lhs,
            'excluded':lhs>l,
            'scope':('All global maxima' if use_global_distance else 'All local maxima')+
                    ' with at least this many points in the stated thin annulus; no regular-angle assumption.'}
