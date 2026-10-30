"""Reusable exact six-neighbor/vertex-cover/quiet-packing certificates."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import sys
from overlap_cost import certified_pi_upper,pair_cost_squared_lower

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'local-defect'))
from packing_constants import packing_less_than
from linear_gap import root_up
from subset_moment import log_interval,outward_decimal


def validate(c,t,b,endpoint_cap,quiet_cap,kappa,boundary=11,cosine_lower=F(3117,5000)):
    c,t,b,kappa,boundary,q=map(F,(c,t,b,kappa,boundary,cosine_lower))
    if c<=0 or c*c<2 or c>=2 or b<=0 or t<=1+b or kappa<=0 or boundary<=0 or not 0<q<1:
        raise ValueError('positive geometric parameters, sqrt(2)<=c<2 and t>1+b required')
    if type(endpoint_cap)is not int or endpoint_cap<1 or type(quiet_cap)is not int or quiet_cap<1:
        raise ValueError('positive integer packing caps required')
    angle=2*certified_pi_upper()/7
    cosine=sum(((-1)**j*angle**(2*j)/factorial(2*j) for j in range(6)),F(0))
    if not angle<1 or not cosine>q:raise ValueError('unproved seven-neighbor cosine comparison')
    U=t+1+b;L=2*(t+b)
    corners=[x*x+y*y-2*q*x*y for x in (c,U) for y in (c,U)]
    if max(corners)>=c*c:raise ValueError('seven neighbors not excluded')
    beta=(c*c+t*t-(1+b)**2)/(2*c*t)
    if c*c<=t*t-(1+b)**2 or beta<=0 or 4*beta*beta<=3:
        raise ValueError('six inflated arcs may cover the query circle')
    if not packing_less_than(U,endpoint_cap+1):raise ValueError('endpoint packing cap too small')
    if (1+2*L/c)**2>=quiet_cap+1:raise ValueError('quiet c-separated packing cap too small')
    a=t+b;margin=boundary/4-a
    if margin<=0 or margin*margin<=2 or 2-a<=0 or (2-a)**2<=F(1,2):
        raise ValueError('all-n boundary annulus bound not proved')
    pair_squared=pair_cost_squared_lower(c)
    if pair_squared<(endpoint_cap*kappa)**2:raise ValueError('pair cost cannot pay vertex-cover population')
    hole=b**4/4
    if hole<quiet_cap*kappa:raise ValueError('hole cost cannot pay quiet population')
    return {'close_cutoff':c,'query_radius':t,'hole_radius':b,'near_radius':U,'greedy_radius':L,
        'endpoint_cap':endpoint_cap,'quiet_cap':quiet_cap,'boundary_coefficient':boundary,
        'linear_coefficient':kappa,'pair_cost_used':endpoint_cap*kappa,'hole_cost':hole,
        'seventh_neighbor_cosine_lower':q,'cap_cosine_lower':beta,
        'parabolic_pair_squared_lower':pair_squared,
        'old_coefficient_ratio':kappa/F(1,2**24),
        'scope':'For every global optimizer, all n>=2: deficit>=kappa*(n-boundary*sqrt(n))_+. Generic parameters are accepted only after every exact comparison succeeds.'}


FAMILY={
    'simple32':(F(191,100),F(11001,10000),F(1,10),16,12,F(1,2**19)),
    'strong37':(F(1913,1000),F(11017,10000),F(127,1250),16,12,F(37,2**24))
}


@lru_cache(maxsize=2)
def family_certificate(name='strong37'):
    if name not in FAMILY:raise ValueError('unknown fixed certificate')
    return validate(*FAMILY[name])


def gap_lower(n,name='strong37',scale=10**12):
    if type(n)is not int or n<2:raise ValueError('integer n>=2 required')
    cert=family_certificate(name)
    return cert['linear_coefficient']*max(F(0),n-cert['boundary_coefficient']*root_up(n,scale))


def gamma_upper_interval(name='strong37',places=30):
    kappa=family_certificate(name)['linear_coefficient']
    lo,hi=log_interval(2,80)
    return outward_decimal((F(3,4)-hi/2-kappa,F(3,4)-lo/2-kappa),places)
