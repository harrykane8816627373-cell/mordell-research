"""Exact constants for the close-pair or empty-hole local alternative.

Geometry uses force-one w=sqrt((n-1)/2)*z: smear radius one,
background radius sqrt(n), pair distances at least sqrt(2).
"""
from fractions import Fraction as F

CLOSE=F(15,8)
QUERY=F(11,10)
HOLE=F(9,200)
INFLATED=1+HOLE
NEAR=QUERY+INFLATED
GREEDY=2*(QUERY+HOLE)
KAPPA=F(1,2**24)
PAIR_COST=F(1,8)*(1-CLOSE/2)**4
HOLE_COST=HOLE**4/4


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def packing_less_than(radius,integer):
    """Certify (1+sqrt(2)*radius)^2 < integer by exact squaring."""
    radius=F(radius)
    if radius<0 or type(integer)is not int or integer<1:
        raise ValueError('nonnegative radius and positive integer required')
    margin=integer-1-2*radius*radius
    return margin>0 and margin*margin>8*radius*radius


def validate():
    # Positive integral x^4(1-x)^4/(1+x^2) = 22/7-pi.
    quotient=[F(4),F(0),F(-4),F(0),F(5),F(-4),F(1)]
    product=[F(0)]*9
    for j,v in enumerate(quotient):product[j]+=v;product[j+2]+=v
    product[0]-=4
    require(product==[0,0,0,0,1,-4,6,-4,1],'Archimedes positive integral identity')
    require(sum((v/F(j+1) for j,v in enumerate(quotient)),F(0))==F(22,7),'pi upper coefficient')
    angle=F(44,49)
    cos_lower=1-angle**2/2+angle**4/24-angle**6/720
    require(0<angle<1 and cos_lower>F(31,50),'cos(2pi/7)>31/50')
    corners=[a*a+b*b-2*a*b*F(31,50) for a in (CLOSE,NEAR) for b in (CLOSE,NEAR)]
    require(max(corners)<CLOSE*CLOSE,'seven neighbors impossible')
    beta=(CLOSE*CLOSE+QUERY*QUERY-INFLATED*INFLATED)/(2*CLOSE*QUERY)
    require(QUERY>INFLATED and CLOSE*CLOSE>QUERY*QUERY-INFLATED*INFLATED,
            'own disk misses query circle and cap angle monotonicity')
    require(beta>0 and 4*beta*beta>3,'each neighbor covers less than pi/3')
    require(packing_less_than(NEAR,17),'endpoint neighborhoods contain <=16 points')
    require(packing_less_than(GREEDY,18),'greedy neighborhoods contain <=17 points')
    require(PAIR_COST==32*KAPPA,'matching pair cost pays 32 point units')
    require(HOLE_COST/F(17)>=KAPPA,'hole cost pays greedy population units')
    boundary_margin=F(11,4)-QUERY-HOLE
    require(boundary_margin>0 and boundary_margin**2>2,'boundary population <11sqrt(n)')
    require((F(2)-QUERY-HOLE)**2>F(1,2),'inner boundary packing annulus positive for n>=4')
    return {'close_pair_distance':CLOSE,'query_circle_radius':QUERY,
        'smoothed_hole_radius':HOLE,'near_radius':NEAR,'greedy_radius':GREEDY,
        'maximum_unmatched_near_neighbors':6,'endpoint_neighborhood_cap':16,
        'greedy_neighborhood_cap':17,'pair_energy_lower':PAIR_COST,
        'one_hole_energy_lower':HOLE_COST,'per_point_lower':KAPPA,
        'cap_cosine_lower':beta,'seventh_neighbor_cosine_comparison':cos_lower,
        'seven_neighbor_corner_distance_squared_upper':max(corners)}
