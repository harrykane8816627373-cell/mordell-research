"""Layer circumference from capacity plus one radial Dirichlet test.

Inputs E are independently proved global gap upper bounds.  The interface
does not infer global optimality from a radius list or from small residuals.
"""
from fractions import Fraction as F
from math import isqrt
from pair_integral import sqrt_bounds,certified_pi_upper,FIXED_COSTS


def validate_n(n):
    if type(n)is not int or n<2:raise ValueError('integer n>=2 required')


def radial_test(radius,a):
    r,a=F(radius),F(a)
    if r<0 or not 0<a<=1:raise ValueError('nonnegative radius and 0<a<=1')
    return (1-r)/max(a,r)


def test_background(a):
    a=F(a)
    if not 0<a<=1:raise ValueError('0<a<=1 required')
    return 1-a+a*a/3


def test_dirichlet(a):
    a=F(a)
    if not 0<a<=1:raise ValueError('0<a<=1 required')
    return (1+1/(a*a))/2


def radial_average_lower(n,proven_gap_upper,a,mesh=10**12):
    validate_n(n);E,a=F(proven_gap_upper),F(a)
    if E<0:raise ValueError('nonnegative proved gap required')
    H=test_dirichlet(a)
    error=sqrt_bounds(E*H,mesh)[1]/n
    lower=test_background(a)-error-F(1,8*n)/a**3
    return {'a':a,'background_average':test_background(a),'dirichlet_norm_squared':H,
            'energy_error_upper':error,'smearing_error_upper':F(1,8*n)/a**3,
            'point_average_lower':lower}


def layer_exclusion(n,maximum_layers,proven_gap_upper,a,cutoff=F(3151,2000)):
    validate_n(n)
    K=maximum_layers;E=F(proven_gap_upper);a=F(a);c=F(cutoff)
    if type(K)is not int or K<1 or 2*K>=n:raise ValueError('integer 1<=K<n/2 required')
    if c not in FIXED_COSTS:raise ValueError('a separately fixed integral cutoff is required')
    test=radial_average_lower(n,E,a)
    root_n_lower=sqrt_bounds(n,10**12)[0]
    radius_sum_over_n=F(K)/root_n_lower+F(K,4*n)/a-test['point_average_lower']/4
    edges=n-2*K
    distance_upper=2*certified_pi_upper()*n*max(F(0),radius_sum_over_n)/edges
    cost=FIXED_COSTS[c]
    available=distance_upper<=c
    excluded=radius_sum_over_n<0 or (available and edges*cost>E)
    return {'n':n,'maximum_layers':K,'proven_gap_upper':E,'radial_test':test,
            'sum_force_radii_over_n_upper':radius_sum_over_n,
            'adjacent_edge_count_lower':edges,'average_adjacent_distance_upper':distance_upper,
            'cutoff_force':c,'accepted_cost_per_pair':cost,
            'cutoff_reached':available,'pair_energy_lower':edges*cost if available else F(0),
            'entire_at_most_K_class_excluded':excluded,
            'scope':'Necessary for global optimizers using same-layer capacity, force radii<=sqrt(n), and the stated proved energy input; not a candidate optimizer test.'}


def elementary_gap_upper(n):
    validate_n(n)
    if n<6:raise ValueError('inherited elementary gap bound requires n>=6')
    root=isqrt(n);root_upper=root if root*root==n else root+1
    return F(n,100)+F(77,6)*root_upper+F(35,4)
