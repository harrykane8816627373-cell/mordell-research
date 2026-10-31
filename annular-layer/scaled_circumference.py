"""Variable-smearing layer circumference; hypotheses are proved separately."""
from fractions import Fraction as F
import sys
from acceptance import ROOT,require
sys.path.insert(0,str(ROOT/'radial-angular'))
from radial_layers import smearing_penalty_upper
sys.path.insert(0,str(ROOT/'joint-variational'))
from affine_pair import sqrt_bounds,pi_bounds


def checked(n,K,E,a,h):
    E,a,h=map(F,(E,a,h))
    if type(n)is not int or n<2 or type(K)is not int or not 1<=K<n:
        raise ValueError('integers n>=2 and 1<=K<n')
    if E<0 or not 0<a<=1 or h<=0:
        raise ValueError('nonnegative proven gap, 0<a<=1, h>0')
    return E,a,h


def certificate(n,K,E,a,h,cutoff,pair_cost_lower):
    """Conditional exclusion from an independently certified P(cutoff) lower."""
    E,a,h=checked(n,K,E,a,h)
    cutoff,cost=F(cutoff),F(pair_cost_lower)
    if cutoff<=0 or cost<=0:raise ValueError('positive cutoff and certified pair cost')
    budget=E+n*smearing_penalty_upper(h)
    H=(1+1/(a*a))/2
    transport=1-a+a*a/3-sqrt_bounds(budget*H)[1]/n-h*h/(8*n*a**3)
    circumference=K*sqrt_bounds(n)[1]+F(K,4)/a-F(n,4)*transport
    # m-1 artificial edge mass is valid for EVERY layer, including m=1,2.
    L=n-K
    distance=2*pi_bounds()[1]*max(F(0),circumference)/(h*L)
    pair_lower=L*cost
    return {'n':n,'maximum_layers':K,'smearing_multiplier':h,
        'smearing_energy_budget_upper':budget,'radial_test_average_lower':transport,
        'distinct_radius_sum_upper':circumference,'edge_mass_lower':L,
        'scaled_mean_distance_upper':distance,'certified_cutoff':cutoff,
        'pair_cost_lower':cost,'pair_sum_lower_if_cutoff_applies':pair_lower,
        'whole_class_excluded':circumference<0 or (distance<=cutoff and pair_lower>budget),
        'scope':'Global radius and layer capacity, normalized actual measurements, independently proved gap upper and pair cost are required. The m-1 mass is a convex lower bound on disjoint layer pairs, not a claim of n-K distinct actual edges.'}
