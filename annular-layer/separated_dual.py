"""Global separation leaves at most one unpaid small singleton layer."""
from fractions import Fraction as F
import sys
from acceptance import ROOT
sys.path.insert(0,str(ROOT/'joint-variational'))
from affine_pair import pi_bounds,sqrt_bounds
from shared_budget import exact_certificate


def strengthened_certificate(packet,K,line,force_radius_lower=None):
    """h=1 only. The global force-scale separation sqrt(2) is a hypothesis.

    A supplied force_radius_lower must hold for ALL actual point centers;
    it is not inferred from a query radius or a small numerical residual.
    """
    alpha,beta=F(line['alpha']),F(line['beta'])
    pl=pi_bounds()[0];s2=sqrt_bounds(2)[0]
    if alpha<0 or beta<=0 or alpha>s2*(pl-1)*beta:
        raise ValueError('certified line with alpha/beta <= sqrt(2)(pi-1)')
    exceptional=1
    if force_radius_lower is not None:
        radius=F(force_radius_lower)
        if radius<0 or 2*pl*beta*radius<alpha:
            raise ValueError('proved minimum center radius must exclude every unpaid singleton')
        exceptional=0
    old=exact_certificate(packet,K,line)
    added_mass=2*K-exceptional
    return {'shared_gap_lower':old['shared_gap_lower']+alpha*added_mass,
        'guaranteed_paid_mass':packet['n']-exceptional,
        'exceptional_singleton_cap':exceptional,'extra_paid_mass_over_prior':added_mass,
        'improvement_over_prior':alpha*added_mass,
        'force_radius_lower_hypothesis':force_radius_lower,
        'scope':'Default h=1; actual global separation sqrt(2), radius/layer capacity, proved full radial Gram and certified pair line are required. A supplied all-center radius lower bound is an additional hypothesis.'}
