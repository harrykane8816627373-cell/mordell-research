"""Exact tests for truncated radial energy transport and layer capacities."""
from fractions import Fraction as F
from pathlib import Path
import json

from radial_transport import (sqrt_lower,energy_mean_lower,affine_layers,pi_layers,
    cdf_mean_from_error,cdf_mean_lower,cdf_affine_layers,finite_family_best,
    root_upper,pi_interval,polar_gap_upper,smearing_penalty_upper)


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def independent_log(q,terms=64):
    q=F(q)
    if q<1:
        a,b=independent_log(1/q,terms)
        return -b,-a
    y=(q-1)/(q+1)
    lo=2*sum((y**(2*j+1)/(2*j+1) for j in range(terms)),F(0))
    return lo,lo+2*y**(2*terms+1)/((2*terms+1)*(1-y*y))


def run():
    roots=0
    for a in range(80):
        for b in range(1,8):
            q=F(a,b);lo=sqrt_lower(q,10000);hi=root_upper(q,2,10000)
            require(lo*lo<=q<=hi*hi,'outward root pair')
            require(hi-lo<=F(1,10000),'adjacent root grid endpoints')
            roots+=1
    # Flux and gradient-norm antiderivatives for w_a=min(1/a,1/R).
    norm_checks=0
    for a in (F(1,10),F(1,3),F(1),F(7,5)):
        for b in (2*a,3*a,10*a):
            # Outside a: integral_a^b R*(w')^2 dR=(a^-2-b^-2)/2.
            norm_partial=(1/a**2-1/b**2)/2
            norm_tail=1/(2*b*b)
            require(norm_partial+norm_tail==1/(2*a*a),'complete Dirichlet norm')
            # Delta w=R^-3 outside plus boundary jump -1/a^2.
            require((1/a-1/b)+1/b-1/a==0,'total radial Laplacian flux vanishes')
            require(1/b**3<=1/a**3,'positive Laplacian density upper bound')
            norm_checks+=1
    # At center zero the disk mean is explicit; the universal 1/8 upper
    # error is derived analytically, not fitted to these finite samples.
    center_checks=0
    for a in (F(1,5),F(1,2),F(1)):
        for r in (a/2,a,2*a,3*a):
            average=1/a if r<=a else (2*r-a)/(r*r)
            require(average<=1/a+r*r/(8*a**3),'centered disk mean comparison')
            require((F(1,1)/4)*(r*r/2)==r*r/8,'quadratic disk mean factor 1/8')
            center_checks+=1
    # Independent radial ring charges: D(ring_R-ring_S)=log(S/R).
    duality_checks=0
    for a in (F(1,5),F(1,2),F(1)):
        for R,S in ((F(1,4),F(1,2)),(F(1,3),F(1)),(F(1),F(2))):
            diff=min(1/a,1/R)-min(1/a,1/S)
            lo,hi=independent_log(S/R)
            require(diff*diff<=lo/(2*a*a),'radial ring Dirichlet duality')
            duality_checks+=1
    affine_capacity_checks=0
    pi_capacity_checks=0
    for t in range(1,6):
        # n=2*t^2 makes sqrt(8n)=4t rational.
        n=2*t*t;C=F(4*t)
        for m in range(1,21):
            R=F(m-1,C) if m>1 else F(0)
            for a in (F(1,10),F(1,3),F(1)):
                w=1/a if R==0 else min(1/a,1/R)
                require(m*w<=C+1/a,'affine capacity includes singleton')
                affine_capacity_checks+=1
        # Separation-only circle capacity; pi>3 gives a safe cutoff.
        n=t*t+2*t+1;root=F(t+1);a=1/(3*root)
        for m in range(1,21):
            R=F(m,3*root) if m>1 else F(0)
            w=1/a if R==0 else min(1/a,1/R)
            require(m*w<=pi_interval()[0]*root,'pi capacity with singleton cutoff')
            pi_capacity_checks+=1
    # Finite exact radial distributions with known Kolmogorov error.
    cdf_checks=0
    for m in range(1,16):
        qs=[F(2*j+1,m) for j in range(m)]
        epsilon=F(1,2*m)
        for a in (F(1,10),F(1,2),F(1),F(7,5)):
            actual_lower=sum((min(1/a,1/root_upper(q,2)) for q in qs),F(0))/m
            bound=cdf_mean_from_error(epsilon,a)
            require(bound<=actual_lower,'clipped CDF transport compared with exact atom data')
            cdf_checks+=1
    # Correct subtraction of orthogonal angular information.
    n=10000;gap=polar_gap_upper(n);a=F(1,10)
    ordinary=energy_mean_lower(n,gap,a)
    angular=energy_mean_lower(n,gap,a,nonradial_delta_lower=100)
    require(angular['mean_lower']>=ordinary['mean_lower'],'angular residual budget improves radial transport')
    require(ordinary['disk_mean_error_upper']==F(1,40),'smear error at n10000 cutoff1/10')
    require(ordinary['radial_energy_budget_upper']==gap,'default smear penalty zero')
    varied=energy_mean_lower(n,gap,a,h=F(19,20))
    require(varied['radial_energy_budget_upper']==gap+n*smearing_penalty_upper(F(19,20)),
            'variable smear energy budget includes penalty')
    require(varied['disk_mean_error_upper']==F(19,20)**2/(4*n*a**3),
            'variable smear center-average error uses h squared')
    examples=[]
    for n in (10000,1000000,100000000,1000000000000):
        gap=polar_gap_upper(n)
        a=root_upper(F(1,n),4)
        direct=affine_layers(n,gap,a)
        baseline=pi_layers(n,gap,a)
        epsilon=cdf_mean_lower(n,gap,a)['cdf_error_upper']
        cdf_a=min(F(7,5),root_upper(2*epsilon,2))
        cdf=cdf_affine_layers(n,gap,cdf_a)
        family=finite_family_best(n,gap)
        require(family['best']['real_layer_lower']>=direct['real_layer_lower'],
                'prescribed finite family includes basic cutoff')
        require(1<=direct['minimum_layers']<=n,'integer necessary count range')
        examples.append({'n':n,'proved_gap_upper':gap,'a':a,
            'direct_affine_minimum_layers':direct['minimum_layers'],
            'direct_real_layer_lower':direct['real_layer_lower'],
            'pi_baseline_minimum_layers':baseline['minimum_layers'],
            'CDF_affine_minimum_layers':cdf['minimum_layers'],
            'CDF_a':cdf_a,'CDF_error_upper':epsilon,
            'finite_family_best_minimum_layers':family['best']['minimum_layers'],
            'finite_family_best_a':family['best']['a'],'finite_family_best_h':family['best']['h']})
    invalid=[lambda:energy_mean_lower(1,1,1),lambda:energy_mean_lower(8,-1,1),
        lambda:energy_mean_lower(8,1,0),lambda:energy_mean_lower(8,1,2),
        lambda:energy_mean_lower(8,1,1,nonradial_delta_lower=2),
        lambda:pi_layers(8,1,F(1,100000)),lambda:cdf_mean_from_error(2,1),
        lambda:finite_family_best(8,1,multipliers=[])]
    for call in invalid:
        try:
            call()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid input accepted')
    return {'passed':True,'arithmetic':'integer/Fraction only; outward roots and logs',
        'outward_root_checks':roots,'Dirichlet_norm_and_flux_checks':norm_checks,
        'centered_disk_mean_checks':center_checks,'radial_duality_checks':duality_checks,
        'affine_layer_capacity_checks':affine_capacity_checks,
        'pi_layer_capacity_checks':pi_capacity_checks,'clipped_CDF_checks':cdf_checks,
        'rejected_invalid_cases':len(invalid),'effective_examples':examples,
        'scope':'Finite arithmetic checks. Distributional Laplacian, disk mean comparison and all-n transport are proved analytically in RADIAL_TRANSPORT.md. Capacity inputs are separately proved; no optimum configuration is classified.'}


if __name__=='__main__':
    report=run()
    Path(__file__).with_name('radial_transport_verification.json').write_text(
        json.dumps(report,default=str,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,indent=2))
