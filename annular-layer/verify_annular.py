"""Restricted affine fans certify population-sensitive layer duals."""
from fractions import Fraction as F
from acceptance import ROOT,require,write_report
from annular_dual import (fan,floor_mesh,pi_bounds,annular_certificate,
                         asymptotic_lower,finite_lower,restricted_line)


def configured(h,price,last,M=2048):
    first=floor_mesh(pi_bounds()[0]/(2*h),100000)
    centers=[first]+[F(i,100) for i in range((first*100).__ceil__(),int(last*100)+1,2)]
    packet=fan(h,centers)
    cert=annular_certificate(packet,price,[F(j*j,M*M) for j in range(M+1)])
    return packet,cert


def run():
    reports=[]
    for h,price,e,last,claim in [(F(109,100),F(12,25),F(1,100),F(33,20),F(513,1000)),
                                (F(51,50),F(69,1000),F(1,800),F(183,100),F(139,250))]:
        packet,cert=configured(h,price,last)
        bound=asymptotic_lower(cert,e)
        require(bound>claim,'strict annular asymptotic layer improvement')
        # Test the two exact construction rules, not a float approximation to P.
        for line in packet['lines']:
            require(restricted_line(packet,line['beta'])==line['alpha'],'fan reproduces every certified knot')
        left,right=packet['lines'][2:4];mid=(left['beta']+right['beta'])/2
        require(restricted_line(packet,mid)==(left['alpha']+right['alpha'])/2,'positive affine mixing')
        top=packet['lines'][-1];beta=2*top['beta']
        require(restricted_line(packet,beta)==top['alpha']+top['beta']*packet['domain_lower'],'restricted-domain tilt')
        require(packet['domain_lower']<=pi_bounds()[0]/(2*h),'actual capacity enforces the full restricted domain')
        # Independent telescoping expansion of the background radial integral.
        rows=cert['bands'];total=rows[-1]['alpha']
        total+=sum(((rows[j]['alpha']-rows[j+1]['alpha'])*rows[j]['b']**2 for j in range(len(rows)-1)),F(0))
        require(total==cert['uniform_radial_integral_lower'],'independent radial summation by parts')
        require(all(a['alpha']>=b['alpha'] for a,b in zip(rows,rows[1:])),'decreasing radial test')
        reports.append({'smearing_multiplier':h,'layer_price':price,'proven_gap_per_n_upper':e,
            'restricted_domain_lower':packet['domain_lower'],'certified_global_lines':packet['lines'],
            'partition':'r_j=(j/2048)^2, j=0,...,2048','number_of_bands':len(rows),
            'uniform_radial_integral_lower':cert['uniform_radial_integral_lower'],
            'asymptotic_layer_lower':bound,'strict_rational_claim':claim,
            'strict_margin_over_claim':bound-claim,
            'input_scope':'elementary effective budget only' if e==F(1,100) else 'retains the previously cited external linear-energy and triangular-trial asymptotic input; no finite starting n'})
    small=annular_certificate(packet,F(69,1000),[0,F(1,3),F(2,3),1])
    v=finite_lower(small,9,3,[1,4,4],layer_counts=[1,1,1])
    shifted=finite_lower(small,9,3,[1,4,4],dirichlet_lower=F(1,20),layer_counts=[1,1,1])
    require(shifted-v==F(1,20),'one independent same-h Dirichlet energy is retained exactly once')
    require(finite_lower(small,9,3,[1,4,4])<=v,'unknown band-layer counts use a safe maximum correction')
    rejected=0
    for action in [lambda:fan(0,[F(3,2)]),lambda:annular_certificate(packet,0,[0,1]),
                   lambda:annular_certificate(packet,1,[0,F(3,4),F(1,2),1]),
                   lambda:finite_lower(small,9,3,[1,4,3]),
                   lambda:finite_lower(small,9,3,[1,4,4],layer_counts=[0,1,2]),
                   lambda:finite_lower(small,9,3,[1,4,4],dirichlet_lower=-1)]:
        try:action()
        except ValueError:rejected+=1
    require(rejected==6,'invalid annular inputs rejected')
    return {'passed':True,'asymptotic_certificates':reports,
        'small_finite_population_interface_lower':v,'small_finite_same_h_energy_shift':F(1,20),
        'invalid_inputs_rejected':rejected,
        'scope':'Finite annular population/layer-count inequality plus fixed 2048-band exact duals. Radial convergence supplies the asymptotic quantifier; the external-input result is explicitly separated. No microscopic lattice or optimality classification.'}


if __name__=='__main__':write_report('annular_verification.json',run())
