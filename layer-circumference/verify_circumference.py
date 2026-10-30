"""Radial-test integrals and layer-capacity ledgers."""
from fractions import Fraction as F
from circumference import (radial_test,test_background,test_dirichlet,
    radial_average_lower,layer_exclusion)
from reporting import require,write_report


def run():
    integrals=[]
    for a in (F(1),F(1,2),F(1,10),F(1,10000)):
        # Directly split the radial density integral 2r dr at r=a.
        inner=(a*a-F(2,3)*a**3)/a
        outer=1-2*a+a*a
        norm_inner=F(1,2)
        norm_outer=1/(2*a*a)
        jump=-1/(a*a)+1/a
        require(inner+outer==test_background(a),'independent background integral')
        require(norm_inner+norm_outer==test_dirichlet(a) and jump<=0,'norm and distributional jump')
        require(radial_test(0,a)==1/a and radial_test(1,a)==0 and radial_test(2,a)==-F(1,2),'origin/boundary/infinity branches')
        integrals.append({'a':a,'background':inner+outer,'norm_squared':norm_inner+norm_outer,'derivative_jump':jump})
    ledgers=0
    rootn=100
    for a in (F(1,2),F(1,10),F(1,100)):
        for m in range(1,51):
            for increment in (F(0),F(1,100),F(1,3)):
                q=F(m-1,4*rootn)+increment
                if q>1:continue
                require(m*radial_test(q,a)<=4*rootn*(1-q)+1/a,'layer capacity applied to radial test');ledgers+=1
    x=radial_average_lower(10**16,F(10**14),F(1,10000))
    require(x['point_average_lower']>F(9998,10000),'small finite radial-test error')
    # Every rejection is an input or range check, never a conclusion about an optimizer.
    invalid=0
    cases=[lambda:radial_test(-1,F(1,2)),lambda:test_background(0),lambda:test_dirichlet(2),
           lambda:radial_average_lower(1,1,F(1,2)),lambda:radial_average_lower(8,-1,F(1,2)),
           lambda:layer_exclusion(8,4,1,F(1,2)),lambda:layer_exclusion(8,1,1,F(1,2),F(3,2))]
    for fn in cases:
        try:fn()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid transport input accepted')
    return {'passed':True,'radial_integrals':integrals,'capacity_ledgers':ledgers,
            'finite_radial_test':x,'invalid_inputs_rejected':invalid,
            'scope':'Ledgers do not assert stationarity or normalization of their synthetic single-layer inputs. Universal radial transport is established by the weak Laplacian proof.'}


if __name__=='__main__':write_report('circumference_verification.json',run())
