"""Rotation selection and cross terms by independent Gaussian-Laplacian integration."""
from fractions import Fraction as F
from gaussian_energy import (standard_probe,center_power,exp_negative,bessel_sum,
    add,sub,scale,mul)
from legacy_projection import cross_mode,dual_lower
from acceptance import require,write_report


def independent_bubble_cross(probe,k,S):
    t=probe['t'];rho=probe['center_squared_radius'];x=rho/(4*t);y=S/(4*t);z=S*rho/(16*t*t)
    def L(order,j):
        v=mul(exp_negative(x),bessel_sum(z,order,y,j))
        return scale(v,S**(order+j+1)/(2*(4*t)**order))
    P=sub(scale(L(k,0),S),L(k,1))
    distance=add(sub(scale(L(k,1),S),L(k,2)),scale(P,rho))
    distance=sub(distance,scale(sub(scale(L(k+1,0),S),L(k+1,1)),rho))
    distance=sub(distance,sub(scale(L(k-1,1),S),L(k-1,2)))
    coefficient=sub(scale(P,1/t),scale(distance,1/(4*t*t)))
    return [scale(coefficient,p) for p in center_power(probe,k)]


def run():
    probes=[standard_probe(),standard_probe(F(1,24),F(3,4))];selection=0;crosses=0
    for probe in probes:
        require(center_power(probe,0)==(0,0),'equal-radius radial cancellation')
        for k in range(1,16):
            if k%4:
                require(center_power(probe,k)==(0,0),'quarter-turn angular selection');selection+=1
        for S in (F(5,4),F(3,2)):
            direct=cross_mode(probe,4,S)[2:];other=independent_bubble_cross(probe,4,S)
            for a,b in zip(direct,other):
                require(max(a[0],b[0])<=min(a[1],b[1]),'gradient/boundary versus Gaussian-Laplacian cross')
                require(max(a[1]-a[0],b[1]-b[0])<F(1,10**25),'cross enclosure precision');crosses+=1
        require(all(v==(0,0) for k in (2,3) for v in cross_mode(probe,k,F(5,4))),'exact compatibility with low angular modes')
    # Synthetic Hilbert-space counterexample to independently adding overlapping probes.
    H=[[(F(1),F(1)),(F(1),F(1))],[(F(1),F(1)),(F(2),F(2))]]
    m=[(F(1),F(1))]*2
    require(dual_lower(H,m,[1,0])==1<F(3,2),'full Gram avoids naive overlapping sum')
    return {'passed':True,'rotation_selection_checks':selection,'independent_cross_checks':crosses,
            'synthetic_true_energy':'1','synthetic_naive_separate_sum':'3/2',
            'scope':'Gaussian probes are orthogonal to radial and k=2,3 sectors. Existing k=4 or other multiples of four must share the full Gram; they cannot be independently added.'}


if __name__=='__main__':write_report('compatibility_verification.json',run())
