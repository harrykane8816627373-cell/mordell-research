"""Accept monotone-factor integral costs using only rational arithmetic."""
from fractions import Fraction as F
from pair_integral import (sqrt_bounds,acos_lower,radial_h_lower,
    pair_integral_lower,fixed_certificates,ROOT)
from reporting import require,write_report
import sys
sys.path.insert(0,str(ROOT/'joint-bounds'))
from subset_moment import log_interval
from overlap_cost import pair_cost_lower


def run():
    roots=0
    for den in (3,7,25,101):
        for num in range(51):
            x=F(num,den);lo,hi=sqrt_bounds(x)
            require(0<=lo<=hi and lo*lo<=x<=hi*hi,'outward root squares');roots+=1
    # Independent cosine series confirms the inverse-angle direction.
    angles=0
    for j in range(21):
        x=F(j,20);a=acos_lower(x)
        cosine_lower=sum(((-1)**k*a**(2*k)/F(__import__('math').factorial(2*k)) for k in range(14)),F(0))
        require(0<=a<2 and cosine_lower>=x,'atan angle lies below arccos');angles+=1
    for j in range(11,21):
        s=F(j,20);ll,lu=log_interval(s)
        require(radial_h_lower(s)<=-lu-F(1,2)+s*s/2,'log lower versus independent atanh enclosure')
    rows=fixed_certificates()
    for row in rows:
        c=row['cutoff_force']
        require(pair_integral_lower(c,64)<row['quadrature_lower'],'coarse integral comparison')
        require(row['accepted_cost_lower']>pair_cost_lower(c),'improvement over parabolic-only lower')
    invalid=0
    for args in [(F(7,5),),(2,),(F(8,5),0),(F(8,5),4096,24,11)]:
        try:pair_integral_lower(*args)
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid integral input accepted')
    return {'passed':True,'root_checks':roots,'independent_cosine_checks':angles,
            'fixed_costs':rows,'bins':4096,'log_terms':24,'atan_terms':12,'sqrt_mesh':10**10,
            'invalid_inputs_rejected':invalid,
            'scope':'Exact lower Riemann sums with separately monotone factors; both ordered pair terms included. Convexity and all-distance extension are proved analytically.'}


if __name__=='__main__':write_report('pair_integral_verification.json',run())
