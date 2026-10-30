"""Global support lines: derivative brackets and an independent tail integral."""
from fractions import Fraction as F
from affine_pair import (tail_interval,acos_upper,acos_lower,sqrt_bounds,
    pi_bounds,fixed_line,support_line,pair_integral_upper,pair_integral_lower)
from acceptance import require,write_report


def run():
    rows=[fixed_line(),support_line(F(15753,10000),F(7877,5000),F(21,250))]
    # Independent radial density integration of the two-disk overlap area.
    density_checks=[]
    pl,pu=pi_bounds()
    for d in (F(3,2),F(9,5)):
        step=(2-d)/256;lo=hi=F(0)
        for j in range(256):
            l=d+j*step;u=l+step
            blo,bhi=sqrt_bounds(1-u*u/4,10**12)
            area_lo=max(F(0),2*acos_lower(u/2,12,10**12)-u*bhi)
            blo,bhi=sqrt_bounds(1-l*l/4,10**12)
            area_hi=2*acos_upper(l/2)-l*blo
            lo+=(u*u-l*l)*area_lo/pu;hi+=(u*u-l*l)*area_hi/pl
        exactlo,exacthi=tail_interval(d)
        require(lo<exactlo<=exacthi<hi,'closed tail inside independent density Riemann bounds')
        density_checks.append({'d':d,'density_integral_lower':lo,'tail_interval':[exactlo,exacthi],'density_integral_upper':hi})
    samples=0
    for line in rows:
        for d in (F(3,2),F(8,5),F(7,4),F(19,10),F(199,100)):
            require(pair_integral_lower(d,512)>=line['alpha']-line['beta']*d,'independent off-contact support checks');samples+=1
    u=pair_integral_upper(F(8999,5000))
    require(pair_integral_lower(F(8999,5000))<u<F(71,100000),'upper integral used for separate-relaxation comparison')
    invalid=0
    for args in [(F(17901,10000),F(8951,5000),F(1,10)),(F(7,5),F(3,2),F(1,10)),(F(3,2),F(3,2),F(1,10))]:
        try:support_line(*args)
        except (ValueError,ArithmeticError):invalid+=1
        else:raise ArithmeticError('invalid support bracket accepted')
    return {'passed':True,'global_support_lines':rows,'independent_tail_density_checks':density_checks,
            'off_contact_checks':samples,'comparison_pair_upper':u,'invalid_inputs_rejected':invalid,
            'scope':'Global validity follows from convexity and the derivative bracket, not from finite distance sampling.'}


if __name__=='__main__':write_report('affine_pair_verification.json',run())
