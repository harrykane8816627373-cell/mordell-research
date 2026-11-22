"""Exact radial dispersion inequality for the half-turn branch."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EPSILON = F(47,1000)


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_bytes())


def verify():
    manifest=read(HERE/'source_manifest.json')
    for base,key in ((HERE,'files'),(ROOT,'dependencies')):
        for name,value in manifest[key].items():
            require(sha256((base/name).read_bytes()).hexdigest()==value, 'fixed input '+name)
    inner=read(ROOT/'ten-point-inner-edge/verification_report.json')
    require(inner['passed'] and inner['normalization_mass']==9 and inner['unresolved']==0
            and inner['applies_to_every_global_maximizer'], 'accepted inner edge interval')
    a_lo,a_hi=(F(d)/4 for d in inner['inner_squared_distance_open_interval'])
    require((a_lo,a_hi)==(F(43,400),F(53,400)), 'inner square-plane point range')
    v_lo,v_hi=(F(9,2)-a_hi)/4,(F(9,2)-a_lo)/4
    radius=v_hi+EPSILON/2
    b_lo=(17-10*a_lo)*a_lo
    b_hi=(17-10*a_hi)*a_hi
    require(0<a_lo<a_hi<F(17,20) and 0<v_lo-EPSILON/2<v_hi<F(17,15),
            'positive mass and Ward bounds throughout')

    # q_4 is one. The first Ward equation bounds q_3, which seeds the outer-force recurrence.
    upper={4:F(1),3:(b_hi+10*EPSILON*radius)/(17-10*v_hi)}
    for j in (2,1):
        lower_denominator=(4-j)*(j+F(11,2)-5*v_hi)
        max_numerator=a_hi*(j+1)*(5*v_hi-j-F(1,2))
        require(lower_denominator>0 and 5*v_lo-j-F(1,2)>0,'definite recurrence coefficients')
        residual=5*EPSILON*(radius+a_hi)*sum(upper[k]*radius**(k-j-1) for k in range(j+1,5))
        upper[j]=(max_numerator*upper[j+1]+residual)/lower_denominator

    # The Ward equation evaluates the linear reciprocal correction exactly.
    lower_radius=v_lo-EPSILON/2
    numerator=b_lo/(5*v_hi**3)-(17-15*v_lo)*upper[3]/(5*v_lo**3)
    error=3*radius*EPSILON**2/(2*lower_radius**4)
    require(numerator>error, 'positive lower bound before multiplication by the root product')
    lower_q1=lower_radius**4*(numerator-error)
    gap=lower_q1-upper[1]
    require(gap>F(1,1000), 'contradiction for almost equal outer radii')
    return {'passed':True,'normalization_mass':9,'assumption':'half-turn symmetry of a global maximizer',
            'outer_square_plane_sum_absolute_deviations_strict_lower':str(EPSILON),
            'full_outer_sum_absolute_squared_radius_deviations_strict_lower':str(2*EPSILON),
            'outer_squared_radius_span_strict_lower':str(EPSILON/2),
            'coefficient_contradiction_margin_strict_lower':'1/1000'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    result=verify()
    message=json.dumps(result,indent=2)+'\n'
    path=HERE/'verification_report.json'
    if args.record_report:
        path.write_text(message,encoding='ascii',newline='\n')
    else:
        require(result==read(path),'frozen report')
    print(message,end='')


if __name__=='__main__':
    main()
