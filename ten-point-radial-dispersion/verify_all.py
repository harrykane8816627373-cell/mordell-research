"""Rational checks for the exclusion of almost equal outer radii."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
EPSILON=F(1,25)


def need(ok, description):
    if not ok:
        raise ArithmeticError(description)


def arithmetic():
    u_lo,u_hi=F(41,400),F(11,80)
    e,d,s=EPSILON,EPSILON/2,5*EPSILON/3
    v_hi,radius=F(11,10),F(53,50)
    coarse_lo=(9-2*u_hi-F(1,20))/8
    need((9-2*u_lo)/8<v_hi and radius**2>v_hi+d,'outer mean and radii')
    need(coarse_lo-3*u_hi-F(1,40)>3*radius/5,'center and first Ward bound')
    v_lo=(9-2*u_hi-s*s/2)/8
    need(v_lo>F(109,100) and v_lo-d>0,'improved outer mean range')
    T=u_hi+s*s/4
    p_upper=radius**2+s*radius+T
    coefficient={8:F(1),7:s}
    remainders={}
    for j in range(7,2,-1):
        E=10*e*p_upper*sum(coefficient[k]*radius**(k-j-1) for k in range(j+1,9))
        remainders[j]=E
        denominator=(9-j)*(j-1)
        coefficient[j-1]=s*coefficient[j]+F((j+1)*(11-j),denominator)*T*coefficient[j+1]+E/denominator
    lower6=8*(10*v_lo-7)/(2*(17-10*v_lo))*(u_lo-s*s/4)-s*s-remainders[7]/12
    need(u_lo-s*s/4>0 and lower6>0,'strict q6 lower')
    linear_error=e*(s*radius+radius**2)/v_lo**3
    quadratic_error=F(15,4)*radius**2*e**2/(v_lo-d)**4
    lower_ratio=lower6/v_hi**2-linear_error-quadratic_error
    need(lower_ratio>0,'positive lower for the approximate reciprocal coefficient')
    reciprocal=(v_lo-d)**4*lower_ratio
    need(coefficient[2]<F(171,1000)<F(173,1000)<reciprocal,'strict reciprocal-coefficient contradiction')
    return {'passed':True,'normalization_mass':9,'applies_to_every_global_maximizer':True,
            'sum_absolute_outer_squared_radius_deviations_strict_lower':str(e),
            'max_outer_squared_radius_deviation_from_mean_strict_lower':str(e/8),
            'sum_squared_outer_squared_radius_deviations_strict_lower':str(e*e/8),
            'outer_squared_radius_span_strict_lower':str(e/4),
            'reciprocity_contradiction_rational_margin':str(reciprocal-coefficient[2]),
            'reciprocity_contradiction_margin_strict_lower':'1/500'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true'); args=parser.parse_args()
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for folder,key in ((HERE,'files'),(ROOT,'dependencies')):
        for name,digest in manifest[key].items():
            need(sha256((folder/name).read_bytes()).hexdigest()==digest,'fixed '+name)
    radial=json.loads((ROOT/'ten-point-refined-radii/verification_report.json').read_bytes())
    pair=json.loads((ROOT/'ten-point-conditioned-inner/verification_report.json').read_bytes())
    need(radial['passed'] and radial['applies_to_every_global_maximizer'] and radial['normalization_mass']==9
         and radial['exact_inner_count']==2 and radial['exact_outer_count']==8
         and 4*F(radial['midpoint']['inner_midpoint_squared_strict_upper'])<F(1,10),'fixed radial and midpoint inputs')
    need(pair['passed'] and pair['applies_to_every_global_maximizer'] and pair['normalization_mass']==9
         and list(map(F,pair['inner_squared_distance_open_interval']))==[F(41,100),F(11,20)],'fixed inner distance input')
    result=json.dumps(arithmetic(),indent=2)+'\n'
    if args.record_report:(HERE/'verification_report.json').write_bytes(result.encode())
    print(result,end='')


if __name__=='__main__':main()
