"""Check the rational geometry and remainder constants for central projection."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_bytes())


def main():
    manifest = read(HERE/'source_manifest.json')
    for folder, entries in ((HERE, manifest['files']), (ROOT, manifest['dependencies'])):
        for name, digest in entries.items():
            require(sha256((folder/name).read_bytes()).hexdigest()==digest, 'fixed input '+name)
    def proof(folder):
        result = read(ROOT/folder/'verification_report.json')
        require(result['passed'] and result['applies_to_every_global_maximizer'], folder+' geometry')
        return result
    bands = proof('ten-point-angular-radii')
    lo, hi = map(F, bands['outer_squared_radius_open_interval'])
    require((bands['normalization_mass'], bands['exact_inner_count'], bands['exact_outer_count'])==(9,2,8),
            'pair populations')
    pairs = proof('ten-point-tight-pairing')['pairing']
    cos_pair = -F(pairs['cyclic_rank_pairing_cosine_strict_upper'])
    abs_inner_m = F(87,1000)
    middle = proof('ten-point-inner-midpoint')
    require(F(middle['inner_midpoint_norm_strict_upper'])==abs_inner_m, 'inner midpoint bound')
    abs_outer_m = F(43,250)
    require(F(pairs['paired_outer_midpoint_squared_strict_upper']) < abs_outer_m**2, 'outer midpoint bound')
    inner_lo, inner_d = map(F, proof('ten-point-inner-edge')['inner_squared_distance_open_interval'])
    _, gap_cos_hi = map(F, proof('ten-point-outer-gap')['outer_successive_angle_cosine_open_interval'])
    r_min, r_max = F(943,1000), F(1128,1000)
    a_min, inner_a_max = F(9359,10000), F(3641,10000)
    d_min, turn_shift = F(4,7), F(3,250)
    require(r_min**2<lo and r_max**2>hi and a_min**2<lo*(1+cos_pair)/2,
            'outer half-difference projection')
    require(inner_a_max**2>inner_d/4 and a_min-inner_a_max>d_min, 'inner-outer separation')
    tan_sq = (1-cos_pair)/(1+cos_pair)
    require(((r_max-r_min)/(r_max+r_min))**2*tan_sq<turn_shift**2, 'angular shift of projected pair')
    require(gap_cos_hi+2*turn_shift<1 and
            2*a_min**2*(1-gap_cos_hi-2*turn_shift)>d_min**2, 'outer-outer separation')
    eta_max = F(43,125)
    require(abs_inner_m<abs_outer_m and 2*abs_outer_m==eta_max, 'difference of two midpoints')
    rho, local_m, long = F(3,8), F(1,30), F(99,100)
    require((eta_max/d_min)**2<rho, 'uniform small perturbation ratio')
    require(a_min**2+inner_lo/4>long**2 and a_min**2>inner_lo/4, 'one long edge per pair')
    def factor(eta_squared):
        return 7*sum((1+eta_squared/length**2)/(length**4*(1-eta_squared/length**2))
                     for length in (d_min,long))
    require(factor(eta_max**2)<150 and factor(2*local_m**2)<74, 'quartic coefficients')
    def sixth_order(eta_squared):
        return 28*sum(1/(length**2*(length**2-eta_squared)**2) for length in (d_min,long))
    require(sixth_order(eta_max**2)<2017 and sixth_order(2*local_m**2)<846,
            'sixth-order midpoint stationarity')
    require(4*r_max/(d_min*(d_min**2-eta_max**2))<38, 'force coefficient')
    print(json.dumps({'passed':True, 'normalization_mass':9, 'applies_to_every_global_maximizer':True,
                      'projected_cross_pair_distance_strict_lower':str(d_min),
                      'one_cross_pair_distance_strict_lower':str(long),
                      'uniform_quartic_energy_coefficient':150,
                      'smaller_midpoint_norm_ceiling':str(local_m),
                      'small_midpoint_quartic_energy_coefficient':74,
                      'midpoint_quadratic_sixth_order_coefficient':2017,
                      'small_midpoint_quadratic_sixth_order_coefficient':846,
                      'quadratic_force_coefficient':38}, indent=2))


if __name__ == '__main__':
    main()
