"""Certify the joint angular and inner-separation Hessian and its two phase tubes."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from phases import make_phases
from tubes import EPS, S, Tube

HERE = Path(__file__).resolve().parent
LOW = F(107, 1000)
HIGH = F(133, 1000)


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_bytes())


def prerequisites(repo):
    manifest = read(HERE / 'source_manifest.json')
    for base, key in ((HERE, 'files'), (repo, 'dependencies')):
        for file, digest in manifest[key].items():
            need(sha256((base/file).read_bytes()).hexdigest() == digest, 'fixed input '+file)

    def proof(name):
        data = read(repo / ('ten-point-'+name) / 'verification_report.json')
        need(data['passed'] and data['applies_to_every_global_maximizer']
             and data['normalization_mass'] == 9 and data['unresolved'] == 0,
             'accepted global prerequisite '+name)
        return data

    radii = proof('angular-radii')
    outer = tuple(map(F, radii['outer_squared_radius_open_interval']))
    need((radii['exact_inner_count'], radii['exact_outer_count']) == (2, 8)
         and outer == (F(89,100), F(1271,1000)), 'two inner points and eight outer points')
    edge = proof('inner-edge')
    w0, w1 = tuple(F(t)/4 for t in edge['inner_squared_distance_open_interval'])
    midpoint = F(proof('inner-midpoint')['inner_midpoint_norm_strict_upper'])
    parallel = F(edge['inner_midpoint_axis_projection_strict_upper'])
    need((w0,w1,midpoint,parallel) == (F(43,400), F(53,400), F(87,1000), F(7,125)),
         'inner pair geometry')
    factors = S/(9-2*w0), S/(9-2*w1)
    need(LOW < w0*factors[0] < w1*factors[1] < HIGH
         and midpoint**2*factors[1] < F(88,1000)**2
         and parallel**2*factors[1] < F(57,1000)**2
         and F(940,1000)**2 < outer[0]*factors[0]
         and outer[1]*factors[1] < F(1130,1000)**2, 'domain in fixed outer scale')
    need(S-2*F(88,1000)**2 > 8*F(109,100), 'positive sum of the radial affine correction')

    low, high = map(F, proof('outer-gap')['outer_successive_angle_cosine_open_interval'])
    need((low,high) == (F(59,100), F(39,50)), 'successive angular gap interval')
    q = F(127,125)
    b = 1/(q+1)
    need((outer[1]+outer[0])**2 < 4*q*q*outer[1]*outer[0]
         and 2-q*q-q*high > 0 and high*q+high**2 < 2, 'curvature monotonicity on the outer shell')

    def curvature(c):
        return (1-q*c)/(q-c)**2

    need(low > F(1,2) and curvature(low) > F(11,5) and curvature(2*low*low-1) > F(3,4),
         'one-step and two-step curvature')
    need(F(283,200)**2 > 2 and 8*b+min((F(11,5)-b)*F(117,200)+2*(F(3,4)-b),
            2*(F(11,5)-b)) > F(547,100), 'outer negative Hessian on the mean-zero space')
    basin = read(repo / 'ten-point-quartic-midpoints/verification_report.json')
    need(basin['passed'] and basin['normalization_mass']==9 and F(basin['half_difference_radius'])==F(7,250),
         'mass-nine global half-difference basin')


def hessian(compiler):
    phase_input = make_phases(512)
    with TemporaryDirectory(prefix='mordell-joint-angular-inner-') as scratch:
        binary = Path(scratch) / 'bounds.exe'
        subprocess.run([compiler, '-std=c++17', '-O3', str(HERE/'bounds.cpp'), '-o', str(binary)],
                       check=True, timeout=120)
        result = subprocess.run([str(binary)], input=phase_input, text=True, check=True,
                                capture_output=True, timeout=600)
    cells, total, denominator = map(int, result.stdout.split())
    need(cells == 23986176 and denominator == 2**20, 'complete fixed grid')
    cost = F(total, denominator)
    need(1/HIGH**2 - 180/(S+2*LOW)**2 - cost > 5, 'joint Schur complement in w')
    # After returning to mass nine, take the weighted Cauchy bound on the sum of
    # the w and mean-zero angle contributions to the half-difference vector.
    radial_lipschitz_sq = 9*(S*S/(4*LOW)+S/2)/(S+2*LOW)**3
    angular_lipschitz_sq = 9*F(113,100)**2/(2*(S+2*LOW))
    need(radial_lipschitz_sq/5+2*angular_lipschitz_sq < F(33,25)**2, 'mass-nine amplitude transfer')
    return {'phase_bins':512, 'successive_bin_index_differences':list(map(int,phase_input.splitlines()[0].split()[1:])),
            'local_hessian_boxes':cells, 'cyclic_schur_sum_upper':str(cost),
            'directional_hessian_weights':{'inner_squared_half_difference':'5','mean_zero_angle':'1/2'}}


def tubes(repo):
    tube = Tube(repo)
    tube.phase.verify_inputs()
    records = {}
    for mode, eps in EPS.items():
        data = tube.data(mode)
        counts = Counter()
        cursor = F(-1)
        for row in read(HERE/(mode+'.json')):
            lo, hi, w, *rest = row
            lo, hi = F(lo), F(hi)
            need(lo == cursor and lo < hi <= 1 and len(rest)==4, 'complete ordered phase interval')
            wanted = rest[-1]
            result = tube.accept(mode, data, lo, hi, list(map(F,[w,*rest[:-1]])))
            need(wanted in ('candidate_basin','energy_excluded') and result == wanted, 'accepted cell')
            counts[wanted] += 1
            cursor = hi
        need(cursor == 1, 'whole phase period')
        records[mode] = {'normalized_shape_tolerance':str(eps), 'phase_cells':sum(counts.values()),
                         'candidate_basin':counts['candidate_basin'], 'energy_excluded':counts['energy_excluded']}
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--repo', type=Path, default=HERE.parent)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    repo = args.repo.resolve()
    prerequisites(repo)
    report = {'passed':True, 'normalization_mass':9, 'outer_scale_constraint':str(S),
              'all_global_maximizers_satisfy_domain':True, 'normalized_w_domain':[str(LOW),str(HIGH)],
              'remaining_shape_parameters':10, 'mass_nine_half_difference_residual_multiplier':'33/25',
              'joint_hessian':hessian(args.compiler), 'full_inner_separation_phase_tubes':tubes(repo),
              'unresolved_inside_the_two_tubes':0}
    output = HERE / 'verification_report.json'
    data = json.dumps(report, indent=2)+'\n'
    if args.record_report:
        output.write_text(data, encoding='ascii', newline='\n')
    else:
        need(report == read(output), 'frozen verification report')
    print(data,end='')


if __name__ == '__main__':
    main()
