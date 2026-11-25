"""Rebuild the four-angle half-turn cover with outward dyadic arithmetic."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import product
import json
from math import comb, floor, ceil, prod
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BITS = 43
DEFAULT_DEPTH = 30
MIDPOINT_CEILING = F(1,18)


def read(path):
    return json.loads(path.read_bytes())


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def load(folder):
    path = ROOT/folder/'verify_all.py'
    spec = importlib.util.spec_from_file_location(folder, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def inputs():
    manifest = read(HERE/'source_manifest.json')
    for folder, entries in ((HERE, manifest['files']), (ROOT, manifest['dependencies'])):
        for name, digest in entries.items():
            require(sha256((folder/name).read_bytes()).hexdigest()==digest, 'fixed input '+name)
    prior = load('ten-point-half-turn-radial-cover')
    prior.prerequisites(ROOT)
    radial = load('ten-point-half-turn-radial-concavity')
    result = read(ROOT/'ten-point-half-turn-radial-concavity/verification_report.json')
    require(result['passed'] and result['cover']['unresolved']==0 and
            F(result['negative_log_radial_hessian_strict_lower'])==F(1,4), 'strict radial concavity')
    a,b = radial.atan_interval(5)
    c,d = radial.atan_interval(239)
    require(F(62831853071795,10**13) < 32*a-8*d < 32*b-8*c < F(62831853071797,10**13),
            'two-pi bracket')
    neighbor = read(ROOT/'ten-point-half-turn-angle-neighborhoods/verification_report.json')
    require(neighbor['passed'] and
            {k:F(v['squared_image_angle_tolerance']) for k,v in neighbor['angle_neighborhoods'].items()}
            == {'candidate': F(11,500), 'saddle': F(3,125)}, 'rigorous angular orbit neighborhoods')
    local = read(ROOT/'ten-point-quartic-midpoints/verification_report.json')
    saddle = read(ROOT/'ten-point-saddle-exclusion/verification_report.json')
    require(local['passed'] and local['normalization_mass']==9 and
            F(local['half_difference_radius'])==F(7,250) and saddle['passed'] and
            F(saddle['mass_nine_excluded_orbit_radius'])==F(21,500), 'two mass-nine orbit tests')
    sixth_order, energy_coefficient, projection_coefficient = check_midpoint_bridge(MIDPOINT_CEILING)
    engine = load('ten-point-half-turn-phase-sweep').Sweep()
    data = [MIDPOINT_CEILING.denominator, sixth_order, energy_coefficient, projection_coefficient]
    def append(interval):
        lo = engine.box.endpoint(interval,False)
        hi = engine.box.endpoint(interval)
        data.extend([floor(lo*(1<<BITS)),ceil(hi*(1<<BITS))])
    for name in ('candidate','saddle'):
        for direction in engine.directions[name][1:]:
            for v in direction:
                append(v)
    for name in ('candidate','saddle'):
        for i in (0,2,3,4,5):
            radius = engine.op.norm2(engine.box.points[name][i])
            lo,hi = (F(43,400),F(53,400)) if i==0 else (F(89,100),F(1271,1000))
            require(lo<engine.box.endpoint(radius,False)<=engine.box.endpoint(radius)<hi, 'anchor radial bands')
            append(radius)
    return ' '.join(map(str,data))+'\n'


def check_midpoint_bridge(mu):
    proof = read(ROOT/'ten-point-pair-perturbation/verification_report.json')
    require(proof['passed'] and proof['applies_to_every_global_maximizer'] and
            proof['small_midpoint_quartic_energy_coefficient']==74 and
            proof['small_midpoint_quadratic_sixth_order_coefficient']==846 and
            F(proof['projected_cross_pair_distance_strict_lower'])==F(4,7) and
            F(proof['one_cross_pair_distance_strict_lower'])==F(99,100), 'projection estimates')
    d, amp, long = F(4,7), F(9359,10000), F(99,100)
    eta_squared = 2*mu**2
    require(mu>0 and mu.numerator==1 and eta_squared<d**2 and mu**2<amp**2, 'chosen midpoint ceiling')
    a, l, u = F(53,400), F(89,100), F(1271,1000)
    cutoff = F(9,2)-a-mu**2
    vertices = [(edges,cutoff-sum(edges)) for edges in product((l,u),repeat=3)]
    require(min(prod(edges)*last for edges,last in vertices if l<=last<=u)>F(43,40)**4,
            'repaired radii belong to the log domain')
    require(1/(amp**2-mu**2)<F(6,5), 'angular and logarithmic repair displacement')
    require(amp**2+F(43,400)>long**2 and
            2*F(1128,1000)*sum(1/(v*(v**2-2*mu**2)) for v in (d,long))<15,
            'one long edge for every pair of half differences')
    cos_hi = F(39,50)+2*F(3,250)
    require(2/(2-2*cos_hi**2)<3 and
            4*a*u/(amp**2-a)**2<3 and 5*3+10*u<28, 'repair Hessian operator norm')
    # At a scalar constrained extremum the coordinates satisfy a depressed
    # cubic. In a three-value pattern the excess multiplicities force a zero.
    fourth_moment = max(F((5-k)**2,25*k)+F(k*k,25*(5-k)) for k in (1,2))
    require(fourth_moment==F(13,20)>max(F(1,2),F(1,4)), 'centered fourth moment')
    edge_moment = 5*fourth_moment+3
    require(edge_moment==F(25,4), 'pairwise fourth moment')
    raw_gradient, displacement = F(84), F(137,100)
    require(raw_gradient**2>15**2*(25*fourth_moment+15) and
            displacement**2>2*F(6,5)**2*fourth_moment,
            'raw gradient and repair displacement')
    require(raw_gradient+28*displacement<123, 'force residual after repair')
    projection = edge_moment*sum((1+eta_squared/v**2)/(v**4*(1-eta_squared/v**2)) for v in (d,long))
    energy_coefficient = ceil(projection+raw_gradient*displacement+14*displacement**2)
    require(4*F(43,400)>F(13,20)**2 and l>F(9,10)**2 and amp>F(9,10),
            'inverse square-root estimate for repaired and unscaled radii')
    require(fourth_moment*(1/(2*amp)+F(1128,1000)*F(3,5))**2<1,
            'representative displacement after repair')
    require(mu**2<F(7,250), 'nonempty candidate orbit test')
    sixth_order = ceil(4*edge_moment*sum(1/(v**2*(v**2-eta_squared)**2) for v in (d,long)))
    return sixth_order, energy_coefficient, ceil(projection)


def replay(compiler, depth, text):
    with TemporaryDirectory(prefix='mordell-four-angle-') as work:
        binary = Path(work)/'cover.exe'
        subprocess.run([compiler, '-std=c++20', '-O3', '-ffp-contract=off', str(HERE/'cover.cpp'),
                        '-o', str(binary)], check=True, timeout=120)
        outputs = [subprocess.run([str(binary),str(depth),*mode], input=text, capture_output=True,
                                  encoding='ascii', check=True, timeout=1200).stdout
                   for mode in ([],['midpoints'])]
    labels = ('root_boxes','visited','bisected','geometry_excluded','energy_excluded',
              'force_excluded','krawczyk_excluded','candidate_basin','saddle_excluded',
              'unresolved','midpoint_stationarity_forces_zero','radial_subcovers','certified_radial_inclusions',
              'energy_radial_localizations','local_hessian_refinements','sharper_radial_inclusions',
              'radial_subcover_nodes','maximum_depth')
    def counts(output):
        numbers = list(map(int,output.split()))
        require(len(numbers)==len(labels) and all(x>=0 for x in numbers), 'counts')
        rows = dict(zip(labels,numbers))
        require(rows['root_boxes']==comb(32,4) and rows['visited']==rows['root_boxes']+2*rows['bisected']
                and rows['visited']==sum(numbers[2:12]) and rows['maximum_depth']<=depth, 'cover accounting')
        return rows
    exact, perturbed = map(counts,outputs)
    sixth_order, energy_coefficient, projection_coefficient = map(int,text.split()[1:4])
    return {'passed': True, 'mass':9, 'assumption':'half-turn symmetry of a global maximizer',
            'depth_limit':depth, 'counts':exact, 'half_turn_classification_complete':exact['unresolved']==0,
            'general_configurations_with_small_midpoints':{'norm_ceiling':str(MIDPOINT_CEILING),
                                                          'sixth_order_coefficient':sixth_order,
                                                          'repaired_energy_coefficient':energy_coefficient,
                                                          'projection_energy_coefficient':projection_coefficient,
                                                          'counts':perturbed,
                                                          'complete':perturbed['unresolved']==0}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--depth', type=int, default=DEFAULT_DEPTH)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    if not 0<=args.depth<=30:
        parser.error('--depth must be between 0 and 30')
    report = replay(args.compiler, args.depth, inputs())
    if args.depth==DEFAULT_DEPTH:
        require(report['half_turn_classification_complete'], 'whole half-turn domain covered')
        require(report['general_configurations_with_small_midpoints']['complete'],
                'whole midpoint neighborhood covered')
    saved = HERE/'verification_report.json'
    output = json.dumps(report, indent=2)+'\n'
    if args.record_report:
        require(args.depth==DEFAULT_DEPTH, 'default depth report')
        saved.write_bytes(output.encode())
    elif args.depth==DEFAULT_DEPTH:
        require(report==read(saved), 'frozen replay counts')
    print(output,end='')


if __name__ == '__main__':
    main()
