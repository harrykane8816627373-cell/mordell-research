"""Replay a finite-depth integer interval cover of the radial half-turn chart."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
DEFAULT_DEPTH = 25


def require(holds, message):
    if not holds:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_bytes())


def prerequisites(repo):
    manifest = read(HERE / 'source_manifest.json')
    for base, group in ((HERE, 'files'), (repo, 'dependencies')):
        for name, digest in manifest[group].items():
            require(sha256((base / name).read_bytes()).hexdigest() == digest, 'fixed input ' + name)

    def proof(folder):
        r = read(repo / folder / 'verification_report.json')
        require(r['passed'] and r['applies_to_every_global_maximizer'] and r['normalization_mass'] == 9,
                'accepted global geometry from ' + folder)
        return r

    radii = proof('ten-point-angular-radii')
    require((radii['exact_inner_count'], radii['exact_outer_count']) == (2, 8)
            and tuple(map(F, radii['outer_squared_radius_open_interval'])) == (F(89,100),F(1271,1000)),
            'fixed shell populations and outer radial limits')
    inner = proof('ten-point-inner-edge')
    require(inner['unresolved']==0 and
            tuple(F(x)/4 for x in inner['inner_squared_distance_open_interval']) == (F(43,400),F(53,400)),
            'fixed central-pair radius limits')
    dispersion = read(repo / 'ten-point-half-turn-dispersion' / 'verification_report.json')
    require(dispersion['passed'] and dispersion['normalization_mass']==9
            and dispersion['assumption']=='half-turn symmetry of a global maximizer'
            and F(dispersion['outer_squared_radius_span_strict_lower']) == F(47,2000),
            'fixed half-turn radial span')
    gap = proof('ten-point-outer-gap')
    lo,hi = map(F, gap['outer_successive_angle_cosine_open_interval'])
    low, high = 2*lo**2-1, 2*hi**2-1
    require((low,high) == (F(-1519,5000),F(271,1250)) and low<0<high,
            'doubled consecutive-angle cosine limits')
    require(1-max(low**2,high**2) > F(19,20)**2, 'adjacent determinant lower bound')
    require(max(2*low**2-1,2*high**2-1) < F(-81,100), 'doubled two-step cosine limit')
    require(1-F(81,100)**2 < (F(5,8)*F(19,20))**2 and F(89,100)>F(17,20),
            'the minimum-radius mass equation has derivative strictly above three eighths')
    minimum = read(repo / 'maximizer-domain' / 'verification_report.json')
    require(minimum['passed'] and minimum['symmetry_assumed'] is False
            and F(minimum['pair_squared_distance_strict_lower']) == F(2,5), 'global pair separation')

    d2 = repo / 'ten-point-d2-cover' / 'verify_all.py'
    spec = importlib.util.spec_from_file_location('halfturn_radial_witness',d2)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    require(verifier.WITNESS == 469765409, 'fixed mass-nine feasible objective threshold')
    verifier.check_witness()

    resultant = read(HERE / 'jacobian_resultant.json')
    out = resultant['output']
    require(resultant['operation_id'] == 'polynomial.multivariate.resultant.compute'
            and out['elimination_variable'] == 'x', 'Jacobian resultant operation')

    def terms(source):
        return {tuple(t['exponents']):F(t['coefficient']['num'])/F(t['coefficient']['den'])
                for t in source['polynomial']['terms']}

    require(out['left']['variables']==['x','l','A','B'] and out['right']['variables']==['x','l','A','B']
            and terms(out['left'])=={(2,0,0,0):1,(0,0,1,0):-1}
            and terms(out['right'])=={(2,0,0,0):1,(1,1,0,0):-2,(0,2,0,0):1,(0,0,0,1):-1},
            'source-bound radical elimination equations')
    polynomial = out['resultant']['value']
    # Product of (A+l²-B)-2*l*x at x and -x, with x²=A.
    expected = {(4,0,0):1,(2,1,0):-2,(2,0,1):-2,(0,2,0):1,(0,1,1):-2,(0,0,2):1}
    require(out['resultant']['kind']=='POLYNOMIAL' and polynomial['variables']==['l','A','B']
            and terms(polynomial)==expected, 'independent radical resultant expansion')


def replay(compiler, depth):
    with TemporaryDirectory(prefix='mordell-half-turn-') as scratch:
        binary = Path(scratch) / 'cover.exe'
        subprocess.run([compiler,'-std=c++20','-O3',str(HERE/'cover.cpp'),'-o',str(binary)],
                       check=True, timeout=120)
        process = subprocess.run([str(binary),str(depth),'0'], check=True, timeout=1200,
                                 capture_output=True, text=True, encoding='ascii')
    rows = {}
    labels = ('visited','split','geometry','force','value','unresolved','multidimensional_force')
    for line in process.stdout.splitlines():
        side, counts = line.split(':')
        require(side in ('-1','1') and side not in rows, 'distinct half-circle charts')
        numbers = list(map(int, counts.split()))
        require(len(numbers) == len(labels) and all(x >= 0 for x in numbers), 'complete count row')
        row = dict(zip(labels,numbers))
        require(row['visited'] == 2*row['split']+1 == sum(numbers[1:]), 'binary subdivision totals')
        rows[side] = row
    require(set(rows) == {'-1','1'}, 'both charts covered to the declared depth')
    unresolved = sum(x['unresolved'] for x in rows.values())
    return {'passed':True, 'normalization_mass':9,
            'scope':'finite-depth partial cover conditional on half-turn symmetry',
            'inner_squared_radius_open_interval':['43/400','53/400'],
            'outer_squared_radius_span_strict_lower':'47/2000',
            'radii_profile_mass_derivative_strict_lower':'3/8', 'max_depth':depth,
            'mass_interval_newton_contraction':True,
            'multidimensional_mean_value_contraction':True,
            'multidimensional_force_exclusions':sum(x['multidimensional_force'] for x in rows.values()),
            'jacobian_radical_resultant_checked':True,
            'nodes':sum(x['visited'] for x in rows.values()), 'unresolved':unresolved,
            'whole_half_turn_branch_closed':unresolved == 0, 'charts':rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=HERE.parent)
    parser.add_argument('--compiler',default='g++')
    parser.add_argument('--depth',type=int,default=DEFAULT_DEPTH)
    parser.add_argument('--record-report',action='store_true')
    args = parser.parse_args()
    prerequisites(args.repo.resolve())
    require(0 <= args.depth <= 35, 'admissible finite depth')
    require(not args.record_report or args.depth == DEFAULT_DEPTH, 'freeze only the declared depth')
    result = replay(args.compiler,args.depth)
    path = HERE / 'verification_report.json'
    output = json.dumps(result,indent=2)+'\n'
    if args.record_report:
        path.write_text(output, encoding='ascii',newline='\n')
    elif args.depth == DEFAULT_DEPTH:
        require(result == read(path), 'matches saved finite-depth report')
    print(output,end='')


if __name__ == '__main__':
    main()
