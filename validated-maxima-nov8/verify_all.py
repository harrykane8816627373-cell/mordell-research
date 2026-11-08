"""Replay Nov8 root and full-Hessian certificates using pinned integer code.

Only the standard library is used. All checks execute with Python -O.
"""
from pathlib import Path
from fractions import Fraction
from math import prod
import argparse
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCIES = {
    'validated-maxima/verify_all.py':'4875024b8eacdc48fefdb83b447655270c2fbf8bcf785ddf53337526d5dcc6e0',
    'validated-maxima/dyadic.py':'94ab5e20011cf2742566e7aaf0af71c41e7fc4254e2320dfb4b3af97fcf07ccc',
    'validated-maxima/PROOF.md':'d3bd80d74bc73873f70fca6ef20f1f48770b5f1a763912e146db7472bd216341',
    'validated-maxima/discover.py':'6243dccedf62b5939b483bc8e498299d3aeff764c1df0e0dba4a9b2e0e15e40d',
    'finite-search-nov5/integer_witnesses.json':'83cc7ee2a03b549f1024e3f683afcdcb13298146293036cc75b54809ff7b2d82',
    'symmetry-breaking/integer_witnesses.json':'73e5f42ef78db59dd0ff3e05227caf75eb3834a02486c509767c8534c3c954d5',
}
EXPECTED_COUNTS = [25,35,42,46]


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def integer_value(points):
    n = len(points)
    need(n >= 2 and all(len(p) == 2 and all(type(q) is int for q in p) for p in points), 'integer seed coordinates')
    sx,sy = sum(p[0] for p in points),sum(p[1] for p in points)
    distances = [(x-u)**2+(y-v)**2 for i,(x,y) in enumerate(points) for u,v in points[i+1:]]
    mass = n*sum(x*x+y*y for x,y in points)-sx*sx-sy*sy
    need(min(distances) > 0 and mass == sum(distances), 'exact centered pair-mass identity')
    exponent = n*(n-1)//2
    return Fraction(n**(2*exponent)*prod(distances),mass**exponent)


def decimal_bounds(low, high, digits=12):
    need(0 < low <= high, 'positive ordered exact ratio interval')
    scale = 10**digits
    a = low.numerator*scale//low.denominator
    b = -((-high.numerator*scale)//high.denominator)
    def fmt(integer):
        return str(integer//scale)+'.'+str(integer%scale).zfill(digits)
    return [fmt(a),fmt(b)]


def load_prior():
    for name, expected in DEPENDENCIES.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, 'pinned dependency '+name)
    # The imported verifier explicitly imports the module named dyadic. Bind that
    # name to the checked file, avoiding any ambient site-package of that name.
    spec = importlib.util.spec_from_file_location('dyadic', ROOT/'validated-maxima/dyadic.py')
    dyadic = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dyadic)
    old = sys.modules.get('dyadic')
    sys.modules['dyadic'] = dyadic
    try:
        spec = importlib.util.spec_from_file_location('nov8_pinned_prior_verifier', ROOT/'validated-maxima/verify_all.py')
        prior = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(prior)
    finally:
        if old is None:
            del sys.modules['dyadic']
        else:
            sys.modules['dyadic'] = old
    need(prior.PRECISION == dyadic.PRECISION == 384, 'fixed exact dyadic precision')
    return prior


def verify_data(data, prior, expected_counts=EXPECTED_COUNTS, arithmetic_checks=True):
    need(data['arithmetic_precision'] == 384 and data['center_exponent'] == 256
         and data['matrix_exponent'] == 192 and data['radius_exponent'] == 80,
         'fixed arithmetic and dyadic lattice parameters')
    counts = [r['n'] for r in data['rows']]
    need(counts == expected_counts, 'complete explicit Nov8 point counts')
    raw = (ROOT/'finite-search-nov5/integer_witnesses.json').read_bytes()
    seeds = {r['n']:r['integer_points'] for r in json.loads(raw)['rows']}
    baseline = {r['n']:r['integer_points'] for r in json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())['rows']}
    results = []
    for row in data['rows']:
        n = row['n']
        need(type(n) is int and row['discovery_integer_points'] == seeds[n], 'exact pinned integer seed n='+str(n))
        need(len(row['discovery_integer_points']) == n, 'actual seed cardinality')
        result = prior.verify_row(row, data)
        low,high = map(Fraction,result['delta_interval'])
        old_value,new_value = integer_value(baseline[n]),integer_value(seeds[n])
        need(low > new_value > old_value, 'certified root exceeds new and old exact integer lower bounds')
        result['root_delta_over_old_integer_value_interval'] = decimal_bounds(low/old_value,high/old_value)
        result['root_delta_over_discovery_integer_value_interval'] = decimal_bounds(low/new_value,high/new_value,22)
        result['certified_root_strictly_exceeds_discovery_integer_value'] = True
        result['full_Lagrangian_Hessian_inertia'] = {'negative':2*n-1,'zero':1,'positive':0} if row['classification'] == 'strict_local_maximum_mod_rotation' else None
        result['sphere_tangent_Lagrangian_Hessian_inertia'] = {'negative':2*n-2,'zero':1,'positive':0} if row['classification'] == 'strict_local_maximum_mod_rotation' else None
        results.append(result)
    return {'passed':True, 'date':'2026-11-08',
            'arithmetic':'outward integer dyadic intervals at 384 bits',
            'checks_survive_optimization':True, 'point_counts':counts,
            'stationary_orbits_certified':len(results),
            'strict_local_maxima_mod_rotation':[r['n'] for r in results if r['classification'] == 'strict_local_maximum_mod_rotation'],
            'unstable_stationary_orbits':[r['n'] for r in results if r['classification'] == 'unstable_stationary_orbit'],
            'independent_arithmetic_checks':prior.independent_checks() if arithmetic_checks else None,
            'global_optimality_claimed':False, 'prior_certificates_rerun_or_modified':False,
            'pinned_dependency_sha256':DEPENDENCIES, 'rows':results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    need(manifest['dependencies'] == DEPENDENCIES, 'transparent pinned dependency set')
    for name, expected in manifest['files'].items():
        need('/' not in name and '\\' not in name and name not in ('.','..'), 'safe local fixed input')
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, 'fixed input '+name)
    prior = load_prior()
    data = json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    report = verify_data(data, prior)
    report['source_sha256'] = manifest['files']
    target = HERE/'verification_report.json'
    if args.record_report:
        target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    else:
        need(report == json.loads(target.read_text(encoding='utf-8')), 'fresh exact results reproduce the frozen report')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
