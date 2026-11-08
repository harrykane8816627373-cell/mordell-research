"""Optional mpmath candidate generation for the pinned Nov8 seed set."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
import time
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PINNED_DISCOVERY_SHA256 = '6243dccedf62b5939b483bc8e498299d3aeff764c1df0e0dba4a9b2e0e15e40d'
PINNED_SEEDS_SHA256 = '83cc7ee2a03b549f1024e3f683afcdcb13298146293036cc75b54809ff7b2d82'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, choices=(25,35,42,46), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = ROOT/'validated-maxima/discover.py'
    if hashlib.sha256(source.read_bytes()).hexdigest() != PINNED_DISCOVERY_SHA256:
        raise ArithmeticError('pinned prior high precision discovery source')
    seed_path = ROOT/'finite-search-nov5/integer_witnesses.json'
    seed_raw = seed_path.read_bytes()
    if hashlib.sha256(seed_raw).hexdigest() != PINNED_SEEDS_SHA256:
        raise ArithmeticError('pinned Nov5 integer seed set')
    spec = importlib.util.spec_from_file_location('nov8_prior_discovery', source)
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    mp.mp.dps = 180
    row = next(r for r in json.loads(seed_raw)['rows'] if r['n'] == args.n)
    start = time.monotonic()
    result = prior.discover(row, 'finite-search-nov5/integer_witnesses.json, n='+str(args.n))
    report = {'date':'2026-11-08', 'arithmetic_precision':384, 'center_exponent':256,
              'matrix_exponent':192, 'radius_exponent':80, 'rows':[result],
              'generation_metadata':{'python':sys.version.split()[0], 'mpmath':mp.__version__,
                                     'decimal_precision':mp.mp.dps,
                                     'prior_discovery_sha256':PINNED_DISCOVERY_SHA256,
                                     'seed_file_sha256':PINNED_SEEDS_SHA256,
                                     'elapsed_seconds':round(time.monotonic()-start,3),
                                     'floating_generation_is_not_acceptance':True}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'n':args.n,'classification_candidate':result['classification'],
                      'newton_steps':result['discovery_newton_steps'],
                      'output':str(args.output),'elapsed_seconds':report['generation_metadata']['elapsed_seconds']}),flush=True)


if __name__ == '__main__':
    main()
