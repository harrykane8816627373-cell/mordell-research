"""Replay the equivalent model, a physical orbit box and the incomplete cover."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-prior', action='store_true', help='caller already replayed radial-gap')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    names = []
    for row in manifest['files']:
        name = row['path']
        rel = PurePosixPath(name)
        need(not rel.is_absolute() and '..' not in rel.parts and ':' not in name
             and '\\' not in name and rel.as_posix() == name, 'safe fixed input')
        data = (ROOT/name).read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
             'fixed bytes '+name)
        names.append(name)
    need(len(names) == len(set(names)), 'unique fixed input names')
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-joint-', dir=ROOT) as temp:
        stage = Path(temp)
        for name in names:
            target = stage/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, target)

        def run(script, *arguments):
            result = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/script), *arguments],
                cwd=stage, env=env, capture_output=True, text=True, encoding='utf-8', timeout=600)
            need(result.returncode == 0, script+': '+result.stdout[-2000:]+result.stderr[-3000:])
            return result.stdout

        if not args.skip_prior:
            run('radial-gap/verify_all.py')
            print('PASS 17 prior radial/input entries', flush=True)
        entries = [('verify_model.py', [], 'model.json'),
                   ('verify_orbit.py', ['joint-closure/certificates/critical_orbit.json'], 'orbit.json'),
                   ('check_cover.py', ['joint-closure/certificates/partial_cover.json.gz'], 'cover.json')]
        results = {}
        for script, arguments, report in entries:
            actual = json.loads(run('joint-closure/'+script, *arguments))
            expected = json.loads((HERE/'reports'/report).read_text(encoding='utf-8'))
            need(actual == expected and actual['passed'] is True, 'deterministic exact report '+report)
            need(actual['full_eight_point_exclusion_proved'] is False, 'no complete exclusion claim')
            results[report] = actual
            print('PASS '+script, flush=True)
        cover = results['cover.json']
        need(cover['complete_cover'] is False and cover['counts']['unresolved'] == 21,
             'all 21 unresolved regions remain visible')
        need(results['orbit.json']['first_four_mode_block_positive_definite'] is True,
             'low-mode versus full-Hessian distinction')
    print(json.dumps({'passed': True, 'new_exact_entries': 3,
        'prior_entries_replayed': 0 if args.skip_prior else 17,
        'polynomial_matrix_degree_upper': 14,
        'unique_exact_physical_stationary_orbit_certified': True,
        'partition_nodes': cover['tree_nodes'], 'excluded_leaves': 29990,
        'unresolved_coordinate_regions': 21, 'complete_cover': False,
        'full_eight_point_exclusion_proved': False}, indent=2))


if __name__ == '__main__':
    main()
