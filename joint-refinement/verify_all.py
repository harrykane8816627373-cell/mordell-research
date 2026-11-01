"""Replay pinned old inputs, new exact arithmetic and the incomplete continuation."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import os
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
    parser.add_argument('--prior-already-replayed', action='store_true',
        help='only when the caller already replayed all joint-closure inputs and certificates')
    parser.add_argument('--record-reports', action='store_true',
        help='write newly computed reports only after every requested replay passed')
    parser.add_argument('--check-fixed-inputs', action='store_true',
        help='check pinned bytes only; do not claim that proof replay ran')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    inputs = {}
    for row in manifest['files']:
        name = row['path']
        path = PurePosixPath(name)
        need(not path.is_absolute() and '..' not in path.parts and ':' not in name
             and '\\' not in name and path.as_posix() == name and name not in inputs,
             'safe unique fixed input')
        data = (ROOT/name).read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
             'pinned source '+name)
        inputs[name] = data
    if args.check_fixed_inputs:
        print(json.dumps({'fixed_inputs_checked':len(inputs),'proof_replay_performed':False}))
        return
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-refinement-', dir=ROOT) as temp:
        stage = Path(temp)
        for name, data in inputs.items():
            path = stage/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)

        def run(script, *arguments):
            result = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/script), *arguments],
                cwd=stage, env=env, capture_output=True, text=True, encoding='utf-8', timeout=1200)
            need(result.returncode == 0, script+': '+result.stdout[-3000:]+result.stderr[-3000:])
            return result.stdout

        if not args.prior_already_replayed:
            run('joint-closure/verify_all.py')
            print('PASS 20 prior entries including the frozen 60021-node tree', flush=True)
        entries = [('verify.py', [], 'arithmetic.json'),
            ('check.py', ['joint-refinement/certificates/continuation.json.gz', '--base-already-replayed'], 'cover.json')]
        results = {}
        for script, arguments, report in entries:
            actual = json.loads(run('joint-refinement/'+script, *arguments))
            need(actual['passed'] is True, 'exact checks passed '+report)
            if not args.record_reports:
                expected = json.loads((stage/'joint-refinement/reports'/report).read_text(encoding='utf-8'))
                need(actual == expected, 'deterministic exact report '+report)
            need(actual['full_eight_point_exclusion_proved'] is False, 'current coverage remains incomplete')
            results[report] = actual
            print('PASS '+script, flush=True)
    cover = results['cover.json']
    if args.record_reports:
        (HERE/'reports').mkdir(exist_ok=True)
        for name, result in results.items():
            (HERE/'reports'/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed': True, 'new_exact_entries': 2,
        'prior_entries_replayed': 0 if args.prior_already_replayed else 20,
        'combined_tree_nodes': cover['combined_tree_nodes'],
        'combined_excluded_leaves': cover['combined_excluded_leaves'],
        'new_Taylor_excluded_leaves': cover['new_Taylor_excluded_leaves'],
        'unresolved_regions': cover['unresolved_regions'],
        'complete_cover': False, 'full_eight_point_exclusion_proved': False}, indent=2))


if __name__ == '__main__':
    main()
