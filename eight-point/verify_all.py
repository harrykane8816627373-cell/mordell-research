"""Verify pinned inputs and replay the current certificate in a fresh copy."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def main():
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    paths = []
    for row in manifest['files']:
        name = row['path']
        rel = PurePosixPath(name)
        require(not rel.is_absolute() and '\\' not in name and ':' not in name
                and '..' not in rel.parts and name == rel.as_posix(), 'safe fixed path')
        data = (ROOT/Path(*rel.parts)).read_bytes()
        require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                'fixed bytes '+name)
        paths.append(name)
    require(len(paths) == len(set(paths)), 'unique fixed inputs')
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-eight-', dir=ROOT) as temp:
        stage = Path(temp)
        for name in paths:
            dest = stage/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, dest)
        result = subprocess.run(
            [sys.executable, '-X', 'int_max_str_digits=0', '-B', '-S', '-O',
             str(stage/'eight-point/verify_tube.py')], cwd=stage, env=env,
            capture_output=True, text=True, encoding='utf-8', timeout=120)
        require(result.returncode == 0, result.stdout[-2000:]+result.stderr[-3000:])
        actual = json.loads((stage/'eight-point/tube_verification.json').read_text(encoding='utf-8'))
        expected = json.loads((HERE/'tube_verification.json').read_text(encoding='utf-8'))
        require(actual == expected and actual['passed'] is True, 'deterministic exact report')
        require(actual['full_eight_point_exclusion_proved'] is False, 'incomplete global coverage')
        require(actual['tube_certificate']['p7_restricted'] is False, 'p7-free scope')
        provenance = json.loads((HERE/'PROVENANCE.json').read_text(encoding='utf-8'))
        require(provenance['full_eight_point_exclusion_proved'] is False, 'provenance coverage boundary')
        print(result.stdout.strip())
    print(json.dumps({'passed': True, 'fixed_inputs': len(paths),
        'legacy_actual_forces_and_complete_spectrum_recomputed': True,
        'new_p7_free_tube_certified': True, 'unsafe_parameter_rejections': 3,
        'full_eight_point_exclusion_proved': False}, indent=2))


if __name__ == '__main__':
    main()
