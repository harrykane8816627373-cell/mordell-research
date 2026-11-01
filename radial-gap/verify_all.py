"""Replay the radial gap, its local-energy input and the supplied unified package."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def safe_relative(name):
    path = PurePosixPath(name)
    need(not path.is_absolute() and '..' not in path.parts and '\\' not in name
         and ':' not in name and path.as_posix() == name, 'safe fixed path')
    return Path(*path.parts)


def unpack_input(archive_path, destination):
    prefix = 'mordell-eight-point-unified-2026-11-01/'
    with zipfile.ZipFile(archive_path) as archive:
        members = archive.namelist()
        need(len(members) == len(set(members)), 'unique attached ZIP members')
        manifest = archive.read(prefix+'manifest.sha256').decode('utf-8-sig')
        checked = set()
        for line in manifest.splitlines():
            if not line.strip():
                continue
            digest, name = line.split(None, 1)
            name = name.lstrip('* ')
            # sha256sum records the original names with one harmless './'.
            if name.startswith('./'):
                name = name[2:]
            safe_relative(name)
            need(name not in checked, 'unique attached manifest rows')
            checked.add(name)
            need(hashlib.sha256(archive.read(prefix+name)).hexdigest() == digest,
                 'attached original hash '+name)
        files = {info.filename[len(prefix):] for info in archive.infolist()
                 if not info.is_dir() and info.filename.startswith(prefix)}
        need(files == checked | {'manifest.sha256'}, 'attached manifest covers all original files')
        for info in archive.infolist():
            need(info.filename.startswith(prefix), 'attached ZIP single expected root')
            if info.is_dir():
                continue
            name = info.filename[len(prefix):]
            target = destination/safe_relative(name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(archive.read(info))
    need(len(checked) == 49, '49 original attached hashes checked')
    return len(checked)


def main():
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    paths = []
    for row in manifest['files']:
        name = row['path']
        path = ROOT/safe_relative(name)
        data = path.read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
             'fixed bytes '+name)
        paths.append(name)
    need(len(paths) == len(set(paths)), 'unique fixed inputs')
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-radial-', dir=ROOT) as temp:
        stage = Path(temp)
        for name in paths:
            target = stage/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, target)

        def run(script, *arguments):
            proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/script), *arguments],
                                  cwd=stage, env=env, capture_output=True, text=True,
                                  encoding='utf-8', timeout=180)
            need(proc.returncode == 0, script+': '+proc.stdout[-2000:]+proc.stderr[-3000:])
            return proc.stdout

        run('eight-point/verify_all.py')
        print('PASS inherited two-entry objective/tube suite, including actual spectra and witnesses', flush=True)
        audit = json.loads((stage/'radial-gap/input_audit.json').read_text(encoding='utf-8'))
        original = stage/'radial-gap/inherited/mordell-eight-point-unified-2026-11-01.zip'
        need(hashlib.sha256(original.read_bytes()).hexdigest() == audit['input_sha256'],
             'supplied original archive kept byte-for-byte')
        hash_count = unpack_input(original, stage/'attached-unified')
        log = run('attached-unified/verification/verify_all.py', '--legacy')
        need('Exact scripted checks: 14/14 passed' in log, '14 supplied standard-library entries passed')
        print('PASS supplied unified archive: 49 pinned hashes, 14 standard-library entries', flush=True)
        run('radial-gap/verify_gap.py')
        result = json.loads((stage/'radial-gap/verification.json').read_text(encoding='utf-8'))
        need(result['passed'] is True and result['full_eight_point_exclusion_proved'] is False,
             'new exact report and coverage boundary')
        need(result['actual_nonstationary_witness']['is_stationary'] is False,
             'witness is not a false stationary counterexample')
        print('PASS new radial/angular constants, actual-root witness and two rejected strengthenings', flush=True)
    print(json.dumps({'passed': True, 'fixed_inputs': len(paths),
        'attached_manifest_entries': hash_count, 'attached_standard_library_entries': 14,
        'inherited_objective_entries': 2, 'new_exact_entries': 1,
        'radial_threshold': '1/1200', 'log_gap_coefficient': '1/32',
        'optional_sympy_root_pairing_tests_replayed': False,
        'all_imported_analytic_theorems_independently_reproved': False,
        'full_eight_point_exclusion_proved': False}, indent=2))


if __name__ == '__main__':
    main()
