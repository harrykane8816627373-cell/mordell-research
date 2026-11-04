"""Replay the pinned, attributed nine-point proof from the full initial domain."""
from pathlib import Path, PurePosixPath
from fractions import Fraction
import argparse, hashlib, json, os, platform
import shutil, subprocess, sys, tempfile, time

HERE = Path(__file__).resolve().parent
UPSTREAM_COMMIT = '31cc934b7dfa0ef2308566a1bce48cbfb2e01559'
COEFFICIENT_SHA256 = '1d1cef17749fcb242e6f09dffe3c26183ae672fe4707b2bad70bb0b12b0e9d88'
# Filled from the byte-for-byte Git-blob-checked import, not from a result report.
UPSTREAM_META_SHA256 = '86e58f21339396e92c00ba7c9ac677e4ca24a5e2c5ba50f7c6eb6a6663a3cbb0'
SCRIPTS = [
    'verify_basics.py', 'verify_fourth_moment_bound.py',
    'verify_discrete_radial_remainder.py', 'verify_pair_bound_definition.py',
    'verify_pair_bound_local_preconditioned.py', 'check_arb_bernstein.py',
    'check_physical_boundary_bernstein.py', 'verify_pair_bound_cover.py',
]
TERMINAL = ('outside', 'moment_outside', 'quartic_outside', 'local',
            'positive', 'cubic_positive', 'arb_positive', 'discrete_radial_outside',
            'physical_correction_positive', 'quadratic_correction_positive',
            'arb_correction_positive')

def need(ok, label):
    if not ok:
        raise ArithmeticError(label)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def read(path):
    return json.loads(path.read_text(encoding='utf-8'))

def inputs():
    meta = (HERE / 'UPSTREAM.json').read_bytes()
    need(digest(meta) == UPSTREAM_META_SHA256, 'fixed attributed upstream metadata')
    upstream = json.loads(meta)
    need(upstream['commit'] == UPSTREAM_COMMIT, 'fixed upstream revision')
    files = {}
    for row in read(HERE / 'source_manifest.json')['files']:
        name = row['path']; p = PurePosixPath(name)
        need(not p.is_absolute() and '..' not in p.parts and ':' not in name
             and '\\' not in name and p.as_posix() == name and name not in files,
             'safe unique fixed input')
        raw = (HERE / name).read_bytes()
        need(len(raw) == row['bytes'] and digest(raw) == row['sha256'], 'fixed input ' + name)
        files[name] = raw
    for row in upstream['files']:
        raw = files[row['path']]
        need(digest(raw) == row['sha256'] and len(raw) == row['bytes'], 'unmodified upstream ' + row['path'])
        blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
        need(blob == row['git_blob_sha'], 'upstream Git blob ' + row['path'])
    need(digest(files['vendor/n9/global_pair_bound_exact_contacts.json']) == COEFFICIENT_SHA256,
         'fixed exact kernel coefficients')
    return files

def validate(stage):
    data = stage / 'vendor/n9'
    cover = read(data / 'global_pair_bound_bernstein_frontier.json')
    need(cover['global_majorant_proved'] is True and cover['frontier'] == []
         and cover['stop_reason'] == 'complete' and cover['strict_auxiliary_counterexample'] is None,
         'fresh full continuous-domain cover is complete')
    need(cover['coefficient_sha256'] == COEFFICIENT_SHA256, 'fresh cover binds the exact kernel')
    need(Fraction(cover['U']) == Fraction(1860281, 1048576)
         and cover['minimum_distance_squared'] == '1/4', 'the full necessary domain')
    counts = cover['counts']
    need(all(type(x) is int and x >= 0 for x in counts.values()), 'exact nonnegative tree counts')
    leaves = sum(counts[k] for k in TERMINAL)
    need(counts['visited'] == 2 * counts['split'] + 1 and leaves == counts['split'] + 1,
         'complete binary partition conservation')
    definition = read(data / 'global_pair_bound_exact_definition_certificate.json')
    need(definition['coefficient_sha256'] == COEFFICIENT_SHA256
         and all(definition[k] is True for k in ['kernel_definition_exact',
                  'all_serialized_entries_checked', 'contact_values_and_gradients_exact'])
         and definition['matrix_count'] == 8 and definition['matrix_dimension'] == 6,
         'fresh exact kernel and contact identities')
    local = read(data / 'global_pair_bound_local_preconditioned_certificate.json')
    need(local['coefficient_sha256'] == COEFFICIENT_SHA256
         and local['all_matrices_positive_definite'] is True and local['local_positivity_proved'] is True
         and [r['contact'] for r in local['records']] == list(range(5)),
         'fresh strict positivity and all interior/boundary contacts')
    need(read(data / 'basic_exact_certificate.json')['candidate_stationary_exact'] is True,
         'fresh candidate algebra')
    moment = read(data / 'fourth_moment_radius_certificate.json')
    need(moment['all_positive_local_maxima_bound_proved'] is True
         and Fraction(moment['upper_squared_radius_x']) == Fraction(1860281,1048576),
         'fresh all-maxima moment bound')
    remainder = read(data / 'discrete_radial_remainder_certificate.json')
    need(remainder['piecewise_polynomial_identities'] is True
         and remainder['seven_variable_cauchy_identity'] is True, 'fresh remaining-radii identities')
    return counts, leaves

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-inputs', action='store_true', help='hash binding only; does not replay the proof')
    parser.add_argument('--record-report', action='store_true', help='save fresh results and log only after success')
    parser.add_argument('--seconds', type=int, default=1800, help='full-cover time budget; unfinished coverage fails')
    args = parser.parse_args(); files = inputs()
    if args.check_inputs:
        need(not args.record_report, 'hash-only mode cannot record proof success')
        print(json.dumps({'fixed_inputs_checked':len(files), 'upstream_commit':UPSTREAM_COMMIT,
                          'proof_replay_performed':False}))
        return
    need(args.seconds > 0, 'positive cover budget')
    env = os.environ.copy(); env['PYTHONUTF8'] = '1'
    env.pop('PYTHONOPTIMIZE', None); env.pop('PYTHONPATH', None)
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='mordell-n9-global-') as temp:
        stage = Path(temp)
        for name, raw in files.items():
            p = stage / name; p.parent.mkdir(parents=True, exist_ok=True); p.write_bytes(raw)
        folder = stage / 'vendor/n9'
        # Saved upstream reports are provenance, never proof inputs. Require every
        # mathematical result below to be generated by this execution.
        fresh_names = ['basic_exact_certificate.json', 'fourth_moment_radius_certificate.json',
                       'discrete_radial_remainder_certificate.json',
                       'global_pair_bound_exact_definition_certificate.json',
                       'global_pair_bound_local_preconditioned_certificate.json',
                       'global_pair_bound_bernstein_frontier.json']
        for name in fresh_names:
            (folder / name).unlink(missing_ok=True)
        log = stage / 'fresh_replay.log'
        with log.open('w', encoding='utf-8') as output:
            for name in SCRIPTS:
                print('REPLAY ' + name + '; assertions enabled', flush=True)
                output.write('\nREPLAY ' + name + '\n'); output.flush()
                command = [sys.executable, '-B', str(folder / name)]
                if name == 'verify_pair_bound_cover.py':
                    command.append(str(args.seconds)) # Deliberately no --resume.
                result = subprocess.run(command, cwd=folder, env=env, stdout=output,
                                        stderr=subprocess.STDOUT, timeout=args.seconds + 120)
                need(result.returncode == 0, 'upstream replay failed: ' + log.read_text(encoding='utf-8')[-5000:])
        counts, leaves = validate(stage)
        # The parent may use -S (the repository's standard-library checker does).
        # Query installed distributions in a normal child that loads this venv.
        versions = subprocess.run([sys.executable, '-B', '-c',
            "import importlib.metadata as m,json; print(json.dumps({n:m.version(n) for n in ['sympy','python-flint','numpy','mpmath']}))"],
            cwd=stage, env=env, capture_output=True, text=True, timeout=30)
        need(versions.returncode == 0, 'numerical dependency versions: ' + versions.stderr[-2000:])
        report = {
            'passed':True, 'proof_replay_performed':True, 'full_nine_point_exclusion_proved':True,
            'upstream_commit':UPSTREAM_COMMIT, 'coefficient_sha256':COEFFICIENT_SHA256,
            'fixed_inputs_checked':len(files), 'upstream_numeric_scripts_replayed':len(SCRIPTS),
            'scripts_replayed':SCRIPTS, 'old_zero_frontier_reused':False,
            'child_assertions_enabled':True, 'unresolved_boxes':0,
            'counts':counts, 'terminal_boxes':leaves, 'maximum_mass8':'8^8',
            'maximum_mass9':'9^36/8^28', 'equality':'origin plus regular octagon; rotation and relabeling',
            'nine_point_stationary_classification_used':False,
            'external_peer_review_claimed':False, 'full_formalization':False,
            'versions':{'python':platform.python_version(), **json.loads(versions.stdout)},
            'elapsed_seconds':round(time.monotonic()-start,3),
            'input_sha256':{name:digest(raw) for name,raw in files.items()},
        }
        if args.record_report:
            out = HERE / 'reports'; out.mkdir(exist_ok=True)
            for name in fresh_names:
                shutil.copyfile(folder / name, out / name)
            shutil.copyfile(log, out / 'fresh_replay.log')
            (out / 'replay.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
        print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
