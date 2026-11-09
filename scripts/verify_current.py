"""Check navigation and small suites; optionally replay the large integer identity."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent.parent


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first3-archive', type=Path)
    parser.add_argument('--joint', action='store_true', help='also replay the partial joint-cover tree')
    parser.add_argument('--refinement', action='store_true', help='also replay the current imported/force closure proof chain')
    parser.add_argument('--global-eight', action='store_true', help='fresh full attributed n8 global proof replay; needs the pinned numerical dependencies')
    parser.add_argument('--global-nine', action='store_true', help='fresh full attributed n9 global proof replay; needs the pinned numerical dependencies')
    parser.add_argument('--jacobian', action='store_true', help='also run native Jacobian checks and regenerate the latest kernel gradient')
    parser.add_argument('--jacobian-python', type=Path, help='optional Jacobian interpreter; implies --jacobian')
    args = parser.parse_args()
    pages = [ROOT/'README.md', ROOT/'configuration-bounds/PROOF_STATUS.md',
             *sorted((ROOT/'docs').glob('*.md')), *sorted((ROOT/'eight-point').glob('*.md')),
             *sorted((ROOT/'algebraic-gap').glob('*.md')), *sorted((ROOT/'radial-gap').glob('*.md')),
             *sorted((ROOT/'joint-closure').glob('*.md')), *sorted((ROOT/'joint-refinement').glob('*.md')),
             *sorted((ROOT/'eight-point-global').glob('*.md')),
             *sorted((ROOT/'nine-point-global').glob('*.md')),
             *sorted((ROOT/'general-structure').glob('*.md')),
             *sorted((ROOT/'ring-rigidity').glob('*.md')),
             *sorted((ROOT/'phase-packing').glob('*.md')),
             *sorted((ROOT/'finite-search-nov5').glob('*.md')),
             *sorted((ROOT/'validated-maxima').glob('*.md')),
             *sorted((ROOT/'four-ring-rigidity').glob('*.md')),
             *sorted((ROOT/'phase-deficit-sharp').glob('*.md')),
             *sorted((ROOT/'family-exclusion-nov8').glob('*.md')),
             *sorted((ROOT/'validated-maxima-nov8').glob('*.md')),
             *sorted((ROOT/'candidate-global-screen').glob('*.md')),
             *sorted((ROOT/'maximizer-domain').glob('*.md')),
             *sorted((ROOT/'candidate-basin').glob('*.md')),
             *sorted((ROOT/'candidate-kernel').glob('*.md')),
             *sorted((ROOT/'twelve-point-kernel-upper').glob('*.md')),
             *sorted((ROOT/'twelve-point-kernel-refinement').glob('*.md')),
             *sorted((ROOT/'jacobian-assist').glob('*.md')),
             *sorted((ROOT/'jacobian-gradient').glob('*.md')),
             *sorted((ROOT/'twelve-point-kernel-refinement-2').glob('*.md')),
             *sorted((ROOT/'twelve-point-subset-domain').glob('*.md')),
             *sorted((ROOT/'stationary-kernel').glob('*.md')),
             *sorted((ROOT/'twelve-point-ward-upper').glob('*.md'))]
    links = 0
    for page in pages:
        text = page.read_text(encoding='utf-8')
        for target in re.findall(r'\[[^\]\n]+\]\(([^)\n]+)\)', text):
            parsed = urlsplit(target.strip('<>'))
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            resolved = (page.parent/unquote(parsed.path)).resolve()
            require(resolved.is_relative_to(ROOT) and resolved.exists(),
                    'missing or escaping local link: '+str(page.relative_to(ROOT))+' -> '+target)
            links += 1
    print('PASS maintained navigation: '+str(len(pages))+' pages, '+str(links)+' local links', flush=True)
    inputs = {}
    # radial-gap includes the two inherited eight-point entries; do not run them twice.
    suites = ('power-moment', 'radial-gap')
    manifests = {}
    for suite in suites:
        manifest_path = ROOT/suite/'source_manifest.json'
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        manifests[suite+'/source_manifest.json'] = manifest_path.read_bytes()
        for row in manifest['files']:
            name = row['path']
            resolved = (ROOT/name).resolve()
            require(resolved.is_relative_to(ROOT), 'manifest input in repository')
            data = resolved.read_bytes()
            require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                    'fixed bytes '+name)
            require(name not in inputs or inputs[name] == data, 'consistent shared dependency')
            inputs[name] = data
    inputs.update(manifests)
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-current-', dir=ROOT) as temp:
        stage = Path(temp)
        for name, data in inputs.items():
            dest = stage/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
        for suite in suites:
            proc = subprocess.run([sys.executable, '-X', 'int_max_str_digits=0', '-B', '-S', '-O',
                str(stage/suite/'verify_all.py')], cwd=stage, env=env,
                capture_output=True, text=True, encoding='utf-8', timeout=180)
            require(proc.returncode == 0, suite+': '+proc.stdout[-2000:]+proc.stderr[-3000:])
            print('PASS '+suite+'/verify_all.py', flush=True)
    proc = subprocess.run([sys.executable, '-B', '-S', '-O',
        str(ROOT/'scripts/verify_ten_point_transition.py')], cwd=ROOT, env=env,
        capture_output=True, text=True, encoding='utf-8', timeout=30)
    require(proc.returncode == 0 and json.loads(proc.stdout)['passed'] is True,
            'exact n10 transition: '+proc.stderr[-2000:])
    print('PASS n10 strict integer witness; no global optimum claimed.', flush=True)
    with tempfile.TemporaryDirectory(prefix='mordell-structure-') as temp:
        stage = Path(temp)
        for source in (ROOT/'general-structure').iterdir():
            if source.is_file() and source.suffix in ('.py', '.md'):
                (stage/source.name).write_bytes(source.read_bytes())
        proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/'verify_all.py')],
            cwd=stage, env=env, capture_output=True, text=True, encoding='utf-8', timeout=180)
        require(proc.returncode == 0 and json.loads(proc.stdout)['passed'] is True,
                'continuous phase structure: '+proc.stdout[-2000:]+proc.stderr[-3000:])
    print('PASS continuous-phase compatibility, strict family gaps and independent exact product checks.', flush=True)
    for suite in ('ring-rigidity', 'phase-packing', 'finite-search-nov5', 'validated-maxima'):
        with tempfile.TemporaryDirectory(prefix='mordell-nov5-') as temp:
            stage = Path(temp)
            shutil.copytree(ROOT/suite, stage/suite, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            if suite == 'finite-search-nov5':
                baseline = stage/'symmetry-breaking'
                baseline.mkdir()
                shutil.copyfile(ROOT/'symmetry-breaking/integer_witnesses.json', baseline/'integer_witnesses.json')
            proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/suite/'verify_all.py')],
                cwd=stage, env=env, capture_output=True, text=True, encoding='utf-8', timeout=180)
            require(proc.returncode == 0 and json.loads(proc.stdout)['passed'] is True,
                    suite+': '+proc.stdout[-2000:]+proc.stderr[-3000:])
        print('PASS '+suite+'/verify_all.py (standard-library exact acceptance).', flush=True)
    nov8_dependencies = {
        'four-ring-rigidity': (),
        'phase-deficit-sharp': ('phase-packing/certificates',),
        'family-exclusion-nov8': ('phase-deficit-sharp', 'symmetry-breaking/integer_witnesses.json'),
        'validated-maxima-nov8': ('validated-maxima/verify_all.py', 'validated-maxima/dyadic.py',
            'validated-maxima/PROOF.md', 'validated-maxima/discover.py',
            'finite-search-nov5/integer_witnesses.json', 'symmetry-breaking/integer_witnesses.json'),
    }
    for suite, dependencies in nov8_dependencies.items():
        with tempfile.TemporaryDirectory(prefix='mordell-nov8-') as temp:
            stage = Path(temp)
            for name in (suite, *dependencies):
                source, dest = ROOT/name, stage/name
                if source.is_dir():
                    shutil.copytree(source, dest, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
                else:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, dest)
            proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/suite/'verify_all.py')],
                cwd=stage, env=env, capture_output=True, text=True, encoding='utf-8', timeout=240)
            require(proc.returncode == 0, suite+': '+proc.stdout[-2000:]+proc.stderr[-3000:])
            # The nonlinear verifier appends a short rejection-count diagnostic.
            result, _ = json.JSONDecoder().raw_decode(proc.stdout.lstrip())
            require(result['passed'] is True, suite+': exact acceptance flag')
        print('PASS '+suite+'/verify_all.py (fixed dependencies and exact Nov8 acceptance).', flush=True)
    candidate_dependencies = {
        'candidate-global-screen': ('validated-maxima/verify_all.py', 'validated-maxima/dyadic.py',
            'validated-maxima/certificates.json', 'validated-maxima/PROOF.md', 'symmetry-breaking/explore.py',
            'symmetry-breaking/integer_witnesses.json', 'finite-search-nov5/integer_witnesses.json',
            'radial-angular/population_search_verification.json'),
        'candidate-basin': ('validated-maxima/verify_all.py', 'validated-maxima/dyadic.py',
            'validated-maxima/certificates.json', 'validated-maxima/PROOF.md'),
        'maximizer-domain': ('eight-point-global', 'nine-point-global', 'finite-n/SUBSET_PROOF.md',
            'symmetry-breaking/integer_witnesses.json'),
        'candidate-kernel': (),
        'twelve-point-kernel-upper': ('candidate-kernel/data/n12_exact_interpretation.json',
            'candidate-kernel/verify_psd.py', 'candidate-kernel/verify_all.py',
            'validated-maxima/verify_all.py', 'validated-maxima/dyadic.py',
            'validated-maxima/certificates.json', 'validated-maxima/PROOF.md',
            'symmetry-breaking/integer_witnesses.json', 'finite-n/SUBSET_PROOF.md'),
        'twelve-point-kernel-refinement': ('twelve-point-kernel-upper/verify_cover.py',
            'twelve-point-kernel-upper/verify_all.py', 'candidate-kernel/verify_all.py',
            'candidate-kernel/verify_psd.py', 'validated-maxima/verify_all.py',
            'validated-maxima/dyadic.py', 'validated-maxima/certificates.json',
            'validated-maxima/PROOF.md', 'symmetry-breaking/integer_witnesses.json',
            'finite-n/SUBSET_PROOF.md'),
        'jacobian-assist': ('candidate-kernel/verify_psd.py',
            'candidate-kernel/data/n12_exact_interpretation.json',
            'twelve-point-kernel-refinement/kernel.json',
            'candidate-kernel/discover_sdp.py', 'validated-maxima/certificates.json'),
        'jacobian-gradient': ('twelve-point-kernel-refinement/kernel.json',
            'twelve-point-kernel-upper/verify_cover.py', 'candidate-kernel/discover_sdp.py'),
        'twelve-point-subset-domain': ('nine-point-global', 'finite-n/SUBSET_PROOF.md',
            'symmetry-breaking/integer_witnesses.json'),
        'stationary-kernel': ('ring-rigidity/THEOREMS.md', 'source-review/revised_main.tex'),
    }
    refined_kernel = 'twelve-point-kernel-refinement-2'
    latest_kernel = 'twelve-point-ward-upper'
    candidate_dependencies[refined_kernel] = candidate_dependencies['twelve-point-kernel-refinement']
    candidate_dependencies[latest_kernel] = (*candidate_dependencies[refined_kernel], 'stationary-kernel/README.md')
    for suite, dependencies in candidate_dependencies.items():
        with tempfile.TemporaryDirectory(prefix='mordell-candidate-') as temp:
            stage = Path(temp)
            for name in (suite, *dependencies):
                source, dest = ROOT/name, stage/name
                if source.is_dir():
                    shutil.copytree(source, dest, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
                else:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, dest)
            timeout = 900 if suite in (refined_kernel, latest_kernel) else 240
            proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/suite/'verify_all.py')],
                cwd=stage, env=env, capture_output=True, text=True, encoding='utf-8', timeout=timeout)
            require(proc.returncode == 0, suite+': '+proc.stdout[-2000:]+proc.stderr[-3000:])
            result, _ = json.JSONDecoder().raw_decode(proc.stdout.lstrip())
            require(result['passed'] is True, suite+': acceptance flag')
        if suite == 'twelve-point-kernel-upper':
            print('PASS fresh full n12 continuous kernel cover and global bound V12 < M12 < (23/20)*V12; no sharp value claimed.', flush=True)
        elif suite == 'twelve-point-kernel-refinement':
            print('PASS fresh full refined n12 cover: 87385 nodes, no pending boxes, M12 < (103/100)*V12; no sharp value claimed.', flush=True)
        elif suite == 'jacobian-assist':
            print('PASS fixed actual Jacobian records and independent rational checks; backend not restarted in this step.', flush=True)
        elif suite == refined_kernel:
            print('PASS refined n12 global upper bound M12 < (301/300)*V12 after full 229113-node replay.', flush=True)
        elif suite == latest_kernel:
            print('PASS Ward n12 global upper bound M12 < (379/378)*V12 after full 268239-node replay.', flush=True)
        elif suite in ('jacobian-gradient', 'twelve-point-subset-domain', 'stationary-kernel'):
            print('PASS '+suite+'/verify_all.py with standard-library exact arithmetic.', flush=True)
        else:
            print('PASS '+suite+'/verify_all.py (candidate structure, domain, local caps or rejected kernel fits; no global sharp bound).', flush=True)
    if args.jacobian or args.jacobian_python:
        interpreter = args.jacobian_python or Path.home()/'.codex/tools/jacobian/0.23.0/venv/Scripts/python.exe'
        require(interpreter.is_file(), 'Jacobian interpreter missing; see jacobian-assist/README.md or supply --jacobian-python')
        proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(ROOT/'jacobian-assist/verify_all.py'),
            '--jacobian-python', str(interpreter)], cwd=ROOT, env=env,
            capture_output=True, text=True, encoding='utf-8', timeout=240)
        require(proc.returncode == 0, 'actual Jacobian re-execution: '+proc.stdout[-2000:]+proc.stderr[-3000:])
        require(json.loads(proc.stdout)['Jacobian_executed_now'] is True, 'actual backend and MCP execution flag')
        print('PASS fresh optional Jacobian native computations and real stdio MCP operation.', flush=True)
        proc = subprocess.run([str(interpreter), '-B', str(ROOT/'scripts/build_ward_gradient.py')], cwd=ROOT, env=env,
            capture_output=True, text=True, encoding='utf-8', timeout=180)
        require(proc.returncode == 0, 'latest native gradient generation: '+proc.stderr[-3000:])
        gradient = json.loads(proc.stdout)
        require(gradient['passed'] is True and gradient['includes_ward_corrections'] is True and
                gradient['global_coefficient_identity_checked'] is True,
                'fresh complete native gradient coefficients')
        print('PASS latest kernel gradient from '+str(gradient['actual_Jacobian_native_gradient_calls'])+
              ' native Jacobian calls.', flush=True)
    if args.joint or args.refinement:
        proc = subprocess.run([sys.executable, '-B', '-S', '-O',
            str(ROOT/'joint-closure/verify_all.py'), '--skip-prior'],
            cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=900)
        require(proc.returncode == 0, 'joint-closure: '+proc.stdout[-2000:]+proc.stderr[-3000:])
        print('PASS joint-closure replay; the stored cover still has 21 unresolved regions.', flush=True)
    else:
        print('NOT RUN joint-closure cover replay: supply --joint (the stored cover remains incomplete).')
    if args.refinement:
        proc = subprocess.run([sys.executable, '-B', '-S', '-O',
            str(ROOT/'joint-refinement/verify_force.py'), '--legacy-inputs-already-replayed'],
            cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=18000)
        require(proc.returncode == 0, 'joint-refinement: '+proc.stdout[-3000:]+proc.stderr[-3000:])
        print('PASS imported/force closure continuation; the stored cover remains incomplete.', flush=True)
    else:
        print('NOT RUN joint-refinement continuation: supply --refinement.')
    if args.first3_archive is not None:
        proc = subprocess.run([sys.executable, '-B', '-S', '-O',
            str(ROOT/'algebraic-gap/verify_all.py'), '--source-archive', str(args.first3_archive.resolve())],
            cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=240)
        require(proc.returncode == 0, 'algebraic-gap: '+proc.stdout[-2000:]+proc.stderr[-3000:])
        print('PASS algebraic-gap/verify_all.py: full inherited identity and quantitative gap', flush=True)
    else:
        print('NOT RUN algebraic-gap full replay: supply --first3-archive (see docs/REPRODUCING.md).')
    global_dir=ROOT/'eight-point-global'
    proc=subprocess.run([sys.executable,'-B','-S','-O',str(global_dir/'verify_global.py'),'--check-inputs'],
        cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',timeout=60)
    require(proc.returncode==0,'global n8 fixed inputs: '+proc.stderr[-3000:])
    accepted=json.loads((global_dir/'accepted_replay.json').read_text(encoding='utf-8'))
    manifest=(global_dir/'source_manifest.json').read_bytes()
    require(hashlib.sha256(manifest).hexdigest()==accepted['source_manifest_sha256'],'accepted global n8 input set unchanged')
    for name,expected in accepted['files_sha256'].items():
        require(hashlib.sha256((global_dir/name).read_bytes()).hexdigest()==expected,'accepted global n8 evidence '+name)
    result=json.loads((global_dir/'reports/replay.json').read_text(encoding='utf-8'))
    require(result['passed'] is True and result['proof_replay_performed'] is True
        and result['full_eight_point_exclusion_proved'] is True and result['unresolved_boxes']==0,
        'accepted complete n8 global replay')
    bindings={row['path']:row['sha256'] for row in json.loads(manifest)['files']}
    require(result['input_sha256']==bindings,'accepted n8 result binds the current inputs')
    proc=subprocess.run([sys.executable,'-B','-S','-O',str(global_dir/'verify_bridge.py')],
        cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',timeout=60)
    require(proc.returncode==0,'global n8 exact bridge: '+proc.stderr[-3000:])
    if args.global_eight:
        proc=subprocess.run([sys.executable,'-B','-S','-O',str(global_dir/'verify_global.py')],
            cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',timeout=4000)
        require(proc.returncode==0,'fresh n8 global replay: '+proc.stdout[-3000:]+proc.stderr[-3000:])
        print('PASS fresh full attributed n8 global replay and exact bridge.',flush=True)
    else:
        print('PASS fixed n8 inputs, accepted replay evidence and exact bridge; full cover reused (use --global-eight to replay).',flush=True)
    global_nine = ROOT/'nine-point-global'
    proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(global_nine/'verify_global.py'), '--check-inputs'],
        cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=60)
    require(proc.returncode == 0, 'global n9 fixed inputs: '+proc.stderr[-3000:])
    accepted = json.loads((global_nine/'accepted_replay.json').read_text(encoding='utf-8'))
    manifest = (global_nine/'source_manifest.json').read_bytes()
    require(hashlib.sha256(manifest).hexdigest() == accepted['source_manifest_sha256'],
            'accepted global n9 input set unchanged')
    for name, expected in accepted['files_sha256'].items():
        require(hashlib.sha256((global_nine/name).read_bytes()).hexdigest() == expected,
                'accepted global n9 evidence '+name)
    result = json.loads((global_nine/'reports/replay.json').read_text(encoding='utf-8'))
    require(result['passed'] is True and result['proof_replay_performed'] is True
            and result['full_nine_point_exclusion_proved'] is True and result['unresolved_boxes'] == 0
            and result['old_zero_frontier_reused'] is False and result['child_assertions_enabled'] is True,
            'accepted complete n9 global replay')
    bindings = {row['path']:row['sha256'] for row in json.loads(manifest)['files']}
    require(result['input_sha256'] == bindings, 'accepted n9 result binds the current inputs')
    counts = result['counts']
    require(counts['visited'] == 2*counts['split']+1 and result['terminal_boxes'] == counts['split']+1,
            'accepted n9 binary partition conservation')
    if args.global_nine:
        proc = subprocess.run([sys.executable, '-B', '-S', '-O', str(global_nine/'verify_global.py')],
            cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=2100)
        require(proc.returncode == 0, 'fresh n9 global replay: '+proc.stdout[-3000:]+proc.stderr[-3000:])
        print('PASS fresh full attributed n9 global replay.', flush=True)
    else:
        print('PASS fixed n9 inputs and accepted fresh full replay evidence; use --global-nine to replay.', flush=True)
    print('PASS requested checks. The n8 and n9 global bounds have accepted full replays; general-n sharp values remain open.')


if __name__ == '__main__':
    main()
