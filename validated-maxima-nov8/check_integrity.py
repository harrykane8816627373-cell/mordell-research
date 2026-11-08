"""Reject false Nov8 root/Hessian claims and changed pinned bytes; -O safe."""
from pathlib import Path
from copy import deepcopy
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from verify_all import load_prior, need, verify_data, DEPENDENCIES

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def run_checks():
    common = json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    prior = load_prior()
    stable = common['rows'][0]
    cases = []
    row = deepcopy(stable)
    row['inverse_numerators'] = [[0 for q in line] for line in row['inverse_numerators']]
    cases.append(('zero Newton preconditioner',row))
    row = deepcopy(stable)
    row['whitener_numerators'] = [[0 for q in line] for line in row['whitener_numerators']]
    cases.append(('singular alleged full Hessian congruence',row))
    row = deepcopy(stable)
    row['center_numerators'][0] += 1 << (common['center_exponent']-3)
    cases.append(('midpoint shifted outside the certified root box',row))
    row = deepcopy(stable)
    row['classification'] = 'unstable_stationary_orbit'
    row['positive_vector_numerators'] = [0 for q in row['center_numerators']]
    cases.append(('false zero positive second-variation direction',row))
    rejected = []
    for label,row in cases:
        try:
            prior.verify_row(row,common)
        except ArithmeticError:
            rejected.append(label)
        else:
            raise ArithmeticError('false certificate accepted: '+label)
    tampered = deepcopy(common)
    tampered['rows'][0]['discovery_integer_points'][0][0] += 1
    try:
        verify_data(tampered,prior,arithmetic_checks=False)
    except ArithmeticError as error:
        need('exact pinned integer seed' in str(error),'integer seed mismatch rejected at provenance check')
        rejected.append('changed discovery seed coordinates')
    else:
        raise ArithmeticError('changed seed accepted')
    with tempfile.TemporaryDirectory(prefix='mordell-nov8-tamper-') as temporary:
        root = Path(temporary)
        stage = root/'validated-maxima-nov8'
        shutil.copytree(HERE,stage)
        altered = stage/'certificates.json'
        original = altered.read_bytes()
        altered.write_bytes(original+b'\n')
        proc = subprocess.run([sys.executable,'-B','-S','-O',str(stage/'verify_all.py')],
                              cwd=root,capture_output=True,text=True,encoding='utf-8',timeout=30)
        need(proc.returncode != 0 and 'fixed input certificates.json' in proc.stderr,'one altered certificate byte rejected')
        rejected.append('one altered byte in the fixed certificate data')
        altered.write_bytes(original)
        for name in DEPENDENCIES:
            path = root/name
            path.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/name,path)
        source = root/'validated-maxima/dyadic.py'
        source.write_bytes(source.read_bytes()+b'\n')
        proc = subprocess.run([sys.executable,'-B','-S','-O',str(stage/'verify_all.py')],
                              cwd=root,capture_output=True,text=True,encoding='utf-8',timeout=30)
        need(proc.returncode != 0 and 'pinned dependency validated-maxima/dyadic.py' in proc.stderr,
             'changed prior interval arithmetic dependency rejected')
        rejected.append('one altered byte in the pinned prior dyadic dependency')
    return {'passed':True,'date':'2026-11-08','optimization_safe':True,
            'false_claims_or_byte_changes_rejected':rejected,
            'prior_files_modified':False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args = parser.parse_args()
    report = run_checks()
    target = HERE/'integrity_report.json'
    if args.record_report:
        target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    else:
        need(report == json.loads(target.read_text(encoding='utf-8')),'fresh false-certificate tests reproduce the frozen report')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
