"""Offline replay of local geometric energy and matched radial hierarchy."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def main():
    started=time.monotonic()
    manifest=json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    paths=[r['path'] for r in manifest['files']]
    require(len(paths)==len(set(paths)),'unique manifest paths')
    for row in manifest['files']:
        rel=Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe fixed input')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'fixed input '+row['path'])
    entries=[('verify_geometry.py','geometry_verification.json'),
             ('verify_defect_ledger.py','defect_ledger_verification.json'),
             ('verify_radial_hierarchy.py','radial_hierarchy_verification.json'),
             ('verify_joint_budget.py','joint_budget_verification.json')]
    runs=[]
    child_env=os.environ.copy()
    child_env['PYTHONUTF8']='1'
    # Prior acceptance stages its own dependencies once more. A short root
    # stage prevents Windows working-directory limits in that nested replay.
    with tempfile.TemporaryDirectory(prefix='ld-',dir=ROOT) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        prior=subprocess.run([sys.executable,'-B','-S','-O','stationary-screening/verify_all.py'],
            cwd=stage,capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env,timeout=180)
        require(prior.returncode==0,'prior suite: '+prior.stdout[-1000:]+prior.stderr[-1200:])
        prior_report=json.loads((stage/'stationary-screening/verification_report.json').read_text(encoding='utf-8'))
        require(prior_report['passed'] is True and len(prior_report['runs'])==4,'four prior entries replayed')
        print('PASS prior stationary-screening suite (4 entries)',flush=True)
        for name,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name],
                cwd=stage/'local-defect',capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env,timeout=180)
            require(proc.returncode==0,name+': '+proc.stdout[-1000:]+proc.stderr[-1600:])
            got=json.loads((stage/'local-defect'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected and got['passed'] is True,'deterministic report '+report)
            runs.append({'entry':name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+name,flush=True)
    geometry=json.loads((HERE/'geometry_verification.json').read_text(encoding='utf-8'))
    require(F(geometry['constants']['per_point_lower'])==F(1,2**24),'explicit linear coefficient')
    ledger=json.loads((HERE/'defect_ledger_verification.json').read_text(encoding='utf-8'))
    require(F(ledger['n144_new_linear_discount_lower'])==F(3,4194304),'finite comparison')
    joint=json.loads((HERE/'joint_budget_verification.json').read_text(encoding='utf-8'))
    require(F(joint['combined_class_energy_lower'])>F(joint['eight_global_gap_strict_upper']),'whole joint moment class excluded')
    result={'passed':True,'date':'2026-10-30','manifest_inputs':len(paths),
        'prior_stationary_screening_entries_replayed':4,'runs':runs,
        'global_deficit_lower':'log(Un/Mn) >= 2^-24*(n-11sqrt(n))_+, every n>=2',
        'all_n_at_least_484':'log(Un/Mn) >= n/2^25',
        'asymptotic_gap_liminf_lower':'1/16777216',
        'conditional_gamma_upper':'3/4-log(2)/2-2^-24; inherited external expansion required only for this interpretation',
        'eight_joint_excluded_class':'sum |z_i|^4>=52/5 and |sum z_i^2|>=9/10',
        'scope':'Finite global bound and necessary moment/distribution constraints; no full optimal configuration or sharp unknown value.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
