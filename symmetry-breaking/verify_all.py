"""Offline exact replay from pinned bytes; no optimizer or optional packages."""
from pathlib import Path
import hashlib
import json
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
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe relative input')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'],'byte count '+row['path'])
        require(hashlib.sha256(data).hexdigest()==row['sha256'],'sha256 '+row['path'])
    entries=[('verify_stationarity.py','stationarity_verification.json'),
             ('verify_pair_stability.py','pair_stability_verification.json'),
             ('verify_witnesses.py','witness_verification.json')]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='mordell-symmetry-breaking-') as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for name,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name],
                cwd=stage/'symmetry-breaking',capture_output=True,text=True,encoding='utf-8',timeout=180)
            require(proc.returncode==0,name+' failed: '+proc.stdout[-1000:]+proc.stderr[-1800:])
            got=json.loads((stage/'symmetry-breaking'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected,'deterministic replay '+report)
            require(got['passed'] is True,'PASS '+name)
            runs.append({'entry':name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+name,flush=True)
    stationarity=json.loads((HERE/'stationarity_verification.json').read_text(encoding='utf-8'))
    witness=json.loads((HERE/'witness_verification.json').read_text(encoding='utf-8'))
    require(witness['strict_complete_envelope_improvements']==list(range(10,101)),
            '91 strictly improved complete prior lower bounds')
    require(stationarity['nonstationary_selected_count']==91,'91 selected nonstationary candidates')
    require(len(stationarity['entire_continuous_phase_families_excluded_by_population'])==23,
            '23 entire population families excluded')
    rows={r['n']:r for r in witness['rows']}
    require(rows[27]['relative_gain_basis_points_floor']>=10412,'n27 gain at least 104.12 percent')
    require(rows[100]['relative_gain_basis_points_floor']>=3512,'n100 gain at least 35.12 percent')
    require(all(rows[n]['new_minimum_layers']>=3 for n in range(39,43)),'39..42 at least 3 radii')
    require(all(rows[n]['new_minimum_layers']>=4 for n in range(80,86)),'80..85 at least 4 radii')
    result={'passed':True,'date':'2026-10-29','manifest_inputs':len(paths),'runs':runs,
        'strict_complete_envelope_improvements':91,'improved_range':[10,100],
        'selected_nonstationary_ring_candidates':91,'entire_population_families_excluded':23,
        'new_global_necessary_layers':{'39..42':3,'80..85':4},
        'n27_gain_basis_points_floor':rows[27]['relative_gain_basis_points_floor'],
        'n100_gain_basis_points_floor':rows[100]['relative_gain_basis_points_floor'],
        'n100_global_deficit_upper':rows[100]['global_optimizer_gap_upper'],
        'scope':'General necessary constraints and constructive lower bounds. No general sharp optimizer, full eight-point proof, exact gamma, or crystallization theorem.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
