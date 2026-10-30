"""Pinned offline replay of convex pair loss and layer-circumference transport."""
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
    paths=[row['path'] for row in manifest['files']]
    require(len(paths)==len(set(paths)),'unique manifest inputs')
    for row in manifest['files']:
        rel=Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe fixed input')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'fixed input '+row['path'])
    prior=[('stationary-screening','verify_imported.py','imported_verification.json'),
           ('stationary-screening','verify_stationary_moments.py','stationary_moments_verification.json'),
           ('stationary-screening','verify_screened_holes.py','screened_holes_verification.json'),
           ('stationary-screening','verify_edge_covering.py','edge_covering_verification.json'),
           ('local-defect','verify_geometry.py','geometry_verification.json'),
           ('local-defect','verify_defect_ledger.py','defect_ledger_verification.json'),
           ('local-defect','verify_radial_hierarchy.py','radial_hierarchy_verification.json'),
           ('local-defect','verify_joint_budget.py','joint_budget_verification.json')]
    previous_overlap=[('overlap-budget','verify_overlap_cost.py','overlap_cost_verification.json'),
             ('overlap-budget','verify_parameters.py','parameters_verification.json'),
             ('overlap-budget','verify_graph_budget.py','graph_budget_verification.json'),
             ('overlap-budget','verify_joint_region.py','joint_region_verification.json')]
    current=[('layer-circumference','verify_pair_integral.py','pair_integral_verification.json'),
             ('layer-circumference','verify_ring_geometry.py','ring_geometry_verification.json'),
             ('layer-circumference','verify_circumference.py','circumference_verification.json'),
             ('layer-circumference','verify_layers.py','layers_verification.json')]
    child_env=os.environ.copy();child_env['PYTHONUTF8']='1'
    runs=[]
    with tempfile.TemporaryDirectory(prefix='lc-',dir=ROOT) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for folder,name,report in prior+previous_overlap+current:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name],
                cwd=stage/folder,capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env,timeout=180)
            require(proc.returncode==0,folder+'/'+name+': '+proc.stdout[-1000:]+proc.stderr[-1600:])
            got=json.loads((stage/folder/report).read_text(encoding='utf-8'))
            expected=json.loads((ROOT/folder/report).read_text(encoding='utf-8'))
            require(got==expected and got['passed'] is True,'deterministic report '+folder+'/'+report)
            runs.append({'entry':folder+'/'+name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+folder+'/'+name,flush=True)
    accepted=json.loads((HERE/'layers_verification.json').read_text(encoding='utf-8'))
    require(F(accepted['finite_layer_coefficient'])==F(5007,10000),'finite layer coefficient')
    require(accepted['finite_all_n_threshold']==10**16,'explicit finite threshold')
    require(accepted['n8_global_short_edge_cap']==12 and accepted['n100_global_short_edge_cap']==101,'short-edge counts')
    result={'passed':True,'date':'2026-10-30','manifest_inputs':len(paths),
        'prior_entries_replayed':12,'new_entries':4,'runs':runs,
        'finite_layer_bound':'K>(5007/10000)*sqrt(n), every n>=10^16; elementary effective gap only',
        'asymptotic_layer_bound_with_inherited_external_energy_input':'liminf K/sqrt(n)>=1/4+22/(25pi)>53/100',
        'force_cutoff':'8/5','n8_global_short_edge_cap':12,'n100_global_short_edge_cap':101,
        'scope':'Necessary radial-layer and point-distribution constraints; no full optimizer classification. The two layer bounds have different energy inputs.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
