"""Pinned offline replay of stronger overlap and graph/Dirichlet budgets."""
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
    current=[('overlap-budget','verify_overlap_cost.py','overlap_cost_verification.json'),
             ('overlap-budget','verify_parameters.py','parameters_verification.json'),
             ('overlap-budget','verify_graph_budget.py','graph_budget_verification.json'),
             ('overlap-budget','verify_joint_region.py','joint_region_verification.json')]
    child_env=os.environ.copy();child_env['PYTHONUTF8']='1'
    runs=[]
    with tempfile.TemporaryDirectory(prefix='ob-',dir=ROOT) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for folder,name,report in prior+current:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name],
                cwd=stage/folder,capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env,timeout=180)
            require(proc.returncode==0,folder+'/'+name+': '+proc.stdout[-1000:]+proc.stderr[-1600:])
            got=json.loads((stage/folder/report).read_text(encoding='utf-8'))
            expected=json.loads((ROOT/folder/report).read_text(encoding='utf-8'))
            require(got==expected and got['passed'] is True,'deterministic report '+folder+'/'+report)
            runs.append({'entry':folder+'/'+name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+folder+'/'+name,flush=True)
    params=json.loads((HERE/'parameters_verification.json').read_text(encoding='utf-8'))
    require(F(params['fixed_parameter_certificates']['strong37']['linear_coefficient'])==F(37,2**24),'37-fold linear coefficient')
    graphs=json.loads((HERE/'graph_budget_verification.json').read_text(encoding='utf-8'))
    require(graphs['short_edge_capacity_examples'][1]['maximum_short_edges']==136,'global n100 edge count')
    joint=json.loads((HERE/'joint_region_verification.json').read_text(encoding='utf-8'))
    require(joint['conditional_maximum_short_edges']==2 and joint['input_short_edge_count']==4,'joint edge class')
    result={'passed':True,'date':'2026-10-30','manifest_inputs':len(paths),'prior_entries_replayed':8,'new_entries':4,'runs':runs,
        'global_deficit_lower':'log(Un/Mn)>=37/2^24*(n-11sqrt(n))_+, every n>=2',
        'improvement_factor_over_previous_linear_coefficient':37,
        'pair_near_contact_leading':'64/(105pi)*(2-d)^(7/2)',
        'force_cutoff':'8/5','n8_global_short_edge_cap':17,'n100_global_short_edge_cap':136,
        'n8_conditional_edge_cap':2,
        'scope':'Improved finite global bound and necessary joint distribution constraints; no full optimizer classification or unknown sharp value.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
