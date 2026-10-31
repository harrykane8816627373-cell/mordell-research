"""Pinned exact replay of annular layer, variable-smearing and joint certificates."""
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

if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)

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
    entries=[
        ('stationary-screening', 'verify_imported.py', 'imported_verification.json'),
        ('stationary-screening', 'verify_stationary_moments.py', 'stationary_moments_verification.json'),
        ('stationary-screening', 'verify_screened_holes.py', 'screened_holes_verification.json'),
        ('stationary-screening', 'verify_edge_covering.py', 'edge_covering_verification.json'),
        ('local-defect', 'verify_geometry.py', 'geometry_verification.json'),
        ('local-defect', 'verify_defect_ledger.py', 'defect_ledger_verification.json'),
        ('local-defect', 'verify_radial_hierarchy.py', 'radial_hierarchy_verification.json'),
        ('local-defect', 'verify_joint_budget.py', 'joint_budget_verification.json'),
        ('overlap-budget', 'verify_overlap_cost.py', 'overlap_cost_verification.json'),
        ('overlap-budget', 'verify_parameters.py', 'parameters_verification.json'),
        ('overlap-budget', 'verify_graph_budget.py', 'graph_budget_verification.json'),
        ('overlap-budget', 'verify_joint_region.py', 'joint_region_verification.json'),
        ('layer-circumference', 'verify_pair_integral.py', 'pair_integral_verification.json'),
        ('layer-circumference', 'verify_ring_geometry.py', 'ring_geometry_verification.json'),
        ('layer-circumference', 'verify_circumference.py', 'circumference_verification.json'),
        ('layer-circumference', 'verify_layers.py', 'layers_verification.json'),
        ('joint-variational', 'verify_affine_pair.py', 'affine_pair_verification.json'),
        ('joint-variational', 'verify_projection.py', 'projection_verification.json'),
        ('joint-variational', 'verify_ledger.py', 'ledger_verification.json'),
        ('joint-variational', 'verify_joint_class.py', 'joint_class_verification.json'),
        ('gaussian-angular', 'verify_kernels.py', 'kernels_verification.json'),
        ('gaussian-angular', 'verify_compatibility.py', 'compatibility_verification.json'),
        ('gaussian-angular', 'verify_multiscale.py', 'multiscale_verification.json'),
        ('gaussian-angular', 'verify_joint_budget.py', 'joint_budget_verification.json'),
        ('annular-layer', 'verify_geometry.py', 'geometry_verification.json'),
        ('annular-layer', 'verify_scaled.py', 'scaled_verification.json'),
        ('annular-layer', 'verify_annular.py', 'annular_verification.json'),
        ('annular-layer', 'verify_joint.py', 'joint_verification.json'),
    ]
    child_env=os.environ.copy();child_env['PYTHONUTF8']='1'
    runs=[]
    with tempfile.TemporaryDirectory(prefix='al-',dir=ROOT) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for folder,name,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name],
                cwd=stage/folder,capture_output=True,text=True,encoding='utf-8',errors='replace',env=child_env,timeout=180)
            require(proc.returncode==0,folder+'/'+name+': '+proc.stdout[-1000:]+proc.stderr[-1600:])
            got=json.loads((stage/folder/report).read_text(encoding='utf-8'))
            expected=json.loads((ROOT/folder/report).read_text(encoding='utf-8'))
            require(got==expected and got['passed'] is True,'deterministic report '+folder+'/'+report)
            runs.append({'entry':folder+'/'+name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+folder+'/'+name,flush=True)
    annular=json.loads((HERE/'annular_verification.json').read_text(encoding='utf-8'))
    joint=json.loads((HERE/'joint_verification.json').read_text(encoding='utf-8'))
    require(all(F(row['strict_margin_over_claim'])>0 for row in annular['asymptotic_certificates']),'two annular asymptotic certificates')
    require(all(F(row['strict_margin'])>0 for row in joint['all_scale_gaussian_contrast_bounds']),'new n100 contrast exclusions')
    result={'passed':True,'date':'2026-10-31','manifest_inputs':len(paths),
        'prior_entries_replayed':24,'new_entries':4,'runs':runs,
        'finite_all_n_layer_bound':'K>0.51*sqrt(n), every n>=10^14',
        'elementary_asymptotic_layer_bound':'liminf K/sqrt(n)>0.513',
        'external_input_asymptotic_layer_bound':'liminf K/sqrt(n)>0.556; no finite starting n',
        'n100_four_layer_p2_strict_upper':'47/20',
        'n100_four_layer_specified_hole_p2_strict_upper':'17/10',
        'n100_all_scale_contrasts':['3/8','73/1000','27/100'],
        'scope':'Necessary layer/distribution bounds and conditional n100 restrictions; no general sharp M_n, n100 five-layer theorem or microscopic classification.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
