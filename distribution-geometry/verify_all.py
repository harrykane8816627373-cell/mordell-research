"""Pinned standard-library replay of four planar distribution methods."""
from fractions import Fraction as F
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
    paths=[row['path'] for row in manifest['files']]
    require(len(paths)==len(set(paths)),'unique manifest paths')
    for row in manifest['files']:
        rel=Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe input path')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],
                'fixed bytes '+row['path'])
    entries=[('verify_microscopic_covering.py',[],'microscopic_covering_verification.json'),
             ('verify_distribution_metric.py',['--write-report'],'distribution_metric_verification.json'),
             ('verify_local_population.py',[],'local_population_verification.json'),
             ('verify_subset_geometry.py',[],'subset_geometry_verification.json')]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='mordell-distribution-',dir=HERE) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for name,args,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name,*args],
                cwd=stage/'distribution-geometry',capture_output=True,text=True,encoding='utf-8',timeout=180)
            require(proc.returncode==0,name+' failed: '+proc.stdout[-1200:]+proc.stderr[-1800:])
            got=json.loads((stage/'distribution-geometry'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected,'deterministic report replay '+report)
            require(got.get('passed') is True or got.get('status')=='PASS','PASS '+name)
            runs.append({'entry':name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+name,flush=True)
    micro=json.loads((HERE/'microscopic_covering_verification.json').read_text(encoding='utf-8'))
    require(micro['population_examples'][0]['lower']==246,'microscopic population lower')
    require(micro['population_examples'][0]['upper']==10200,'microscopic population upper')
    transport=json.loads((HERE/'distribution_metric_verification.json').read_text(encoding='utf-8'))
    require(F(transport['asymptotic_upper_outward'])<F(981,1000),'conditional narrow asymptotic constant')
    require(F(transport['prior_witness_W1_upper_outward'])<F('0.111199952674'),'finite n=100 W1 bound')
    require(transport['wider_cross_strip_excluded_with_stored_witness'] is True,'nonradial class exclusion')
    local=json.loads((HERE/'local_population_verification.json').read_text(encoding='utf-8'))
    require(F(local['two_hole_certificate']['delta_lower'])>F(16705,12),'two disjoint holes excluded')
    strip=json.loads((HERE/'subset_geometry_verification.json').read_text(encoding='utf-8'))
    require(strip['old_full_pair_slacks_checked']==19900,'full pair witness count')
    require(strip['nonempty_full_witness']['point_count']==200,'full strip-class witness')
    result={'passed':True,'date':'2026-10-30','manifest_inputs':len(paths),'runs':runs,
        'bulk_covering':'For every global optimizer: min_i|z_i-u| < 6/sqrt(n-1), all |u|<=1 and n>=64.',
        'effective_planar_distribution_rate':'W1(eta_z,mu0)=O(n^-1/2), effective with the elementary prior gap.',
        'n100_W1_upper_outward':transport['prior_witness_W1_upper_outward'],
        'asymptotic_W1_upper_outward_external_inputs_retained':transport['asymptotic_upper_outward'],
        'all_center_bounded_radius_disk_discrepancy_rate':'O(n^-1/3)',
        'full_background_covering_rate':'O(n^-1/4), Coulomb x scale',
        'thin_strip_example':'n=200: any length 27/10 width 1/100 rectangle has fewer than 20 local-maximum points.',
        'scope':'Necessary distribution and geometry results; no complete optimizer, unknown sharp value or crystallization.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
