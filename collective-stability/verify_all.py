"""Replay all four exact entries in a temporary copy of pinned inputs."""
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
    paths=[r['path'] for r in manifest['files']]
    require(len(paths)==len(set(paths)),'unique fixed input paths')
    for row in manifest['files']:
        rel=Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe relative path')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],
                'fixed byte identity '+row['path'])
    entries=[('verify_ring_population.py',[],'ring_population_verification.json'),
             ('verify_radial_quantization.py',[],'radial_quantization_verification.json'),
             ('verify_local_energy.py',['--write-report'],'local_energy_verification.json'),
             ('verify_block_stability.py',[],'block_stability_verification.json')]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='mordell-collective-stability-',dir=HERE) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for entry,args,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',entry,*args],
                cwd=stage/'collective-stability',capture_output=True,text=True,encoding='utf-8',timeout=120)
            require(proc.returncode==0,entry+' failed: '+proc.stdout[-1200:]+proc.stderr[-1800:])
            got=json.loads((stage/'collective-stability'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected,'deterministic report '+report)
            require(got.get('passed') is True or got.get('status')=='PASS','PASS '+entry)
            runs.append({'entry':entry,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+entry,flush=True)
    radial=json.loads((HERE/'radial_quantization_verification.json').read_text(encoding='utf-8'))
    require(F(radial['certificate']['headline_coefficient'])==F(13,30),'new asymptotic headline')
    require(F(radial['positive_asymptotic_margin'])>0,'strict asymptotic margin')
    ring=json.loads((HERE/'ring_population_verification.json').read_text(encoding='utf-8'))
    rows={r['n']:r for r in ring['rows']}
    require(rows[1000]['minimum_layers']==11 and rows[10000]['minimum_layers']==32,'finite population examples')
    block=json.loads((HERE/'block_stability_verification.json').read_text(encoding='utf-8'))
    require(block['three_point_example']['new_whole_box']['excludes_entire_box'],'three-point whole box')
    require(block['six_point_example']['new_whole_box']['excludes_entire_box'],'six-point whole box')
    local=json.loads((HERE/'local_energy_verification.json').read_text(encoding='utf-8'))
    require(F(local['rational_shape_local_plus_j0'])==F(local['rational_shape_angular_delta'])+
            F(local['rational_shape_global_j0_bonus']),'orthogonal local and global j0 add')
    result={'passed':True,'date':'2026-10-29','manifest_inputs':len(paths),'runs':runs,
        'all_n_layer_bound':'K > n/(pi sqrt(n-1)) > sqrt(n)/pi, n>=2',
        'asymptotic_layer_bound':'liminf K_n/sqrt(n) >= 13/30',
        'asymptotic_margin':radial['positive_asymptotic_margin'],
        'finite_population_examples':{'1000':11,'10000':32},
        'independent_integer_compositions':ring['independent_ordered_compositions'],
        'whole_coordinate_boxes_excluded':2,
        'scope':'Necessary layer, population and stability constraints and conditional local energy bounds. No general optimal configuration, full eight-point sharp bound, exact gamma or crystallization.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
