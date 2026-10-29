"""Pinned offline replay of population capacity, transport and local energy."""
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
    require(len(paths)==len(set(paths)),'unique manifest paths')
    for row in manifest['files']:
        rel=Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe fixed input path')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],
                'fixed bytes '+row['path'])
    entries=[('verify_angular_capacity.py',[],'angular_capacity_verification.json'),
             ('verify_radial_transport.py',[],'radial_transport_verification.json'),
             ('verify_effective_layers.py',[],'effective_layers_verification.json'),
             ('verify_overlapping_energy.py',['--write-report'],'overlapping_energy_verification.json')]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='mordell-population-transport-',dir=HERE) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for name,args,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name,*args],
                cwd=stage/'population-transport',capture_output=True,text=True,encoding='utf-8',timeout=180)
            require(proc.returncode==0,name+' failed: '+proc.stdout[-1200:]+proc.stderr[-1800:])
            got=json.loads((stage/'population-transport'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected,'deterministic report replay '+report)
            require(got.get('passed') is True or got.get('status')=='PASS','PASS '+name)
            runs.append({'entry':name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+name,flush=True)
    angular=json.loads((HERE/'angular_capacity_verification.json').read_text(encoding='utf-8'))
    require(angular['eight_point_annulus_certificate']['excluded'] is True,'entire eight-point thin band excluded')
    require(F(angular['eight_point_necessary_radial_variance_strict_lower'])==F(39601,10**8),
            'eight-point positive radial variance')
    radial=json.loads((HERE/'radial_transport_verification.json').read_text(encoding='utf-8'))
    rows={r['n']:r for r in radial['effective_examples']}
    counts={10000:46,10**6:487,10**8:4958,10**12:499575}
    require(all(rows[n]['finite_family_best_minimum_layers']==k for n,k in counts.items()),
            'exact finite parameter family necessary layers')
    effective=json.loads((HERE/'effective_layers_verification.json').read_text(encoding='utf-8'))
    require(effective['rows'][0]['minimum_layers']==45,'closed finite corollary')
    overlap=json.loads((HERE/'overlapping_energy_verification.json').read_text(encoding='utf-8'))
    require(F(overlap['neighborhood_delta_lower_outward'][0])>F(1,3),
            'nonempty crossing-interface neighborhood energy certificate')
    result={'passed':True,'date':'2026-10-29','manifest_inputs':len(paths),'runs':runs,
        'all_local_maximum_layer_bound':'K > n/(1+2sqrt(2(n-1)))',
        'global_effective_layer_bound':'K > sqrt(n)/2-(11/20)n**(1/4), all n>=10000',
        'global_asymptotic_layer_bound':'liminf K_n/sqrt(n) >= 1/2',
        'external_linear_asymptotics_required_for_new_layer_theorem':False,
        'finite_family_minimum_layers':counts,
        'eight_point_strict_radial_variance_lower':'39601/100000000',
        'scope':'Necessary geometric/energy constraints and exact finite certificates; no complete general optimizer, eight-point sharp bound or exact gamma.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
