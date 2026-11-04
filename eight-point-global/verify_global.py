"""Replay the attributed fixed n8 global certificate in a disposable copy."""
from pathlib import Path,PurePosixPath
from fractions import Fraction
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time

HERE=Path(__file__).resolve().parent
UPSTREAM_COMMIT='31cc934b7dfa0ef2308566a1bce48cbfb2e01559'
UPSTREAM_META_SHA256='a4a765c6ceede97c56ee51d1a233c5dd7e5ca6cbc06b45cd5657729b28210d2a'
COEFFICIENT_SHA256='405ebbccbe70b6da7c69f13e1fdbb130625f5aba98898729cb2cdf16153995d9'

def need(ok,label):
    if not ok:raise ArithmeticError(label)
def digest(raw):return hashlib.sha256(raw).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))

def inputs():
    raw=(HERE/'UPSTREAM.json').read_bytes()
    need(digest(raw)==UPSTREAM_META_SHA256,'fixed attributed upstream metadata')
    upstream=json.loads(raw)
    need(upstream['commit']==UPSTREAM_COMMIT,'fixed upstream revision')
    manifest=read(HERE/'source_manifest.json');files={}
    for row in manifest['files']:
        name=row['path'];p=PurePosixPath(name)
        need(not p.is_absolute() and '..' not in p.parts and ':' not in name and '\\' not in name
            and p.as_posix()==name and name not in files,'safe unique fixed input')
        value=(HERE/name).read_bytes()
        need(len(value)==row['bytes'] and digest(value)==row['sha256'],'fixed input '+name)
        files[name]=value
    for row in upstream['files']:
        need(row['path'] in files and digest(files[row['path']])==row['sha256']
            and len(files[row['path']])==row['bytes'],'unmodified attributed upstream file '+row['path'])
    need(digest(files['vendor/n8/data/global_pair_bound_exact_contacts.json'])==COEFFICIENT_SHA256,
        'fixed exact kernel coefficients')
    return files

def validate(stage,files):
    data=stage/'vendor/n8/data'
    final=read(data/'FINAL_CERTIFICATE.json');cover=read(data/'global_pair_bound_bernstein_frontier.json')
    need(final['global_majorant_proved'] is True and final['all_required_computational_checks_passed'] is True
        and final['unresolved_boxes']==0 and cover['global_majorant_proved'] is True and cover['frontier']==[],
        'fresh full continuous-domain cover is complete')
    need(final['coefficient_sha256']==cover['coefficient_sha256']==COEFFICIENT_SHA256,
        'fresh reports bind the exact coefficients')
    counts=cover['counts']
    need(counts==final['cover_counts'] and all(type(v) is int and v>=0 for v in counts.values()),'fresh exact tree counts')
    need(counts['visited']==2*counts['split']+1
        and counts['outside']+counts['local']+counts['positive']==counts['split']+1,'complete binary partition conservation')
    need(Fraction(cover['U'])==Fraction(965267,524288)
        and cover['minimum_distance_squared']=='2/7','the full verified domain contains the existing subset domain')
    for name,flag in [('global_pair_bound_exact_definition_certificate.json','kernel_definition_exact'),
            ('global_pair_bound_psd_certificate.json','exact_auxiliary_matrices_positive_definite'),
            ('global_pair_bound_local_certificate.json','local_positivity_proved')]:
        result=read(data/name)
        need(result[flag] is True and result['coefficient_sha256']==COEFFICIENT_SHA256,'fresh exact/positivity/contact result '+name)
    need(read(data/'global_pair_bound_exact_definition_certificate.json')['contact_values_and_gradients_exact'] is True,
        'exact contact derivatives')
    for name,sha in final['files_sha256'].items():
        need('vendor/n8/'+name in files and digest(files['vendor/n8/'+name])==sha,
            'fresh result binds its unchanged upstream script '+name)
    need(len(final['files_sha256'])==6,'all six original numerical scripts are bound')
    return final,counts

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-inputs',action='store_true',help='hash binding only; does not execute the proof')
    parser.add_argument('--record-report',action='store_true',help='save the fresh reports and log only after success')
    args=parser.parse_args();files=inputs()
    if args.check_inputs:
        need(not args.record_report,'hash-only mode cannot record proof success')
        print(json.dumps({'fixed_inputs_checked':len(files),'upstream_commit':UPSTREAM_COMMIT,'proof_replay_performed':False}))
        return
    env=os.environ.copy();env['PYTHONUTF8']='1'
    env.pop('PYTHONOPTIMIZE',None);env.pop('PYTHONPATH',None)
    start=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='mordell-n8-global-') as temp:
        stage=Path(temp)
        for name,raw in files.items():
            path=stage/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
        bridge=subprocess.run([sys.executable,'-B','-S','-O',str(stage/'verify_bridge.py')],
            cwd=stage,env=env,capture_output=True,text=True,encoding='utf-8',timeout=60)
        need(bridge.returncode==0,'exact bridge check: '+bridge.stderr[-2000:])
        bridge_result=json.loads(bridge.stdout);need(bridge_result['passed'] is True,'exact bridge identities')
        log=stage/'fresh_replay.log'
        print('REPLAY all six fixed upstream scripts; assertions enabled; start from the full domain',flush=True)
        # The child deliberately receives neither -O nor -S. Its own explicit
        # __debug__ guard rejects disabled assertions, even if this driver uses -O.
        with log.open('w',encoding='utf-8') as output:
            result=subprocess.run([sys.executable,'-B',str(stage/'vendor/n8/verify_all.py')],
                cwd=stage/'vendor/n8',env=env,stdout=output,stderr=subprocess.STDOUT,timeout=3600)
        need(result.returncode==0,'upstream full replay failed: '+log.read_text(encoding='utf-8')[-5000:])
        final,counts=validate(stage,files)
        report={'passed':True,'proof_replay_performed':True,'full_eight_point_exclusion_proved':True,
            'upstream_commit':UPSTREAM_COMMIT,'coefficient_sha256':COEFFICIENT_SHA256,
            'fixed_inputs_checked':len(files),'upstream_numeric_scripts_replayed':6,
            'old_zero_frontier_reused':False,'child_assertions_enabled':True,'unresolved_boxes':0,
            'counts':counts,'terminal_boxes':counts['outside']+counts['local']+counts['positive'],
            'exact_bridge':bridge_result,'maximum_mass7':'7^7','maximum_mass8':'8^28/7^21',
            'equality':'origin plus regular heptagon; rotation and relabeling',
            'nonzero_s2_s3_remaining_system_empty':True,'seven_point_classification_used_by_bridge':False,
            'external_peer_review_claimed':False,'full_formalization':False,
            'versions':final['versions'],'elapsed_seconds':round(time.monotonic()-start,3),
            'input_sha256':{name:digest(raw) for name,raw in files.items()}}
        if args.record_report:
            out=HERE/'reports';out.mkdir(exist_ok=True)
            for path in (stage/'vendor/n8/data').glob('*.json'):
                if path.name!='global_pair_bound_exact_contacts.json':shutil.copyfile(path,out/path.name)
            shutil.copyfile(log,out/'fresh_replay.log')
            (out/'replay.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
        print(json.dumps(report,indent=2))

if __name__=='__main__':main()
