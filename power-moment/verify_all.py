"""Replay the four targeted exact power-coordinate verifiers offline."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile,time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)

def require(ok,label):
    if not ok:raise ArithmeticError(label)

def main():
    started=time.monotonic()
    manifest=json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    files=manifest['files'];paths=[row['path'] for row in files]
    require(len(paths)==len(set(paths)),'unique fixed inputs')
    for row in files:
        p=Path(row['path'])
        require(not p.is_absolute() and not p.drive and '..' not in p.parts,'safe input')
        data=(ROOT/p).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'fixed bytes '+row['path'])
    entries=[('verify_algebra.py','algebra_verification.json'),
             ('verify_structure.py','structure_verification.json'),
             ('verify_isolation.py','isolation_verification.json'),
             ('verify_regions.py','regions_verification.json')]
    env=os.environ.copy();env['PYTHONUTF8']='1';runs=[]
    with tempfile.TemporaryDirectory(prefix='powers-',dir=ROOT) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for script,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',script],
                cwd=stage/'power-moment',capture_output=True,text=True,encoding='utf-8',errors='replace',env=env,timeout=120)
            require(proc.returncode==0,script+': '+proc.stdout[-1000:]+proc.stderr[-1600:])
            got=json.loads((stage/'power-moment'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected and got['passed'] is True,'deterministic report '+report)
            runs.append({'entry':script,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+script,flush=True)
    status=json.loads((HERE/'regions_verification.json').read_text(encoding='utf-8'))
    provenance=json.loads((HERE/'PROVENANCE.json').read_text(encoding='utf-8'))
    require(status['full_eight_point_exclusion_proved'] is False and provenance['full_eight_point_exclusion_proved'] is False,'partial certificates are not a full-domain proof')
    result={'passed':True,'date':'2026-10-31','manifest_inputs':len(files),'targeted_exact_entries':4,
        'runs':runs,'general_n_powers':'finite Newton closure, exact Hankel objective and lacunary/cyclic-symmetry structure',
        'n8_local_power_radius':'1/10, force-one coordinates, reference p7=896',
        'nonzero_p2_p3_whole_boxes_excluded':status['entire_independent_moment_boxes_excluded'],
        'full_eight_point_exclusion_proved':False,
        'remaining_gap':'No exhaustive inclusion of the complete physical stationary stable target-superlevel set in the certified regions.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
