"""Verify the published linear-asymptotic constants and scalar identities."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile,time
BASE=Path(__file__).resolve().parent

def main():
 if sys.flags.optimize:raise RuntimeError('Run without -O')
 manifest=json.loads((BASE/'source_manifest.json').read_text(encoding='utf8'))
 for row in manifest['files']:
  p=(BASE/row['path']).resolve()
  if not p.is_relative_to(BASE):raise ValueError('Unsafe source path')
  if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:raise ArithmeticError('Source hash mismatch: '+row['path'])
 workbase=(BASE/'.verification-work').resolve();workbase.mkdir(exist_ok=True)
 report={'status':'RUNNING','date':'2026-10-27','jobs':[],'full_eight_point_bound_proved':False,
  'scope':'Exact constants and scalar coefficient identities. Analytic proofs and explicitly cited external theorems require reading README.md.'}
 result=BASE/'verification_report.json'
 def save():result.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
 start=time.monotonic();save()
 try:
  with tempfile.TemporaryDirectory(prefix='linear-',dir=workbase) as d:
   if not Path(d).resolve().is_relative_to(workbase):raise ValueError('Unsafe temporary target')
   work=Path(d)/'proofs';shutil.copytree(BASE,work,ignore=shutil.ignore_patterns('.verification-work','verification_report.json','__pycache__'))
   for name in ['verify_gamma_bounds.py','verify_finite_bound.py']:
    print('RUN '+name,flush=True);tick=time.monotonic();env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);env['PYTHONIOENCODING']='utf8'
    r=subprocess.run([sys.executable,'-B','-S',str(work/name)],cwd=work,env=env,capture_output=True,text=True,encoding='utf8',errors='replace')
    report['jobs'].append({'script':name,'returncode':r.returncode,'status':'PASS' if r.returncode==0 else 'FAIL','seconds':round(time.monotonic()-tick,3)});save()
    if r.returncode:
     print((r.stdout+'\n'+r.stderr)[-5000:],flush=True);raise RuntimeError('Acceptance failed: '+name)
    cert='gamma_bounds_certificate.json' if name=='verify_gamma_bounds.py' else 'finite_bound_certificate.json'
    if json.loads((work/cert).read_text(encoding='utf8'))!=json.loads((BASE/cert).read_text(encoding='utf8')):raise ArithmeticError('Regenerated certificate differs: '+cert)
    report['jobs'][-1]['published_certificate_matches_replay']=True;save()
    print('PASS '+name,flush=True)
  report['status']='PASS';report['seconds']=round(time.monotonic()-start,3);save();print(json.dumps({'status':'PASS','jobs':len(report['jobs']),'seconds':report['seconds']}),flush=True)
 except Exception as e:
  report['status']='FAIL';report['error']=str(e);save();raise
if __name__=='__main__':main()
