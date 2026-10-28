"""Replay this update with Python's standard library and exact arithmetic.

All generated files are written in a temporary directory. Historical long-running
six/seven point source checks are not rerun or represented as new checks here.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,shutil,subprocess,sys,tempfile,time,zipfile

ROOT=Path(__file__).resolve().parent

def require(ok,label):
    if not ok:raise ArithmeticError(label)

def main():
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    started=time.monotonic()
    manifest=json.loads((ROOT/'source_manifest.json').read_text(encoding='utf-8'))
    for item in manifest['files']:
        b=(ROOT/item['path']).read_bytes()
        require(len(b)==item['bytes'],'size '+item['path'])
        require(hashlib.sha256(b).hexdigest()==item['sha256'],'sha256 '+item['path'])
    updates=json.loads((ROOT/'status_update_manifest.json').read_text(encoding='utf-8'))
    with zipfile.ZipFile(ROOT/'status_updates.zip') as z:
        require(z.testzip() is None,'status ZIP CRC')
        require(sorted(z.namelist())==sorted(x['path']for x in updates['files']),'status ZIP members')
        for item in updates['files']:
            b=z.read(item['path'])
            require(len(b)==item['bytes'] and hashlib.sha256(b).hexdigest()==item['sha256'],'status snapshot '+item['path'])
    entries=[('verify_rational_ascent.py',[],'rational_ascent_checks.json'),
             ('verify_phase_graph.py',['--baseline','.'],'phase_graph_comparison.json'),
             ('verify_radial_distribution.py',[],'radial_distribution_verification.json'),
             ('ascent_phase_independent_check.py',[],'ascent_phase_independent_checks.json'),
             ('compare_ascent_bounds.py',[],'ascent_envelope_comparison.json'),
             ('verify_hermite.py',[],None)]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='mordell-bounds-') as td:
        stage=Path(td)
        for p in ROOT.iterdir():
            if p.suffix in ('.py','.json') and p.name!='verify_all.py':shutil.copyfile(p,stage/p.name)
        for name,args,report in entries:
            t=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name,*args],cwd=stage,capture_output=True,text=True,encoding='utf-8',timeout=180)
            require(proc.returncode==0,name+' failed: '+proc.stderr[-2500:])
            if report:
                got=json.loads((stage/report).read_text(encoding='utf-8'))
                expected=json.loads((ROOT/report).read_text(encoding='utf-8'))
                require(got==expected,'replay mismatch '+report)
                require(got.get('passed') is True,'pass field '+report)
            else:
                got=json.loads((stage/'hermite_result.json').read_text(encoding='utf-8'))
                require(got['status']=='CERTIFIED' and got['corrected_gap_lower']=='3/100','Hermite certificate')
            runs.append({'entry':name,'passed':True,'seconds':round(time.monotonic()-t,3)})
            print('PASS '+name,flush=True)
        # Rebuild the complete old baseline rather than trusting a stored comparison.
        code="import json;from compare_envelopes import audit;from fractions import Fraction as F;from pathlib import Path;a=audit(30);old=json.loads(Path('old_envelope_n2_30.json').read_text());assert all(F(x['final_exact'])==F(y['value']) for x,y in zip(a['rows'],old['rows']));assert len(a['rows'])==len(old['rows'])==29;print('29 complete restored baseline values reproduced')"
        proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-c',code],cwd=stage,capture_output=True,text=True,encoding='utf-8',timeout=120)
        require(proc.returncode==0,'complete baseline replay '+proc.stderr[-2000:])
        print('PASS complete restored baseline n=2..30',flush=True)
    phase=json.loads((ROOT/'phase_graph_comparison.json').read_text(encoding='utf-8'))
    ascent=json.loads((ROOT/'ascent_envelope_comparison.json').read_text(encoding='utf-8'))
    p=phase['strict_complete_envelope_improvements'];a=ascent['ascent_improves_complete_old_plus_new_phase_at']
    require(len(p)==13 and len(a)==13 and not set(p)&set(a),'disjoint finite improvements')
    row26=next(x for x in ascent['rows']if x['n']==26)
    require(row26['ascent_gain_basis_points_floor_over_old_and_phase']>=9830,'n26 certified gain')
    require(Fraction(90155,84693)**12>2,'n30 certified gain')
    report={'passed':True,'date':'2026-10-28','manifest_inputs':len(manifest['files']),
        'status_snapshot_files':len(updates['files']),'runs':runs,'complete_baseline_values_reproduced':29,
        'phase_improved_n':p,'ascent_further_improved_n':a,'combined_distinct_improved_n':sorted(set(p)|set(a)),
        'known_limits':'No unknown finite-n global sharpness, no global convergence from stationarity, no new unconditional upper bound. Regional upper bounds require their explicit radial hypotheses.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (ROOT/'verification_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
