"""Pinned offline replay of imported input and three new combinations."""
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
    require(len(paths)==len(set(paths)),'unique pinned paths')
    for row in manifest['files']:
        rel=Path(row['path'])
        require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe pinned path')
        data=(ROOT/rel).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],
                'fixed input '+row['path'])
    entries=[('verify_imported.py','imported_verification.json'),
             ('verify_stationary_moments.py','stationary_moments_verification.json'),
             ('verify_screened_holes.py','screened_holes_verification.json'),
             ('verify_edge_covering.py','edge_covering_verification.json')]
    runs=[]
    with tempfile.TemporaryDirectory(prefix='mordell-stationary-screening-',dir=HERE) as tmp:
        stage=Path(tmp)
        for rel in paths:
            target=stage/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,target)
        for name,report in entries:
            before=time.monotonic()
            proc=subprocess.run([sys.executable,'-X','int_max_str_digits=0','-B','-S','-O',name],
                cwd=stage/'stationary-screening',capture_output=True,text=True,encoding='utf-8',timeout=180)
            require(proc.returncode==0,name+' failed: '+proc.stdout[-1000:]+proc.stderr[-1600:])
            got=json.loads((stage/'stationary-screening'/report).read_text(encoding='utf-8'))
            expected=json.loads((HERE/report).read_text(encoding='utf-8'))
            require(got==expected and got['passed'] is True,'deterministic PASS report '+report)
            runs.append({'entry':name,'passed':True,'seconds':round(time.monotonic()-before,3)})
            print('PASS '+name,flush=True)
    moments=json.loads((HERE/'stationary_moments_verification.json').read_text(encoding='utf-8'))
    require(F(moments['eight_point_certificate']['p2_modulus_strict_upper'])==F(8,7),'new eight-point p2 bound')
    require(moments['eliminated_higher_moment_checks']==21,'quartic lower-moment interface')
    holes=json.loads((HERE/'screened_holes_verification.json').read_text(encoding='utf-8'))
    require(F(holes['n100_second_moment_screened_energy'])<F(holes['n100_certified_global_gap_upper'])
            <F(holes['n100_stationary_dipole_screened_energy']),'actual n100 gap crosses new conditional projection')
    cover=json.loads((HERE/'edge_covering_verification.json').read_text(encoding='utf-8'))
    row=next(r for r in cover['whole_background_examples'] if r['n']==10**6)
    require(F(row['whole_background_covering_upper_outward'])<F(1165,100000),'million-point whole background cover')
    result={'passed':True,'date':'2026-10-30','manifest_inputs':len(paths),'runs':runs,
        'eight_point_global_moment_region':'|P2(z)|^2+(27/100)|P3(z)|^2<64/49; |P2(z)|<8/7',
        'new_nonlinear_interface':'Eliminating modes 4 and 5 leaves positive |P2|^4 and |P2*P3|^2 terms.',
        'screened_holes':'Full Schur conditioning on matched second moment and stationary compact modes.',
        'whole_background_covering_rate':'O(sqrt(log(n)/n)), original z scale, every global optimizer, n>=64',
        'whole_background_million_point_radius_upper':row['whole_background_covering_upper_outward'],
        'scope':'Necessary moment, local energy and distribution theorems. Complete eight-point sharp bound, general optimizer and crystallization remain open.',
        'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'verification_report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
