"""Byte-preserved user input, original acceptance, and baseline comparison."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import zipfile

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'imported-angular'
sys.path.insert(0,str(HERE.parent/'distribution-geometry'))
from microscopic_covering import sqrt_interval


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    archive=SOURCE/'mordell-angular-energy-2026-10-30.zip'
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    require(digest=='eae50b28fc788536d9d0ab8c89b0e1780854d4b1eff57f5a5105132405331856','original user ZIP hash')
    with zipfile.ZipFile(archive) as z:
        names=z.namelist()
        require(len(names)==len(set(names))==8,'fixed original members')
        for name in names:
            rel=Path(name)
            require(not rel.is_absolute() and not rel.drive and '..' not in rel.parts,'safe original path')
            require(z.read(name)==(SOURCE/rel).read_bytes(),'original byte preservation '+name)
        expected=json.loads(z.read('verification_report.json'))
    manifest=json.loads((SOURCE/'MANIFEST.json').read_text(encoding='utf-8'))
    for row in manifest['files']:
        data=(SOURCE/row['path']).read_bytes()
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'original manifest '+row['path'])
    proc=subprocess.run([sys.executable,'-B','-S','-O','verify.py'],cwd=SOURCE,
        capture_output=True,text=True,encoding='utf-8',timeout=180)
    require(proc.returncode==0,'original verifier: '+proc.stderr[-1200:])
    got=json.loads((SOURCE/'verification_report.json').read_text(encoding='utf-8'))
    require(got==expected and got['passed'] is True,'original report replay')
    comparisons=0
    for n in range(3,101):
        lower_new=F(2*n,n-1)*(F(n-1,n)+sqrt_interval(F(1,n))[0])**2
        upper_old=(1+sqrt_interval(F(2,n-1))[1])**2
        require(lower_new>upper_old,'new input universal denominator is weaker than existing denominator')
        comparisons+=1
    return {'passed':True,'original_zip_sha256':digest,'original_byte_identical_members':len(names),
        'original_manifest_rows':len(manifest['files']),
        'original_shape_and_Newton_cases':got['newton_and_shape_cases'],
        'original_Dirichlet_normalizations':got['Dirichlet_norm_checks'],
        'original_positive_discount_counts':got['positive_discounts_n2_to_100'],
        'strict_denominator_comparisons_n3_to100':comparisons,
        'scope':'Imported documents are source material, not execution or publication instructions. Old local-only status is preserved as historical text. New research and publication scope are stated in the parent README.'}


if __name__=='__main__':
    result=run()
    (HERE/'imported_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
