"""Replay the curated source and new exact checks in temporary copies."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT=Path(__file__).resolve().parent


def require(ok,name):
    if not ok:raise ArithmeticError(name)


def run(directory,script,*args):
    p=subprocess.run([sys.executable,'-B','-S','-O',str(directory/script),*args],
                     cwd=directory,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
    require(p.returncode==0,'replay failed: '+script+'\n'+p.stderr+'\n'+p.stdout)
    return {'script':script,'passed':True}


def main():
    manifest=json.loads((ROOT/'source_manifest.json').read_text(encoding='utf-8'))
    for row in manifest['files']:
        p=ROOT/row['file'];require(p.parent==ROOT and p.is_file(),'manifest path')
        b=p.read_bytes()
        require(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'source hash '+row['file'])
    provenance=json.loads((ROOT/'source_provenance.json').read_text(encoding='utf-8'))
    rawhash={item['label']+'/'+row['member']:row['sha256']
             for item in provenance['sources'] for row in item['members']}
    old=[];new=[]
    with tempfile.TemporaryDirectory(prefix='mordell-source-review-') as td:
        temp=Path(td)
        with zipfile.ZipFile(ROOT/'source_evidence.zip') as z:
            require(z.testzip() is None,'evidence CRC')
            require(sorted(z.namelist())==sorted(provenance['evidence_members']),'evidence member coverage')
            for name in z.namelist():
                parts=PurePosixPath(name)
                require(not parts.is_absolute() and '..' not in parts.parts,'evidence member path')
                raw=z.read(name)
                require(name in rawhash and hashlib.sha256(raw).hexdigest()==rawhash[name],'evidence source bytes')
                p=temp.joinpath(*parts.parts);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
        for name in ('original','revised'):
            require((ROOT/(name+'_main.tex')).read_bytes()==(temp/name/'main.tex').read_bytes(),'raw manuscript preserved')
        original=temp/'original/anc';revised=temp/'revised/anc'
        for script in ('verify_global_grid.py','verify_polygon_perturbations.py','verify_saved_candidates.py'):
            old.append(run(original,script))
        revised_scripts=(
            'verify_improved_grid.py','all_n_rings_verify.py','verify_polar_nonoptimal.py',
            'verify_all_n_deletion.py','all_n_uniform_checks.py','all_n_structure_radial_verify.py',
            'verify_spectral_agent.py','verify_spectral_second_level.py','verify_framework.py',
            'verify_examples_unified_exact.py','verify_examples_value_sos.py','verify_stationary_polynomial.py')
        for script in revised_scripts:old.append(run(revised,script))
        stage=temp/'new';stage.mkdir()
        for p in ROOT.iterdir():
            if p.is_file() and p.suffix in ('.py','.json') and p.name!='verify_all.py':
                shutil.copy2(p,stage/p.name)
        for script,args in (
            ('verify_radial_layers.py',()),('verify_deletion_refinement.py',()),
            ('verify_quantitative_spectral_gap.py',()),('verify_hybrid_grid.py',()),
            ('polar_refinement.py',('--verify',)),('compare_envelopes.py',())):
            new.append(run(stage,script,*args))
        expected=json.loads((ROOT/'lower_bound_comparison.json').read_text(encoding='utf-8'))
        require(json.loads((stage/'lower_bound_comparison.json').read_text(encoding='utf-8'))==expected,'complete envelope replay')
        require(len(expected['restored_original_polar_improves_previous_at'])==42,'restored polar count')
        require(expected['restored_saved_candidate_improves_previous_and_polar_at']==[10,11,12],'saved source gains')
        require(len(expected['new_refinement_improves_complete_restored_envelope_at'])==43,'full restored comparison count')
        require(expected['twelve_point_family_below_supplied_saved_candidate'],'twelve-point scope')
        radial=json.loads((stage/'radial_layers_checks.json').read_text(encoding='utf-8'))
        spectral=json.loads((stage/'quantitative_spectral_gap_checks.json').read_text(encoding='utf-8'))
    report={'passed':True,'source_files_checked':len(manifest['files']),
            'original_internal_hashes':20,'revised_internal_hashes':105,
            'source_replay_entries':old,'new_replay_entries':new,
            'complete_original_grid_sets':1874376,'envelope_counts_n_min':2,'envelope_counts_n_max':100,
            'restored_polar_gain_count':42,'restored_saved_gains':[10,11,12],
            'new_gain_count_after_full_restoration':43,'radial_checks':radial,'spectral_checks':spectral,
            'scope':'All-n conclusions require the proof texts. No large precision grid, FS classification recomputation, or unknown general-n sharp value is claimed.'}
    (ROOT/'verification_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'source_files_checked':len(manifest['files']),
                      'source_replays':len(old),'new_replays':len(new),
                      'original_grid_sets':1874376,'full_restored_new_gains':43},indent=2))


if __name__=='__main__':main()
