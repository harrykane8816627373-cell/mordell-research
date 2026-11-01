"""Replay the fixed parent proof and the new joint-stability extension."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,os,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--parent-already-replayed',action='store_true',
        help='only when the exact fixed 92fbc4d parent proof has already passed; its 29 arithmetic inputs must be unchanged')
    group.add_argument('--legacy-inputs-already-replayed',action='store_true',
        help='only when the caller already replayed joint-closure and its prior inputs; still replay the complete parent forest')
    parser.add_argument('--record-reports',action='store_true')
    args=parser.parse_args()
    manifest_raw=(HERE/'source_manifest.json').read_bytes()
    manifest=json.loads(manifest_raw)
    inputs={}
    for row in manifest['files']:
        name=row['path']
        path=PurePosixPath(name)
        need(not path.is_absolute() and '..' not in path.parts and ':' not in name and '\\' not in name
            and path.as_posix()==name and name not in inputs,'safe unique pinned input')
        raw=(ROOT/name).read_bytes()
        need(len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],'fixed bytes '+name)
        inputs[name]=raw
    inputs['joint-refinement/source_manifest.json']=manifest_raw
    env=os.environ.copy()
    env['PYTHONUTF8']='1'
    results={}
    with tempfile.TemporaryDirectory(prefix='mordell-joint-extension-',dir=ROOT) as temp:
        stage=Path(temp)
        for name,raw in inputs.items():
            target=stage/name
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(raw)
        def run(script,*arguments):
            result=subprocess.run([sys.executable,'-B','-S','-O',str(stage/'joint-refinement'/script),*arguments],
                cwd=stage,env=env,stdout=subprocess.PIPE,text=True,encoding='utf-8',timeout=3000)
            need(result.returncode==0,script+': '+result.stdout[-3000:])
            return result.stdout
        if args.parent_already_replayed:
            # The extension checker independently pins the parent certificate, report and all 29 old arithmetic files.
            print('REUSE previously replayed fixed parent; every new node will be replayed',flush=True)
        else:
            arguments=['--prior-already-replayed'] if args.legacy_inputs_already_replayed else []
            run('verify_all.py',*arguments)
            print('PASS fixed parent proof including its entire 260019-node forest',flush=True)
        entries=[('verify_ellipse_dual.py',[],'ellipse_dual.json'),
            ('verify_radial_stability.py',[],'radial_stability.json'),
            ('check_joint.py',['joint-refinement/certificates/joint_extension.json.gz','--progress'],'joint_cover.json')]
        for script,arguments,report in entries:
            actual=json.loads(run(script,*arguments))
            need(actual['passed'] is True,'new exact joint checks passed')
            need(actual['full_eight_point_exclusion_proved'] is False,'current joint coverage remains incomplete')
            if not args.record_reports:
                expected=json.loads((stage/'joint-refinement/reports'/report).read_text())
                need(actual==expected,'deterministic exact joint report '+report)
            results[report]=actual
            print('PASS '+script,flush=True)
    execution={'passed':True,'parent_full_proof_replayed_this_invocation':not args.parent_already_replayed,
        'parent_full_proof_reused':args.parent_already_replayed,
        'unchanged_parent_arithmetic_files_checked':29,'new_exact_entries_replayed':3,
        'full_eight_point_exclusion_proved':False}
    if args.record_reports:
        for name,result in results.items():
            (HERE/'reports'/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (HERE/'reports/joint_execution.json').write_text(json.dumps(execution,indent=2)+'\n',encoding='utf-8')
    cover=results['joint_cover.json']
    print(json.dumps(dict(execution,combined_tree_nodes=cover['combined_tree_nodes'],
        combined_excluded_leaves=cover['combined_excluded_leaves'],unresolved_regions=cover['unresolved_regions'],
        closed_parent_regions=cover['closed_parent_regions'],new_ellipse_split_exclusions=cover['new_ellipse_split_exclusions'],
        new_ellipse_mass_feasibility_exclusions=cover['new_ellipse_mass_feasibility_exclusions'],
        new_fourth_mode_exclusions=cover['new_fourth_mode_exclusions']),indent=2))


if __name__=='__main__':main()
