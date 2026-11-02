"""Replay the accepted parent chain, attachment and complete force-extension proof."""
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
        help='only when the exact e846b8e parent proof chain is already accepted and its 36 arithmetic inputs are unchanged')
    group.add_argument('--legacy-inputs-already-replayed',action='store_true',
        help='caller has replayed joint-closure and its prior inputs; still replay both later parent forests')
    parser.add_argument('--record-reports',action='store_true')
    args=parser.parse_args()
    manifest_raw=(HERE/'source_manifest.json').read_bytes()
    manifest=json.loads(manifest_raw)
    inputs={}
    for row in manifest['files']:
        name=row['path'];p=PurePosixPath(name)
        need(not p.is_absolute() and '..' not in p.parts and ':' not in name and '\\' not in name
            and p.as_posix()==name and name not in inputs,'safe unique pinned input')
        raw=(ROOT/name).read_bytes()
        need(len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],'fixed input '+name)
        inputs[name]=raw
    inputs['joint-refinement/source_manifest.json']=manifest_raw
    env=os.environ.copy();env['PYTHONUTF8']='1'
    results={}
    with tempfile.TemporaryDirectory(prefix='mordell-force-proof-',dir=ROOT) as temp:
        stage=Path(temp)
        for name,raw in inputs.items():
            p=stage/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
        def run(script,*arguments):
            result=subprocess.run([sys.executable,'-B','-S','-O',str(stage/'joint-refinement'/script),*arguments],
                cwd=stage,env=env,stdout=subprocess.PIPE,text=True,encoding='utf-8',timeout=3600)
            need(result.returncode==0,script+': '+result.stdout[-3000:])
            return result.stdout
        if args.parent_already_replayed:print('REUSE accepted e846b8e parent; all 36 arithmetic inputs must match',flush=True)
        else:
            arguments=['--legacy-inputs-already-replayed'] if args.legacy_inputs_already_replayed else []
            run('verify_joint.py',*arguments)
            print('PASS complete accepted parent proof chain',flush=True)
        entries=[('audit_six.py',['joint-refinement/inherited/mordell-eight-point-six-closures-2026-11-01.zip'],'force_import.json'),
            ('verify_force_enclosure.py',[],'force_arithmetic.json'),
            ('check_force.py',['joint-refinement/certificates/force_extension.json.gz','--progress'],'force_cover.json')]
        for script,arguments,name in entries:
            result=json.loads(run(script,*arguments))
            need(result['passed'] is True and result['full_eight_point_exclusion_proved'] is False,'honest new partial proof status')
            if not args.record_reports:
                expected=json.loads((stage/'joint-refinement/reports'/name).read_text())
                need(result==expected,'deterministic exact report '+name)
            results[name]=result
            print('PASS '+script,flush=True)
    execution={'passed':True,'parent_forests_replayed_this_invocation':not args.parent_already_replayed,
        'parent_full_chain_replayed_this_invocation':not (args.parent_already_replayed or args.legacy_inputs_already_replayed),
        'inherited_dependency_proof_reused':args.parent_already_replayed or args.legacy_inputs_already_replayed,
        'parent_chain_reused':args.parent_already_replayed,'unchanged_parent_arithmetic_inputs':36,
        'new_exact_entries_replayed':3,'independent_attachment_replayer_uses_assert':False,
        'full_eight_point_exclusion_proved':False}
    if args.record_reports:
        for name,result in results.items():(HERE/'reports'/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (HERE/'reports/force_execution.json').write_text(json.dumps(execution,indent=2)+'\n',encoding='utf-8')
    cover=results['force_cover.json']
    print(json.dumps(dict(execution,total_closed_parent_regions=cover['total_closed_parent_regions'],
        newly_closed_beyond_upload=cover['newly_closed_beyond_upload'],remaining_parent_regions=cover['remaining_parent_regions'],
        open_leaf_boxes=cover['open_leaf_boxes'],combined_tree_nodes=cover['combined_tree_nodes'],
        combined_excluded_leaves=cover['combined_excluded_leaves']),indent=2))


if __name__=='__main__':main()
