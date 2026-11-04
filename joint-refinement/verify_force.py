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
    group.add_argument('--checkpoint-already-replayed',action='store_true',
        help='reuse exactly the accepted 768598a force prefix and parent chain; replay all descendants and check 50 unchanged arithmetic inputs')
    group.add_argument('--accepted-prefix-already-replayed',action='store_true',
        help='reuse exactly the accepted 8be9f91 force prefix and parent chain; replay all descendants and check 61 unchanged arithmetic inputs')
    group.add_argument('--nov2-prefix-already-replayed',action='store_true',
        help='reuse exactly the accepted 19dbbd3 force prefix and parent chain; replay all descendants and check 70 unchanged arithmetic inputs')
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
                cwd=stage,env=env,stdout=subprocess.PIPE,text=True,encoding='utf-8',
                timeout=10800 if script=='check_force.py' else 3600)
            need(result.returncode==0,script+': '+result.stdout[-3000:])
            return result.stdout
        if args.nov2_prefix_already_replayed:
            print('REUSE accepted 19dbbd3 prefix (94,399 processed nodes) and parent chain; all 70 arithmetic inputs must match',flush=True)
        elif args.accepted_prefix_already_replayed:
            print('REUSE accepted 8be9f91 prefix (60,791 processed nodes) and parent chain; all 61 arithmetic inputs must match',flush=True)
        elif args.checkpoint_already_replayed:
            print('REUSE accepted 768598a prefix (25,790 processed nodes) and parent chain; all 50 arithmetic inputs must match',flush=True)
        elif args.parent_already_replayed:print('REUSE accepted e846b8e parent; all 36 arithmetic inputs must match',flush=True)
        else:
            arguments=['--legacy-inputs-already-replayed'] if args.legacy_inputs_already_replayed else []
            run('verify_joint.py',*arguments)
            print('PASS complete accepted parent proof chain',flush=True)
        entries=[('audit_six.py',['joint-refinement/inherited/mordell-eight-point-six-closures-2026-11-01.zip'],'force_import.json'),
            ('verify_force_enclosure.py',[],'force_arithmetic.json'),
            ('verify_edge_force.py',[],'edge_arithmetic.json'),
            ('verify_subset_budget.py',[],'subset_arithmetic.json'),
            ('verify_target_radius.py',[],'target_radius.json'),
            ('verify_schwarz2_radius.py',[],'schwarz2_radius.json'),
            ('verify_ellipse_schur_radius.py',[],'ellipse_schur_radius.json'),
            ('verify_nov2_prefix.py',[],'nov2_prefix.json'),
            ('verify_translated_subset.py',[],'translated_subset.json'),
            ('check_force.py',['joint-refinement/certificates/force_extension.json.gz','--progress']+
                (['--reuse-nov2-prefix'] if args.nov2_prefix_already_replayed else
                    ['--reuse-accepted-prefix'] if args.accepted_prefix_already_replayed else
                    ['--reuse-checkpoint'] if args.checkpoint_already_replayed else []),'force_cover.json')]
        for script,arguments,name in entries:
            result=json.loads(run(script,*arguments))
            need(result['passed'] is True,'successful exact proof entry')
            if name!='force_cover.json':need(result['full_eight_point_exclusion_proved'] is False,'arithmetic lemma alone is not global exclusion')
            else:need(result['full_eight_point_exclusion_proved'] is (result['open_leaf_boxes']==0),
                'complete exclusion requires zero unresolved coordinate boxes')
            if not args.record_reports:
                expected=json.loads((stage/'joint-refinement/reports'/name).read_text())
                need(result==expected,'deterministic exact report '+name)
            results[name]=result
            print('PASS '+script,flush=True)
    cover=results['force_cover.json']
    reused=args.parent_already_replayed or args.checkpoint_already_replayed or args.accepted_prefix_already_replayed or args.nov2_prefix_already_replayed
    prefix_nodes=94399 if args.nov2_prefix_already_replayed else 60791 if args.accepted_prefix_already_replayed else 25790 if args.checkpoint_already_replayed else 0
    execution={'passed':True,'parent_forests_replayed_this_invocation':not reused,
        'parent_full_chain_replayed_this_invocation':not (reused or args.legacy_inputs_already_replayed),
        'inherited_dependency_proof_reused':reused or args.legacy_inputs_already_replayed,
        'parent_chain_reused':reused,'unchanged_parent_arithmetic_inputs':36,
        'checkpoint_prefix_reused':args.checkpoint_already_replayed,
        'checkpoint_processed_nodes_reused':25790 if args.checkpoint_already_replayed else 0,
        'unchanged_checkpoint_arithmetic_inputs':50 if args.checkpoint_already_replayed else 0,
        'accepted_prefix_reused':args.accepted_prefix_already_replayed,
        'accepted_prefix_processed_nodes_reused':60791 if args.accepted_prefix_already_replayed else 0,
        'unchanged_accepted_prefix_arithmetic_inputs':61 if args.accepted_prefix_already_replayed else 0,
        'nov2_prefix_reused':args.nov2_prefix_already_replayed,
        'nov2_prefix_processed_nodes_reused':94399 if args.nov2_prefix_already_replayed else 0,
        'unchanged_nov2_prefix_arithmetic_inputs':70 if args.nov2_prefix_already_replayed else 0,
        'force_processed_nodes_replayed_this_invocation':cover['processed_nodes']-prefix_nodes,
        'new_exact_entries_replayed':10,'independent_attachment_replayer_uses_assert':False,
        'full_eight_point_exclusion_proved':cover['full_eight_point_exclusion_proved']}
    if args.record_reports:
        for name,result in results.items():(HERE/'reports'/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (HERE/'reports/force_execution.json').write_text(json.dumps(execution,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(execution,total_closed_parent_regions=cover['total_closed_parent_regions'],
        newly_closed_beyond_upload=cover['newly_closed_beyond_upload'],remaining_parent_regions=cover['remaining_parent_regions'],
        open_leaf_boxes=cover['open_leaf_boxes'],combined_tree_nodes=cover['combined_tree_nodes'],
        combined_excluded_leaves=cover['combined_excluded_leaves']),indent=2))


if __name__=='__main__':main()
