"""Extend a force proof forest using edge-correlated remainder bounds."""
from pathlib import Path
from collections import Counter
import argparse,gzip,json,time
import audit_six as A
import edge_force as K

E=K.E
CHECKPOINT_SHA256='2f2af13b0164ded1fc5b199874c76f790b0e3c713fa4b432b36d86525d22b30b'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes',type=int,default=20000)
    parser.add_argument('--roots',type=int,nargs='*')
    parser.add_argument('--leaves',type=int,nargs='*',help='extend only these open leaf IDs from the resumed forest')
    parser.add_argument('--narrow-first',action='store_true',help='process the narrowest selected original open boxes first')
    parser.add_argument('--resume',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--phase-remainder',action='store_true',help='retain complex remainder intervals instead of only their moduli')
    parser.add_argument('--hessian-gram',action='store_true',help='phase remainders plus complex two-column full-Hessian Gram cuts')
    parser.add_argument('--subset-budget',action='store_true',help='shared strict-subset mass and centroid budgets with the accepted force kernel')
    parser.add_argument('--subset-fast',action='store_true',help='one shared-budget/force composition pass per node')
    parser.add_argument('--target-radius',action='store_true',help='apply the target-value/M7/Schwarz radius bound 119/80')
    parser.add_argument('--schwarz2-radius',action='store_true',help='apply the centered second-order Schwarz radius bound 231/160')
    parser.add_argument('--ellipse-schur-radius',action='store_true',help='apply the old-ellipse Schur radius bound 23/16')
    parser.add_argument('--translated-subset',action='store_true',help='fixed-centre strict-subset bounds after the 23/16 kernel')
    parser.add_argument('--sweep-only',action='store_true',help='try every selected open leaf once; retain unresolved original boxes without splitting')
    args=parser.parse_args()
    data=json.loads(gzip.decompress(args.resume.read_bytes()))
    A.need(data['parent_sha256']==A.PARENT_SHA256 and data['attachment_sha256']==A.ZIP_SHA256,'same force frontier model')
    data['checkpoint_sha256']=CHECKPOINT_SHA256
    data['accepted_prefix_sha256']='173bf7f6ccf804a687afbb122af5b9130625f455ef2c104e1e3f2dbef390184f'
    import recent_force
    recent_force.descendant(data,recent_force.load_recent()[0])
    A.need(args.nodes>=0,'nonnegative finite search budget')
    A.need(args.roots is None or set(args.roots)<={r for r,i in data['roots']},'requested original root IDs exist')
    A.need(args.leaves is None or bool(args.leaves) and set(args.leaves)<={i for i,b,r in data['pending']},
        'requested leaf IDs must actually be open in the resumed forest')
    selected=set(args.roots) if args.roots else {r for r,i in data['roots']}
    held=[];pending=[]
    for i,b,r in data['pending']:
        chosen=r in selected and (args.leaves is None or i in args.leaves)
        (pending if chosen else held).append((i,[tuple(p) for p in b],r))
    if args.narrow_first:pending.sort(key=lambda row:max(hi-lo for lo,hi in row[1]),reverse=True)
    initial_selected={i for i,b,r in pending};origins={i:i for i,b,r in pending}
    E.setup(40);settings=E.model('global_maximum')
    kernel=K;engine='joint_force_edge'
    if args.phase_remainder or args.hessian_gram:
        import phase_force
        kernel=phase_force;engine='joint_force_phase'
    if args.hessian_gram:
        import gram_force
        kernel=gram_force;engine='joint_force_gram'
    if args.subset_budget:
        import subset_force
        kernel=subset_force;engine='joint_force_subset'
    if args.subset_fast:
        import subset_force
        class OnePass:
            contract=staticmethod(lambda original,model:subset_force.contract(original,model,passes=1))
        kernel=OnePass;engine='joint_force_subset_fast'
    if args.target_radius:
        import target_radius
        kernel=target_radius;engine='joint_force_target_radius'
    if args.schwarz2_radius:
        import schwarz2_radius
        kernel=schwarz2_radius;engine='joint_force_schwarz2_radius'
    if args.ellipse_schur_radius:
        import ellipse_schur_radius
        kernel=ellipse_schur_radius;engine='joint_force_ellipse_schur_radius'
    if args.translated_subset:
        import translated_subset
        kernel=translated_subset;engine='joint_force_translated_subset'
    counts=Counter();start=time.monotonic()
    retained=[]
    for number in range(args.nodes):
        if not pending:break
        index,original,root=pending.pop()
        box,reason,aux=kernel.contract(original,settings)
        if box is None:
            record={'kind':'excluded','engine':engine,'reason':reason};counts[reason]+=1
        else:
            if args.sweep_only:
                retained.append((index,original,root))
                continue
            dim=max(range(13),key=lambda k:box[k][1]-box[k][0])
            lo,hi=box[dim];cut=(lo+hi)//2
            A.need(lo<cut<hi,'exact proper grid partition')
            left,right=list(box),list(box);left[dim],right[dim]=(lo,cut),(cut,hi)
            li,ri=len(data['nodes']),len(data['nodes'])+1
            data['nodes'] += [{'kind':'unknown'},{'kind':'unknown'}]
            origins[li]=origins[ri]=origins[index]
            record={'kind':'split','engine':engine,'dim':dim,'cut':cut,'left':li,'right':ri}
            pending += [(ri,right,root),(li,left,root)];counts['split']+=1
        data['nodes'][index]=record
        if (number+1)%500==0:print(json.dumps({'processed':number+1,'selected_open':len(pending),
            'closed_selected_old_boxes':len(initial_selected-{origins[i] for i,b,r in retained+pending}),
            'counts':dict(counts),'seconds':time.monotonic()-start}),flush=True)
    pending=retained+pending
    data['pending']=held+pending
    data['complete_cover']=not data['pending'];data['full_eight_point_exclusion_proved']=not data['pending']
    args.output.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),compresslevel=9,mtime=0))
    print(json.dumps({'new_processed':sum(counts.values()),'counts':dict(counts),
        'closed_selected_initial_leaf_ids':sorted(initial_selected-{origins[i] for i,b,r in pending}),
        'selected_roots_complete':not pending,'remaining_parent_regions':sorted({r for i,b,r in data['pending']}),
        'open_leaf_boxes':len(data['pending']),'seconds':time.monotonic()-start,
        'full_eight_point_exclusion_proved':not data['pending']},indent=2))


if __name__=='__main__':main()
