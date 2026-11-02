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
    parser.add_argument('--resume',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--phase-remainder',action='store_true',help='retain complex remainder intervals instead of only their moduli')
    parser.add_argument('--hessian-gram',action='store_true',help='phase remainders plus complex two-column full-Hessian Gram cuts')
    parser.add_argument('--sweep-only',action='store_true',help='try every selected open leaf once; retain unresolved original boxes without splitting')
    args=parser.parse_args()
    data=json.loads(gzip.decompress(args.resume.read_bytes()))
    A.need(data['parent_sha256']==A.PARENT_SHA256 and data['attachment_sha256']==A.ZIP_SHA256,'same force frontier model')
    data['checkpoint_sha256']=CHECKPOINT_SHA256
    selected=set(args.roots) if args.roots else {r for r,i in data['roots']}
    held=[];pending=[]
    for i,b,r in data['pending']:
        (pending if r in selected else held).append((i,[tuple(p) for p in b],r))
    E.setup(40);settings=E.model('global_maximum')
    kernel=K;engine='joint_force_edge'
    if args.phase_remainder or args.hessian_gram:
        import phase_force
        kernel=phase_force;engine='joint_force_phase'
    if args.hessian_gram:
        import gram_force
        kernel=gram_force;engine='joint_force_gram'
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
            record={'kind':'split','engine':engine,'dim':dim,'cut':cut,'left':li,'right':ri}
            pending += [(ri,right,root),(li,left,root)];counts['split']+=1
        data['nodes'][index]=record
        if (number+1)%500==0:print(json.dumps({'processed':number+1,'selected_open':len(pending),
            'counts':dict(counts),'seconds':time.monotonic()-start}),flush=True)
    pending=retained+pending
    data['pending']=held+pending
    data['complete_cover']=not data['pending'];data['full_eight_point_exclusion_proved']=not data['pending']
    args.output.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),compresslevel=9,mtime=0))
    print(json.dumps({'new_processed':sum(counts.values()),'counts':dict(counts),
        'selected_roots_complete':not pending,'remaining_parent_regions':sorted({r for i,b,r in data['pending']}),
        'open_leaf_boxes':len(data['pending']),'seconds':time.monotonic()-start,
        'full_eight_point_exclusion_proved':not data['pending']},indent=2))


if __name__=='__main__':main()
