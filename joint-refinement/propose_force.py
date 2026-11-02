"""Compose all uploaded trees and extend selected unresolved parent regions."""
from pathlib import Path
from collections import deque,Counter
import argparse,gzip,hashlib,json,time
import audit_six as A
import force_contract as K

HERE=Path(__file__).resolve().parent
ATTACHMENT=HERE/'inherited/mordell-eight-point-six-closures-2026-11-01.zip'
E=K.E


def initialize(path=ATTACHMENT):
    archive,frontier,frontier_sha=A.input_data(path)
    nodes=[{'kind':'unknown'} for _ in frontier]
    roots=[[root,i] for i,(root,b) in enumerate(frontier)]
    pending={root:[(i,[tuple(p) for p in b],root)] for i,(root,b) in enumerate(frontier)}
    names=[(root,'certificates/closed_'+str(root)+'.json.gz') for root in A.CLOSED]
    names += [(root,'certificates/open_'+str(root)+'_'+str(n)+'.json.gz') for root,n in A.PARTIAL]
    for root,name in names:
        data=json.loads(gzip.decompress(archive.read(name)))
        queue=deque(pending[root])
        A.need(data['original_roots']==[root] and data['input_file_sha256']==frontier_sha,'bound uploaded root')
        for row in data['proof_nodes']:
            i,box,parent=queue.popleft()
            A.need(row[1]==root and row[2]==[list(v) for v in box],'exact uploaded prefix box')
            if row[0]=='X':record={'kind':'excluded','engine':'uploaded','reason':row[3]}
            elif row[0]=='S':
                dim,cut=row[3:]
                A.need(type(dim) is int and 0<=dim<13 and box[dim][0]<cut<box[dim][1],'proper uploaded split')
                left,right=list(box),list(box)
                left[dim],right[dim]=(box[dim][0],cut),(cut,box[dim][1])
                li,ri=len(nodes),len(nodes)+1
                nodes += [{'kind':'unknown'},{'kind':'unknown'}]
                record={'kind':'split','engine':'uploaded','dim':dim,'cut':cut,'left':li,'right':ri}
                queue.extend([(li,left,root),(ri,right,root)])
            else:raise ArithmeticError('known uploaded step')
            nodes[i]=record
        A.need([(r,[list(p) for p in b]) for i,b,r in queue]==[(r,b) for r,b in data['pending']],
            'uploaded unresolved leaves preserved verbatim')
        pending[root]=list(queue)
    ordered=[entry for root,b in frontier for entry in pending[root]]
    return {'parent_sha256':A.PARENT_SHA256,'attachment_sha256':A.ZIP_SHA256,'precision_bits':40,
        'scope':'global_maximum','kernel_version':1,'roots':roots,'nodes':nodes,'pending':ordered,
        'complete_cover':not ordered,'full_eight_point_exclusion_proved':not ordered}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes',type=int,default=20000)
    parser.add_argument('--roots',nargs='*',type=int)
    parser.add_argument('--resume',type=Path)
    parser.add_argument('--global-endpoint',action='store_true',help='use the inherited endpoint distance bound in force remainders')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data=json.loads(gzip.decompress(args.resume.read_bytes())) if args.resume else initialize()
    A.need(data['parent_sha256']==A.PARENT_SHA256 and data['attachment_sha256']==A.ZIP_SHA256
        and data['kernel_version']==1,'same proposed force proof model')
    E.setup(40)
    settings=E.model('global_maximum')
    selected=set(args.roots) if args.roots else {i for i,j in data['roots']}
    held=[]
    pending=[]
    for i,b,r in data['pending']:
        row=(i,[tuple(p) for p in b],r)
        (pending if r in selected else held).append(row)
    counts=Counter()
    start=time.monotonic()
    for number in range(args.nodes):
        if not pending:break
        i,original,root=pending.pop()
        box,reason,aux=K.contract(original,settings,endpoint_separation=settings['separation'] if args.global_endpoint else None)
        engine='joint_force_global' if args.global_endpoint else 'joint_force'
        if box is None:
            record={'kind':'excluded','engine':engine,'reason':reason}
            counts[reason]+=1
        else:
            dim=max(range(13),key=lambda k:box[k][1]-box[k][0])
            lo,hi=box[dim]
            cut=(lo+hi)//2
            if not lo<cut<hi:raise ArithmeticError('grid resolution exhausted')
            left,right=list(box),list(box)
            left[dim],right[dim]=(lo,cut),(cut,hi)
            li,ri=len(data['nodes']),len(data['nodes'])+1
            data['nodes'] += [{'kind':'unknown'},{'kind':'unknown'}]
            record={'kind':'split','engine':engine,'dim':dim,'cut':cut,'left':li,'right':ri}
            pending += [(ri,right,root),(li,left,root)]
            counts['split']+=1
        data['nodes'][i]=record
        if (number+1)%500==0:print(json.dumps({'processed':number+1,'selected_open':len(pending),
            'counts':dict(counts),'seconds':time.monotonic()-start}),flush=True)
    data['pending']=held+pending
    data['complete_cover']=not data['pending']
    data['full_eight_point_exclusion_proved']=not data['pending']
    args.output.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),compresslevel=9,mtime=0))
    print(json.dumps({'new_processed':sum(counts.values()),'counts':dict(counts),'remaining_leaf_boxes':len(data['pending']),
        'remaining_parent_regions':sorted({r for i,b,r in data['pending']}),
        'selected_roots_complete':not pending,'full_eight_point_exclusion_proved':not data['pending'],
        'seconds':time.monotonic()-start},indent=2))


if __name__=='__main__':main()
