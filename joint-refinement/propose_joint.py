"""Extend exactly the 19 open regions of the accepted 92fbc4d snapshot."""
from pathlib import Path
from collections import Counter
import argparse,gzip,hashlib,json,time
import ellipse_dual as J

HERE=Path(__file__).resolve().parent
PARENT=HERE/'certificates/continuation.json.gz'
PARENT_SHA256='20b66183f129822133dc0f8ab2d5471cd7c2a035e914dc4eb301f43d6d0fff6e'
E=J.E


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes',type=int,default=10000)
    parser.add_argument('--resume',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    raw=PARENT.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PARENT_SHA256:
        raise ArithmeticError('fixed accepted parent snapshot')
    parent=json.loads(gzip.decompress(raw))
    E.setup(parent['precision_bits'])
    settings=E.model(parent['settings']['scope'])
    if args.resume:
        data=json.loads(gzip.decompress(args.resume.read_bytes()))
        if data['parent_sha256']!=PARENT_SHA256 or data['settings']!=parent['settings'] or data['precision_bits']!=E.BITS:
            raise ArithmeticError('same exact joint continuation scope')
        nodes,roots=data['nodes'],data['roots']
        pending=[(i,[tuple(v) for v in b]) for i,b in data['pending']]
    else:
        roots=[[i,j] for j,(i,b) in enumerate(parent['pending'])]
        nodes=[{'kind':'unknown'} for _ in roots]
        pending=[(j,[tuple(v) for v in b]) for j,(i,b) in enumerate(parent['pending'])]
    counts=Counter()
    start=time.monotonic()
    for number in range(args.nodes):
        if not pending:break
        index,original=pending.pop()
        box,reason,aux=J.contract(original,settings)
        if box is None:
            record={'kind':'excluded','reason':reason}
            counts[reason]+=1
        else:
            dim=max(range(13),key=lambda k:box[k][1]-box[k][0])
            lo,hi=box[dim]
            cut=(lo+hi)//2
            if not lo<cut<hi:raise ArithmeticError('exact grid resolution exhausted')
            left,right=list(box),list(box)
            left[dim],right[dim]=(lo,cut),(cut,hi)
            li,ri=len(nodes),len(nodes)+1
            nodes += [{'kind':'unknown'},{'kind':'unknown'}]
            record={'kind':'split','dim':dim,'cut':cut,'left':li,'right':ri}
            pending += [(ri,right),(li,left)]
            counts['split']+=1
        nodes[index]=record
        if (number+1)%2000==0:
            print(json.dumps({'processed':number+1,'open':len(pending),'counts':dict(counts),
                'seconds':time.monotonic()-start}),flush=True)
    data={'parent_sha256':PARENT_SHA256,'precision_bits':E.BITS,'settings':parent['settings'],
        'roots':roots,'nodes':nodes,'pending':pending,'complete_cover':not pending,
        'full_eight_point_exclusion_proved':not pending}
    args.output.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),compresslevel=9,mtime=0))
    print(json.dumps({'counts':dict(counts),'open':len(pending),'seconds':time.monotonic()-start,
        'full_eight_point_exclusion_proved':not pending},indent=2))


if __name__=='__main__':main()
