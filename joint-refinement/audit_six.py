"""Independent assert-free replay of the user's six-closure attachment."""
from pathlib import Path
from collections import Counter,deque
import gzip,hashlib,json,zipfile
import force_enclosure as F

HERE=Path(__file__).resolve().parent
E=F.E
ZIP_SHA256='87c08af482404117068462cf7151464f4685509523bb64f0e998ae7ca48204fb'
PARENT_SHA256='44562c49f23e2fb322b91a0b45279cd4da5581644b07be6bf2b1d7c052f2bc2b'
CLOSED=(9976,10012,10014,10016,10019,10020)
PARTIAL=((7954,1800),(9814,4000))


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def input_data(path):
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==ZIP_SHA256,'fixed original user attachment')
    parent_raw=(HERE/'certificates/joint_extension.json.gz').read_bytes()
    need(hashlib.sha256(parent_raw).hexdigest()==PARENT_SHA256,'fixed e846b8e parent frontier')
    parent=json.loads(gzip.decompress(parent_raw))
    archive=zipfile.ZipFile(path)
    frontier_raw=archive.read('frontier21.json')
    frontier=json.loads(frontier_raw)
    need(frontier==parent['pending'] and len(frontier)==21,'all twenty-one uploaded boxes equal the exact parent boxes in order')
    E.setup(parent['precision_bits'])
    return archive,frontier,hashlib.sha256(frontier_raw).hexdigest()


def excludes(box,reason):
    q=E.positions(box)
    if reason in ('force_affine','force_taylor'):
        return F.exclusion(F.enclose(box),reason.split('_')[1]) is not None
    if reason=='radius':return any(E.cnorm(z)[0]*32>49*E.S for z in q)
    if reason=='separation':return any(E.cnorm(E.cs(q[i],q[j]))[1]*2<=E.S for i,j in E.PAIRS)
    if reason=='ellipse':
        s2,s3=E.ci(0),E.ci(0)
        for z in q:
            z2=E.cm(z,z)
            s2=E.ca(s2,z2)
            s3=E.ca(s3,E.cm(z2,z))
        return 175*E.cnorm(s2)[0]+54*E.cnorm(s3)[0]>=175*E.S
    raise ArithmeticError('unimplemented imported proof reason '+str(reason))


def replay(data,frontier,frontier_sha,expected_root):
    need(data['input_file_sha256']==frontier_sha and data['precision_bits']==40
        and data['algorithm_version']=='affine','bound imported algorithm and exact input')
    need(data['original_roots']==[expected_root],'this imported certificate covers exactly its named old region')
    original=dict((i,[tuple(p) for p in b]) for i,b in frontier)[expected_root]
    pending=deque([original])
    counts=Counter()
    for node in data['proof_nodes']:
        need(bool(pending),'no proof steps after a region is already closed')
        box=pending.popleft()
        need(node[1]==expected_root and [list(p) for p in box]==node[2],'exact next whole box')
        if node[0]=='X':
            need(excludes(box,node[3]),'independent whole-box exclusion '+str(node[3]))
            counts[node[3]]+=1
        elif node[0]=='S':
            dim,cut=node[3:]
            need(type(dim) is int and 0<=dim<13 and type(cut) is int and box[dim][0]<cut<box[dim][1],
                'proper exact imported split')
            left,right=list(box),list(box)
            left[dim],right[dim]=(box[dim][0],cut),(cut,box[dim][1])
            pending.extend([left,right])
            counts['split']+=1
        else:raise ArithmeticError('unknown imported proof step')
    recorded=[(root,[tuple(p) for p in box]) for root,box in data['pending']]
    need(recorded==[(expected_root,b) for b in pending],'every unresolved imported box is retained')
    need(dict(counts)==data['counts'] and data['open_count']==len(pending),'exact imported accounting')
    need(data['selected_roots_complete'] is (not pending) and data['full_exclusion_proved'] is False,
        'a selected root closure is not a proof of the complete eight-point problem')
    need(data['closed_roots']==([] if pending else [expected_root]),'correct closed-root claim')
    return {'passed':True,'root':expected_root,'processed_nodes':len(data['proof_nodes']),
        'counts':dict(sorted(counts.items())),'pending_boxes':len(pending),'selected_scope_complete':not pending}


def audit(path):
    archive,frontier,sha=input_data(path)
    reports=[]
    for root in CLOSED:
        data=json.loads(gzip.decompress(archive.read('certificates/closed_'+str(root)+'.json.gz')))
        result=replay(data,frontier,sha,root)
        need(result['selected_scope_complete'],'each of the six named closure proofs is complete')
        reports.append(result)
    partial=[]
    for root,n in PARTIAL:
        data=json.loads(gzip.decompress(archive.read('certificates/open_'+str(root)+'_'+str(n)+'.json.gz')))
        result=replay(data,frontier,sha,root)
        need(not result['selected_scope_complete'],'partial input trees are not silently closed')
        partial.append(result)
    return {'passed':True,'attachment_sha256':ZIP_SHA256,'parent_sha256':PARENT_SHA256,
        'all_uploaded_boxes_equal_parent_boxes':True,'closed_parent_regions':list(CLOSED),
        'remaining_parent_regions':sorted(set(i for i,b in frontier)-set(CLOSED)),
        'complete_proof_nodes':sum(r['processed_nodes'] for r in reports),'closure_proofs':reports,
        'partial_proofs':partial,'uses_assert_for_acceptance':False,'uses_floating_point_for_acceptance':False,
        'full_parent_arithmetic_replayed':False,'full_eight_point_exclusion_proved':False}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('attachment',type=Path)
    args=parser.parse_args()
    print(json.dumps(audit(args.attachment),indent=2))
