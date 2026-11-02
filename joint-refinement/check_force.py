"""Replay the complete uploaded-prefix plus joint-force descendant forest."""
from pathlib import Path,PurePosixPath
from collections import Counter
import argparse,gzip,hashlib,json,sys
import audit_six as A
import force_contract as K
from propose_force import initialize,ATTACHMENT

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
E=K.E


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def parent_data():
    meta=json.loads((HERE/'certificates/force_parent.json').read_text())
    need(meta['parent_commit']=='e846b8e9be902a636787974fba1e219cddc8f495'
        and meta['parent_certificate_sha256']==A.PARENT_SHA256,'fixed accepted force parent identity')
    names=set()
    for row in meta['arithmetic_inputs']:
        name=row['path']
        path=PurePosixPath(name)
        need(not path.is_absolute() and '..' not in path.parts and ':' not in name and '\\' not in name
            and name not in names and path.as_posix()==name,'safe unique legacy arithmetic input')
        names.add(name)
        raw=(ROOT/name).read_bytes()
        need(len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],
            'unchanged accepted parent arithmetic '+name)
    need(len(names)==meta['arithmetic_input_count'],'all declared legacy arithmetic inputs checked')
    raw=(HERE/'certificates/joint_extension.json.gz').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==A.PARENT_SHA256,'fixed parent certificate bytes')
    report_raw=(HERE/'reports/joint_cover.json').read_bytes()
    need(hashlib.sha256(report_raw).hexdigest()==meta['parent_report_sha256'],'fixed original parent report')
    report=json.loads(report_raw)
    need(report['passed'] is True and report['extension_sha256']==A.PARENT_SHA256
        and report['unresolved_regions']==21 and report['full_eight_point_exclusion_proved'] is False,
        'parent proof remains accepted and incomplete')
    parent=json.loads(gzip.decompress(raw))
    E.setup(parent['precision_bits'])
    settings=E.model('global_maximum')
    need(parent['settings']==json.loads(json.dumps(settings)),'same physical normalization and global-maximum scope')
    return parent,report,settings,meta


def check(path,progress=False):
    parent,old_report,settings,meta=parent_data()
    initial=initialize()
    raw=path.read_bytes()
    data=json.loads(gzip.decompress(raw))
    need(data['parent_sha256']==A.PARENT_SHA256 and data['attachment_sha256']==A.ZIP_SHA256
        and data['precision_bits']==40 and data['scope']=='global_maximum' and data['kernel_version']==1,
        'bound force model, attachment and parent')
    need(data['roots']==initial['roots'],'all original twenty-one parent regions are covered once')
    nodes=data['nodes']
    need(len(nodes)>=len(initial['nodes']),'entire imported forest retained')
    imported_processed=0
    for i,node in enumerate(initial['nodes']):
        if node['kind']!='unknown':
            need(nodes[i]==node,'every imported proof step is preserved verbatim in the composed graph')
            imported_processed+=1
    pending={}
    for i,box,root in data['pending']:
        need(type(i) is int and i not in pending and type(root) is int,'unique typed pending leaf')
        need(len(box)==13 and all(len(p)==2 and all(type(v) is int for v in p) and p[0]<=p[1] for p in box),
            'exact pending coordinate box')
        pending[i]=([tuple(p) for p in box],root)
    complete=not pending
    need(data['complete_cover'] is complete and data['full_eight_point_exclusion_proved'] is complete,
        'open leaves prohibit complete exclusion')
    old={i:[tuple(p) for p in b] for i,b in parent['pending']}
    frontier={root:{'parent_region_id':root,'excluded_leaves':0,'open_leaves':[]} for root in old}
    stack=[(j,old[root],root) for root,j in data['roots']]
    seen=set()
    counts=Counter()
    engines=Counter()
    newtons=0
    while stack:
        i,original,root=stack.pop()
        need(type(i) is int and 0<=i<len(nodes) and i not in seen,'disjoint acyclic force forest')
        seen.add(i)
        if progress and len(seen)%2000==0:print('Replayed force nodes:',len(seen),'of',len(nodes),file=sys.stderr,flush=True)
        node=nodes[i]
        if node['kind']=='unknown':
            need(i in pending and pending[i]==(original,root),'every unresolved coordinate region stays explicit')
            frontier[root]['open_leaves'].append(i)
            counts['unresolved']+=1
            continue
        need(i not in pending and node['engine'] in ['uploaded','joint_force','joint_force_global'],'known processed engine')
        engines[node['engine']]+=1
        if node['engine']=='uploaded':
            box,reason=original,None
            if node['kind']=='excluded':
                need(A.excludes(box,node['reason']),'independent exact uploaded whole-box rejection')
                reason=node['reason']
        else:
            endpoint=settings['separation'] if node['engine']=='joint_force_global' else None
            box,reason,aux=K.contract(original,settings,endpoint_separation=endpoint)
        if node['kind']=='excluded':
            need(reason==node['reason'] and (node['engine']=='uploaded' or box is None),'verified entire feasible leaf exclusion')
            counts[node['engine']+'/'+reason]+=1
            frontier[root]['excluded_leaves']+=1
        elif node['kind']=='split':
            need(box is not None,'new split follows only sound contractions')
            dim,cut=node['dim'],node['cut']
            need(type(dim) is int and 0<=dim<13 and type(cut) is int and box[dim][0]<cut<box[dim][1],
                'proper complete exact binary split')
            left,right=list(box),list(box)
            left[dim],right[dim]=(box[dim][0],cut),(cut,box[dim][1])
            stack += [(node['left'],left,root),(node['right'],right,root)]
            counts['split']+=1
        else:raise ArithmeticError('unknown force proof node')
    need(len(seen)==len(nodes) and counts['unresolved']==len(pending),'every proof node and open leaf accounted for')
    rows=[]
    for root in old:
        row=frontier[root]
        row['open_leaves'].sort()
        row['closed']=not row['open_leaves']
        rows.append(row)
    closed=sorted(row['parent_region_id'] for row in rows if row['closed'])
    remaining=sorted(set(old)-set(closed))
    excluded=sum(v for k,v in counts.items() if k not in ['split','unresolved'])
    return {'passed':True,'parent_commit':meta['parent_commit'],'parent_certificate_sha256':A.PARENT_SHA256,
        'attachment_sha256':A.ZIP_SHA256,'force_certificate_sha256':hashlib.sha256(raw).hexdigest(),
        'scope':'global_maximum','parent_frontier_regions':21,'force_tree_nodes':len(nodes),
        'processed_nodes':len(nodes)-len(pending),'imported_processed_nodes':imported_processed,
        'new_processed_nodes':engines['joint_force']+engines['joint_force_global'],
        'endpoint_distance_processed_nodes':engines['joint_force_global'],'counts':dict(sorted(counts.items())),
        'uploaded_six_closed_regions':list(A.CLOSED),'total_closed_parent_regions':len(closed),
        'closed_parent_ids':closed,'newly_closed_beyond_upload':sorted(set(closed)-set(A.CLOSED)),
        'remaining_parent_regions':len(remaining),'remaining_parent_ids':remaining,'open_leaf_boxes':len(pending),
        'combined_tree_nodes':old_report['combined_tree_nodes']-21+len(nodes),
        'combined_excluded_leaves':old_report['combined_excluded_leaves']+excluded,
        'frontier':rows,'unchanged_parent_arithmetic_inputs':meta['arithmetic_input_count'],
        'parent_arithmetic_replayed_by_this_checker':False,'complete_cover':complete,'full_eight_point_exclusion_proved':complete}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    parser.add_argument('--progress',action='store_true')
    args=parser.parse_args()
    print(json.dumps(check(args.certificate,args.progress),indent=2))
