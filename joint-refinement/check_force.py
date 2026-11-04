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


def check(path,progress=False,reuse_checkpoint=False,reuse_recent=False,reuse_nov2=False):
    parent,old_report,settings,meta=parent_data()
    initial=initialize()
    raw=path.read_bytes()
    data=json.loads(gzip.decompress(raw))
    need(data['parent_sha256']==A.PARENT_SHA256 and data['attachment_sha256']==A.ZIP_SHA256
        and data['precision_bits']==40 and data['scope']=='global_maximum' and data['kernel_version']==1,
        'bound force model, attachment and parent')
    need(data['roots']==initial['roots'],'all original twenty-one parent regions are covered once')
    nodes=data['nodes']
    checkpoint=None
    if 'checkpoint_sha256' in data or reuse_checkpoint:
        import checkpoint_force as P
        checkpoint=P.load_checkpoint()
        P.descendant(data,checkpoint[0])
    recent=None
    if 'accepted_prefix_sha256' in data or reuse_recent:
        import recent_force as N
        recent=N.load_recent()
        N.descendant(data,recent[0])
    nov2=None
    if 'nov2_prefix_sha256' in data or reuse_nov2:
        import nov2_force as V
        nov2=V.load_current()
        V.descendant(data,nov2[0])
    need(sum([reuse_checkpoint,reuse_recent,reuse_nov2])<=1,'select only one reused prefix')
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
    accepted_frontier={}
    if recent:
        accepted_frontier={i:{'accepted_leaf_id':i,'parent_region_id':r,'excluded_leaves_since_prefix':0,'open_leaves':[]}
            for i,b,r in recent[0]['pending']}
    nov2_frontier={}
    if nov2:
        nov2_frontier={i:{'accepted_leaf_id':i,'parent_region_id':r,'excluded_leaves_since_prefix':0,'open_leaves':[]}
            for i,b,r in nov2[0]['pending']}
    stack=[(j,old[root],root,None,None) for root,j in data['roots']]
    seen=set()
    counts=Counter()
    engines=Counter()
    if reuse_checkpoint or reuse_recent or reuse_nov2:
        base,accepted,_=nov2 if reuse_nov2 else recent if reuse_recent else checkpoint
        seen={i for i,node in enumerate(base['nodes']) if node['kind']!='unknown'}
        counts=Counter({k:v for k,v in accepted['counts'].items() if k!='unresolved'})
        engines=Counter(n['engine'] for n in base['nodes'] if n['kind']!='unknown')
        stack=[(i,[tuple(p) for p in b],r,i if reuse_recent or reuse_nov2 else None,i if reuse_nov2 else None)
            for i,b,r in base['pending']]
        for row in accepted['frontier']:
            frontier[row['parent_region_id']]['excluded_leaves']=row['excluded_leaves']
        if reuse_nov2:
            for row in accepted['accepted_frontier']:
                accepted_frontier[row['accepted_leaf_id']]['excluded_leaves_since_prefix']=row['excluded_leaves_since_prefix']
    while stack:
        i,original,root,origin,nov2_origin=stack.pop()
        if i in accepted_frontier:origin=i
        if i in nov2_frontier:nov2_origin=i
        need(type(i) is int and 0<=i<len(nodes) and i not in seen,'disjoint acyclic force forest')
        seen.add(i)
        if progress and len(seen)%2000==0:print('Replayed force nodes:',len(seen),'of',len(nodes),file=sys.stderr,flush=True)
        node=nodes[i]
        if node['kind']=='unknown':
            need(i in pending and pending[i]==(original,root),'every unresolved coordinate region stays explicit')
            frontier[root]['open_leaves'].append(i)
            if recent:
                need(origin in accepted_frontier,'every current open leaf descends from one accepted open region')
                accepted_frontier[origin]['open_leaves'].append(i)
            if nov2:
                need(nov2_origin in nov2_frontier,'every current open leaf descends from one Nov 2 open region')
                nov2_frontier[nov2_origin]['open_leaves'].append(i)
            counts['unresolved']+=1
            continue
        need(i not in pending and node['engine'] in ['uploaded','joint_force','joint_force_global','joint_force_edge','joint_force_phase','joint_force_gram','joint_force_subset','joint_force_subset_fast','joint_force_target_radius','joint_force_schwarz2_radius','joint_force_ellipse_schur_radius','joint_force_translated_subset'],'known processed engine')
        engines[node['engine']]+=1
        if node['engine']=='uploaded':
            box,reason=original,None
            if node['kind']=='excluded':
                need(A.excludes(box,node['reason']),'independent exact uploaded whole-box rejection')
                reason=node['reason']
        elif node['engine']=='joint_force_edge':
            import edge_force
            box,reason,aux=edge_force.contract(original,settings)
        elif node['engine']=='joint_force_phase':
            import phase_force
            box,reason,aux=phase_force.contract(original,settings)
        elif node['engine']=='joint_force_gram':
            import gram_force
            box,reason,aux=gram_force.contract(original,settings)
        elif node['engine']=='joint_force_subset':
            import subset_force
            box,reason,aux=subset_force.contract(original,settings)
        elif node['engine']=='joint_force_subset_fast':
            import subset_force
            box,reason,aux=subset_force.contract(original,settings,passes=1)
        elif node['engine']=='joint_force_target_radius':
            import target_radius
            box,reason,aux=target_radius.contract(original,settings)
        elif node['engine']=='joint_force_schwarz2_radius':
            import schwarz2_radius
            box,reason,aux=schwarz2_radius.contract(original,settings)
        elif node['engine']=='joint_force_ellipse_schur_radius':
            import ellipse_schur_radius
            box,reason,aux=ellipse_schur_radius.contract(original,settings)
        elif node['engine']=='joint_force_translated_subset':
            import translated_subset
            box,reason,aux=translated_subset.contract(original,settings)
        else:
            endpoint=settings['separation'] if node['engine']=='joint_force_global' else None
            box,reason,aux=K.contract(original,settings,endpoint_separation=endpoint)
        if node['kind']=='excluded':
            need(reason==node['reason'] and (node['engine']=='uploaded' or box is None),'verified entire feasible leaf exclusion')
            counts[node['engine']+'/'+reason]+=1
            frontier[root]['excluded_leaves']+=1
            if origin is not None:accepted_frontier[origin]['excluded_leaves_since_prefix']+=1
            if nov2_origin is not None:nov2_frontier[nov2_origin]['excluded_leaves_since_prefix']+=1
        elif node['kind']=='split':
            need(box is not None,'new split follows only sound contractions')
            dim,cut=node['dim'],node['cut']
            need(type(dim) is int and 0<=dim<13 and type(cut) is int and box[dim][0]<cut<box[dim][1],
                'proper complete exact binary split')
            left,right=list(box),list(box)
            left[dim],right[dim]=(box[dim][0],cut),(cut,box[dim][1])
            stack += [(node['left'],left,root,origin,nov2_origin),(node['right'],right,root,origin,nov2_origin)]
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
    accepted_rows=[]
    for row in accepted_frontier.values():
        row['open_leaves'].sort();row['closed']=not row['open_leaves'];accepted_rows.append(row)
    closed_accepted=sorted(row['accepted_leaf_id'] for row in accepted_rows if row['closed'])
    nov2_rows=[]
    for row in nov2_frontier.values():
        row['open_leaves'].sort();row['closed']=not row['open_leaves'];nov2_rows.append(row)
    closed_nov2=sorted(row['accepted_leaf_id'] for row in nov2_rows if row['closed'])
    return {'passed':True,'parent_commit':meta['parent_commit'],'parent_certificate_sha256':A.PARENT_SHA256,
        'attachment_sha256':A.ZIP_SHA256,'force_certificate_sha256':hashlib.sha256(raw).hexdigest(),
        'scope':'global_maximum','parent_frontier_regions':21,'force_tree_nodes':len(nodes),
        'processed_nodes':len(nodes)-len(pending),'imported_processed_nodes':imported_processed,
        'new_processed_nodes':sum(v for k,v in engines.items() if k!='uploaded'),
        'endpoint_distance_processed_nodes':engines['joint_force_global'],'counts':dict(sorted(counts.items())),
        'edge_correlated_processed_nodes':engines['joint_force_edge'],
        'complex_remainder_processed_nodes':engines['joint_force_phase'],
        'pair_Gram_processed_nodes':engines['joint_force_gram'],
        'subset_budget_processed_nodes':engines['joint_force_subset'],
        'subset_budget_fast_processed_nodes':engines['joint_force_subset_fast'],
        'target_radius_processed_nodes':engines['joint_force_target_radius'],
        'schwarz2_radius_processed_nodes':engines['joint_force_schwarz2_radius'],
        'ellipse_schur_radius_processed_nodes':engines['joint_force_ellipse_schur_radius'],
        'translated_subset_processed_nodes':engines['joint_force_translated_subset'],
        'nov2_prefix_sha256':data.get('nov2_prefix_sha256'),
        'nov2_prefix_processed_nodes':nov2[1]['processed_nodes'] if nov2 else 0,
        'processed_since_nov2_prefix':len(nodes)-len(pending)-nov2[1]['processed_nodes'] if nov2 else None,
        'newly_closed_since_nov2_prefix':sorted(set(closed)-set(nov2[1]['closed_parent_ids'])) if nov2 else [],
        'nov2_frontier_regions':len(nov2_rows),'closed_nov2_frontier_ids':closed_nov2,
        'remaining_nov2_frontier_regions':len(nov2_rows)-len(closed_nov2),'nov2_frontier':nov2_rows,
        'accepted_prefix_sha256':data.get('accepted_prefix_sha256'),
        'accepted_prefix_processed_nodes':recent[1]['processed_nodes'] if recent else 0,
        'processed_since_accepted_prefix':len(nodes)-len(pending)-recent[1]['processed_nodes'] if recent else None,
        'newly_closed_since_accepted_prefix':sorted(set(closed)-set(recent[1]['closed_parent_ids'])) if recent else [],
        'accepted_frontier_regions':len(accepted_rows),'closed_accepted_frontier_ids':closed_accepted,
        'remaining_accepted_frontier_regions':len(accepted_rows)-len(closed_accepted),
        'accepted_frontier':accepted_rows,
        'checkpoint_sha256':data.get('checkpoint_sha256'),
        'checkpoint_processed_nodes':checkpoint[1]['processed_nodes'] if checkpoint else 0,
        'processed_since_checkpoint':len(nodes)-len(pending)-checkpoint[1]['processed_nodes'] if checkpoint else None,
        'newly_closed_since_checkpoint':sorted(set(closed)-set(checkpoint[1]['closed_parent_ids'])) if checkpoint else [],
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
    parser.add_argument('--reuse-checkpoint',action='store_true',
        help='only when the exact archived 768598a proof is already accepted; checks all 50 unchanged arithmetic inputs')
    parser.add_argument('--reuse-accepted-prefix',action='store_true',
        help='only when the exact archived 8be9f91 proof is already accepted; checks all 61 unchanged arithmetic inputs')
    parser.add_argument('--reuse-nov2-prefix',action='store_true',
        help='only when the exact archived 19dbbd3 proof is already accepted; checks all 70 unchanged arithmetic inputs')
    args=parser.parse_args()
    print(json.dumps(check(args.certificate,args.progress,args.reuse_checkpoint,args.reuse_accepted_prefix,args.reuse_nov2_prefix),indent=2))
