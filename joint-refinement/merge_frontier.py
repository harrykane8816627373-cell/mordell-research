"""Compose complete accepted-frontier descendants; arithmetic must be replayed separately."""
from pathlib import Path
import argparse,gzip,hashlib,json
import recent_force as P


def collect(data,base):
    P.descendant(data,base)
    nbase=len(base['nodes']);seen=set();groups={}
    for root,box,parent in base['pending']:
        stack=[root];members=[];open_ids=[]
        while stack:
            i=stack.pop()
            P.need(type(i) is int and 0<=i<len(data['nodes']) and i not in seen
                and (i>=nbase or i==root),'unshared descendant forest for every accepted open box')
            seen.add(i);members.append(i);node=data['nodes'][i]
            if node['kind']=='unknown':open_ids.append(i)
            elif node['kind']=='split':stack.extend([node['left'],node['right']])
            else:P.need(node['kind']=='excluded','known proposal node kind')
        groups[root]={'nodes':members,'closed':not open_ids,'parent':parent,'original_box':box}
    P.need(len(seen)==len(data['nodes'])-base_count(base),'every appended proposal node is reachable')
    return groups


def base_count(base):
    return sum(node['kind']!='unknown' for node in base['nodes'])


def merge(paths,output):
    base,accepted,meta=P.load_recent()
    proposals=[];choices={}
    for path in paths:
        raw=path.read_bytes();data=json.loads(gzip.decompress(raw))
        P.need('accepted_prefix_sha256' not in data or data['accepted_prefix_sha256']==P.SHA256,
            'a proposal may not relabel a different declared accepted prefix')
        data['accepted_prefix_sha256']=P.SHA256
        groups=collect(data,base)
        source={'path':path.name,'sha256':hashlib.sha256(raw).hexdigest()}
        proposals.append((data,groups,source))
        for root,group in groups.items():
            if group['closed'] and (root not in choices or len(group['nodes'])<choices[root][0]):
                choices[root]=(len(group['nodes']),len(proposals)-1)
    result=json.loads(json.dumps(base));result['accepted_prefix_sha256']=P.SHA256
    selected=[]
    for root,original,parent in base['pending']:
        if root not in choices:continue
        count,index=choices[root];data,groups,source=proposals[index];group=groups[root]
        extra=sorted(i for i in group['nodes'] if i!=root)
        offset=len(result['nodes']);mapping={root:root}
        mapping.update({old:offset+j for j,old in enumerate(extra)})
        result['nodes'].extend({'kind':'unknown'} for _ in extra)
        for i in [root]+extra:
            node=dict(data['nodes'][i])
            if node['kind']=='split':node.update(left=mapping[node['left']],right=mapping[node['right']])
            result['nodes'][mapping[i]]=node
        selected.append({'accepted_leaf_id':root,'parent_region_id':parent,'proposed_processed_nodes':count,
            'proposal_sha256':source['sha256']})
    result['pending']=[row for row in result['pending'] if row[0] not in choices]
    result['complete_cover']=not result['pending'];result['full_eight_point_exclusion_proved']=not result['pending']
    P.descendant(result,base)
    raw=gzip.compress((json.dumps(result,separators=(',',':'))+'\n').encode(),compresslevel=9,mtime=0)
    output.write_bytes(raw)
    return {'arithmetic_replay_performed':False,'accepted_prefix_sha256':P.SHA256,
        'proposed_certificate_sha256':hashlib.sha256(raw).hexdigest(),'proposed_certificate_bytes':len(raw),
        'complete_old_frontier_subtrees_selected':len(selected),'selected':selected,
        'new_proposed_processed_nodes':sum(row['proposed_processed_nodes'] for row in selected),
        'old_open_boxes_retained_without_added_partial_descendants':len(result['pending']),
        'new_partial_descendants_imported':0,'proposal_sources':[v[2] for v in proposals],
        'full_eight_point_exclusion_proved':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('proposals',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    print(json.dumps(merge(args.proposals,args.output),indent=2))
