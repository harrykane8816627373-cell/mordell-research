"""Compose complete accepted-frontier descendants; arithmetic must be replayed separately."""
from pathlib import Path
import argparse,gzip,hashlib,json
import recent_force as P


def collect(data,base,prefix=P):
    prefix.descendant(data,base)
    nbase=len(base['nodes']);seen=set();groups={}
    pending={row[0]:row for row in data['pending']}
    P.need(len(pending)==len(data['pending']),'unique proposed pending IDs')
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
        for i in open_ids:P.need(i in pending and pending[i][2]==parent,'every open proposed descendant stays explicit')
        groups[root]={'nodes':members,'closed':not open_ids,'parent':parent,'original_box':box,
            'open_ids':open_ids,'processed':sum(data['nodes'][i]['kind']!='unknown' for i in members)}
    P.need(len(seen)==len(data['nodes'])-base_count(base),'every appended proposal node is reachable')
    P.need(set(pending)=={i for group in groups.values() for i in group['open_ids']},'no missing or surplus proposed open leaves')
    return groups


def base_count(base):
    return sum(node['kind']!='unknown' for node in base['nodes'])


def merge(paths,output,nov2=False,include_partial=False):
    prefix=P;key='accepted_prefix_sha256'
    if nov2:
        import nov2_force as prefix
        key=prefix.IDENTITY_KEY
    base,accepted,meta=prefix.load_current() if nov2 else prefix.load_recent()
    proposals=[];choices={}
    for path in paths:
        raw=path.read_bytes();data=json.loads(gzip.decompress(raw))
        P.need(key not in data or data[key]==prefix.SHA256,
            'a proposal may not relabel a different declared accepted prefix')
        data[key]=prefix.SHA256
        groups=collect(data,base,prefix)
        source={'path':path.name,'sha256':hashlib.sha256(raw).hexdigest()}
        proposals.append((data,groups,source))
        for root,group in groups.items():
            if not group['closed'] and not (include_partial and group['processed']):continue
            score=(0,len(group['nodes'])) if group['closed'] else (1,-group['processed'])
            if root not in choices or score<choices[root][0]:
                choices[root]=(score,len(proposals)-1)
    result=json.loads(json.dumps(base));result[key]=prefix.SHA256
    selected=[];new_pending=[]
    for root,original,parent in base['pending']:
        if root not in choices:continue
        score,index=choices[root];data,groups,source=proposals[index];group=groups[root]
        extra=sorted(i for i in group['nodes'] if i!=root)
        offset=len(result['nodes']);mapping={root:root}
        mapping.update({old:offset+j for j,old in enumerate(extra)})
        result['nodes'].extend({'kind':'unknown'} for _ in extra)
        for i in [root]+extra:
            node=dict(data['nodes'][i])
            if node['kind']=='split':node.update(left=mapping[node['left']],right=mapping[node['right']])
            result['nodes'][mapping[i]]=node
        proposed_pending={row[0]:row for row in data['pending']}
        for i in group['open_ids']:
            row=proposed_pending[i];new_pending.append([mapping[i],row[1],row[2]])
        selected.append({'accepted_leaf_id':root,'parent_region_id':parent,'proposed_processed_nodes':group['processed'],
            'proposal_sha256':source['sha256'],'closed':group['closed'],'open_descendants':len(group['open_ids'])})
    result['pending']=[row for row in result['pending'] if row[0] not in choices]+new_pending
    result['complete_cover']=not result['pending'];result['full_eight_point_exclusion_proved']=not result['pending']
    prefix.descendant(result,base)
    raw=gzip.compress((json.dumps(result,separators=(',',':'))+'\n').encode(),compresslevel=9,mtime=0)
    output.write_bytes(raw)
    return {'arithmetic_replay_performed':False,'accepted_prefix_identity_key':key,'accepted_prefix_sha256':prefix.SHA256,
        'proposed_certificate_sha256':hashlib.sha256(raw).hexdigest(),'proposed_certificate_bytes':len(raw),
        'complete_old_frontier_subtrees_selected':sum(row['closed'] for row in selected),'selected':selected,
        'partial_old_frontier_subtrees_selected':sum(not row['closed'] for row in selected),
        'new_proposed_processed_nodes':sum(row['proposed_processed_nodes'] for row in selected),
        'old_open_boxes_retained_without_added_partial_descendants':len(base['pending'])-len(selected),
        'new_partial_descendants_imported':len(new_pending),'proposal_sources':[v[2] for v in proposals],
        'full_eight_point_exclusion_proved':not result['pending']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('proposals',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--nov2-prefix',action='store_true',help='compose only closed descendants of accepted 19dbbd3')
    parser.add_argument('--include-partial',action='store_true',help='retain exact open descendants of partial proposals; never mark them closed')
    args=parser.parse_args()
    print(json.dumps(merge(args.proposals,args.output,args.nov2_prefix,args.include_partial),indent=2))
