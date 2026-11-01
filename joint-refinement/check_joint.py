"""Replay every new descendant of the fixed, previously accepted 19-box frontier."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import argparse,gzip,hashlib,json,sys
import ellipse_dual as J
from propose_joint import PARENT,PARENT_SHA256

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
E=J.E


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def accepted_parent():
    raw=PARENT.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==PARENT_SHA256,'fixed accepted parent certificate')
    meta=json.loads((HERE/'certificates/accepted_parent.json').read_text())
    need(meta['certificate_sha256']==PARENT_SHA256,'accepted-parent identity')
    for row in meta['unchanged_arithmetic_inputs']:
        data=(ROOT/row['path']).read_bytes()
        need(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],
            'unchanged parent arithmetic input '+row['path'])
    report_raw=(HERE/'reports/cover.json').read_bytes()
    need(hashlib.sha256(report_raw).hexdigest()==meta['cover_report_sha256'],'original parent replay report')
    report=json.loads(report_raw)
    need(report['passed'] is True and report['continuation_sha256']==PARENT_SHA256
        and report['unresolved_regions']==19 and report['full_eight_point_exclusion_proved'] is False,
        'accepted parent still explicitly has nineteen open regions')
    parent=json.loads(gzip.decompress(raw))
    E.setup(parent['precision_bits'])
    settings=E.model(parent['settings']['scope'])
    need(parent['settings']==json.loads(json.dumps(settings)),'same exact physical scope')
    return parent,report,settings


def check(path,progress=False):
    parent,old_report,settings=accepted_parent()
    raw=path.read_bytes()
    data=json.loads(gzip.decompress(raw))
    need(data['parent_sha256']==PARENT_SHA256 and data['precision_bits']==E.BITS
        and data['settings']==parent['settings'],'same accepted parent and interval model')
    old={i:[tuple(v) for v in box] for i,box in parent['pending']}
    roots=data['roots']
    need(len(roots)==19 and {i for i,j in roots}==set(old)
        and len({j for i,j in roots})==19,'each of the nineteen old open boxes occurs once')
    pending={}
    for i,box in data['pending']:
        need(type(i) is int and i not in pending,'unique open leaf index')
        pending[i]=[tuple(v) for v in box]
    complete=not pending
    need(data['complete_cover'] is complete and data['full_eight_point_exclusion_proved'] is complete,
        'open leaves prevent a full exclusion claim')
    rows={i:{'parent_region_id':i,'excluded_descendant_leaves':0,'open_descendants':[]} for i in old}
    stack=[(j,old[i],i) for i,j in roots]
    nodes=data['nodes']
    seen=set()
    counts=Counter()
    split_bounds=[]
    parameter_x=[]
    parameter_y=[]
    centers=Counter()
    while stack:
        i,original,root=stack.pop()
        need(type(i) is int and 0<=i<len(nodes) and i not in seen,'acyclic forest without shared or omitted nodes')
        seen.add(i)
        if progress and len(seen)%1000==0:
            print('Replayed joint extension nodes:',len(seen),'of',len(nodes),file=sys.stderr,flush=True)
        node=nodes[i]
        if node['kind']=='unknown':
            need(i in pending and pending[i]==original,'every open box remains explicit')
            counts['unresolved']+=1
            rows[root]['open_descendants'].append(i)
            continue
        need(i not in pending,'processed nodes cannot also be pending')
        box,reason,aux=J.contract(original,settings)
        if node['kind']=='excluded':
            need(box is None and reason==node['reason'],'exact whole-leaf joint exclusion')
            counts[reason]+=1
            rows[root]['excluded_descendant_leaves']+=1
            if reason=='ellipse_split_log':
                split_bounds.append(aux['log_Delta_over_target_upper'])
                parameter_x.append(aux['mass_parameters'][0])
                parameter_y.append(aux['mass_parameters'][1])
                centers[aux['center_source']]+=1
        elif node['kind']=='split':
            need(box is not None,'partition follows valid joint contractions')
            dim,cut=node['dim'],node['cut']
            need(type(dim) is int and 0<=dim<13 and type(cut) is int and box[dim][0]<cut<box[dim][1],
                'proper exact interval split')
            left,right=list(box),list(box)
            left[dim],right[dim]=(box[dim][0],cut),(cut,box[dim][1])
            stack += [(node['left'],left,root),(node['right'],right,root)]
            counts['split']+=1
        else:raise ArithmeticError('unknown extension node kind')
    need(len(seen)==len(nodes) and counts['unresolved']==len(pending),'all new nodes and pending leaves accounted for')
    excluded=sum(v for k,v in counts.items() if k not in ['split','unresolved'])
    frontier=[]
    for i in old:
        row=rows[i]
        row['open_descendants'].sort()
        row['closed']=not row['open_descendants']
        frontier.append(row)
    return {'passed':True,'parent_commit':'92fbc4ddbdca4e37f081c061d1a2c49c25af4fed',
        'parent_certificate_sha256':PARENT_SHA256,'extension_sha256':hashlib.sha256(raw).hexdigest(),
        'scope':settings['scope'],'precision_bits':E.BITS,'parent_frontier_regions':19,
        'extension_nodes':len(nodes),'processed_extension_nodes':len(nodes)-len(pending),
        'extension_counts':dict(sorted(counts.items())),'additional_excluded_leaves':excluded,
        'combined_tree_nodes':old_report['combined_tree_nodes']-19+len(nodes),
        'combined_excluded_leaves':old_report['combined_excluded_leaves']+excluded,
        'closed_parent_regions':sum(row['closed'] for row in frontier),
        'unclosed_parent_regions':sum(not row['closed'] for row in frontier),
        'unresolved_regions':len(pending),'frontier':frontier,
        'new_ellipse_split_exclusions':len(split_bounds),'ellipse_center_counts':dict(sorted(centers.items())),
        'largest_ellipse_split_log_upper':str(F(max(split_bounds),E.S)) if split_bounds else None,
        'real_mass_parameter_range':[str(F(min(parameter_x),E.S)),str(F(max(parameter_x),E.S))] if parameter_x else None,
        'imag_mass_parameter_range':[str(F(min(parameter_y),E.S)),str(F(max(parameter_y),E.S))] if parameter_y else None,
        'new_ellipse_mass_feasibility_exclusions':counts['ellipse_mass_split'],
        'new_fourth_mode_exclusions':sum(counts[k] for k in ['Hessian_inner_particle','Hessian_radial_mass','Hessian_radial_budget']),
        'complete_cover':complete,'full_eight_point_exclusion_proved':complete}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    parser.add_argument('--progress',action='store_true')
    args=parser.parse_args()
    print(json.dumps(check(args.certificate,args.progress),indent=2))
