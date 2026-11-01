"""Verify prior-snapshot identity and summarize closures after full cover replay."""
from pathlib import Path
import gzip
import hashlib
import json

HERE=Path(__file__).resolve().parent


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def summarize():
    metadata=json.loads((HERE/'certificates/prior_frontier.json').read_text())
    raw=(HERE/'certificates/continuation.json.gz').read_bytes()
    current=json.loads(gzip.decompress(raw))
    cover=json.loads((HERE/'reports/cover.json').read_text())
    need(cover['passed'] is True and cover['continuation_sha256']==hashlib.sha256(raw).hexdigest(),
        'prerequisite complete cover replay report for the same certificate')
    frontier=metadata['pending']
    old_ids={i for i,b in frontier}
    need(len(old_ids)==len(frontier) and metadata['node_count']-len(frontier)==metadata['processed_node_count'],
        'prior frontier bookkeeping')
    prefix=[dict(node) for node in current['nodes'][:metadata['node_count']]]
    for i in old_ids:
        prefix[i]={'kind':'unknown'}
    reconstructed={key:current[key] for key in ['base_sha256','precision_bits','settings','roots']}
    reconstructed.update({'nodes':prefix,'pending':frontier,
        'complete_cover':False,'full_eight_point_exclusion_proved':False})
    encoded=(json.dumps(reconstructed,separators=(',',':'))+'\n').encode('utf-8')
    need(hashlib.sha256(encoded).hexdigest()==metadata['uncompressed_sha256'],
        'all old processed nodes, root settings and frontier boxes reconstruct the exact prior snapshot')
    seen,all_open,rows=set(),set(),[]
    for index,box in frontier:
        stack,opened,excluded=[index],[],0
        while stack:
            i=stack.pop()
            need(type(i) is int and 0<=i<len(current['nodes']) and i not in seen,
                'disjoint acyclic subtrees of old open regions')
            need(i>=metadata['node_count'] or i in old_ids,'new subtrees do not reuse old processed nodes')
            seen.add(i)
            node=current['nodes'][i]
            if node['kind']=='split':
                stack += [node['left'],node['right']]
            elif node['kind']=='unknown':
                opened.append(i)
                all_open.add(i)
            else:
                need(node['kind']=='excluded','known accepted node kind')
                excluded+=1
        rows.append({'prior_frontier_id':index,'closed':not opened,
            'excluded_descendant_leaves':excluded,'open_descendants':sorted(opened)})
    need(all_open=={i for i,b in current['pending']}
        and len(seen)==len(current['nodes'])-metadata['processed_node_count'],
        'every appended node and open leaf remains accounted for')
    return {'passed':True,'proof_replay_performed':False,'requires_complete_cover_replay':True,
        'prior_commit':metadata['prior_commit'],'prior_certificate_sha256':metadata['compressed_sha256'],
        'prior_snapshot_reconstructed_exactly':True,'current_certificate_sha256':hashlib.sha256(raw).hexdigest(),
        'prior_regions':len(rows),'closed_prior_regions':sum(r['closed'] for r in rows),
        'unclosed_prior_regions':sum(not r['closed'] for r in rows),'current_open_descendants':len(all_open),
        'scope':'Coordinate regions; no stationary point classification or covered-volume percentage.',
        'rows':rows,'full_eight_point_exclusion_proved':cover['full_eight_point_exclusion_proved']}


if __name__=='__main__':
    print(json.dumps(summarize(),indent=2))
