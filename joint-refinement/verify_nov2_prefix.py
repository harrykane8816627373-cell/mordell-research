"""Check the accepted Nov 2 prefix binding and reject forged continuations."""
from pathlib import Path
import copy,gzip,json,tempfile
import nov2_force as P
from check_force import check


def verify():
    base,accepted,meta=P.load_current()
    baseline=copy.deepcopy(base);baseline[P.IDENTITY_KEY]=P.SHA256
    rejected=0
    with tempfile.TemporaryDirectory(prefix='mordell-prefix-test-') as temp:
        path=Path(temp)/'candidate.json.gz'
        def save(data):
            path.write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))
        save(baseline)
        result=check(path,reuse_nov2=True)
        P.need(result['processed_since_nov2_prefix']==0 and result['closed_nov2_frontier_ids']==[]
            and result['open_leaf_boxes']==14 and result['remaining_parent_ids']==list(range(10))
            and result['closed_accepted_frontier_ids']==accepted['closed_accepted_frontier_ids']
            and result['counts']==accepted['counts'],'zero-new-step reuse retains the complete accepted accounting')
        changes=[
            lambda d:d.update(nov2_prefix_sha256='0'*64),
            lambda d:d.update(scope='user_physical_system'),
            lambda d:d.update(precision_bits=41),
            lambda d:d['roots'].pop(),
            lambda d:d['nodes'][next(i for i,n in enumerate(base['nodes']) if n['kind']!='unknown')].update(kind='unknown'),
            lambda d:d['pending'].pop(),
            lambda d:d['pending'][0][1][0].__setitem__(0,d['pending'][0][1][0][0]+1),
            lambda d:d['pending'].append(copy.deepcopy(d['pending'][0])),
            lambda d:d.update(complete_cover=True,full_eight_point_exclusion_proved=True),
        ]
        for change in changes:
            forged=copy.deepcopy(baseline);change(forged);save(forged)
            try:check(path,reuse_nov2=True)
            except ArithmeticError:rejected+=1
            else:raise ArithmeticError('forged Nov 2 continuation accepted')
    return {'passed':True,'accepted_commit':meta['accepted_commit'],
        'accepted_processed_nodes':94399,'unchanged_arithmetic_inputs':70,
        'accepted_open_boxes':14,'forged_prefix_variants_rejected':rejected,
        'zero_new_step_accounting_checked':True,'prior_arithmetic_replayed_by_this_test':False,
        'full_eight_point_exclusion_proved':False}


if __name__=='__main__':print(json.dumps(verify(),indent=2))
