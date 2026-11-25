"""Replay the sharp ten-point bound and its equality certificate."""

import argparse
import json

import antipodal
import cover

model=cover.model


def verification():
    contacts,kernel=model.exact_kernel()
    limit=100000
    inner=cover.inner_cover(kernel,contacts,limit)
    names={1:'mixed',2:'outer_step_one',3:'outer_step_two',4:'outer_step_three'}
    general={name:cover.cover_class(cls,kernel,contacts,limit) for cls,name in names.items()}
    opposite=antipodal.verify(kernel,contacts,limit)
    covers={'inner':inner,**general,'outer_antipodal':opposite}
    model.need([result['class'] for result in covers.values()]==list(range(6)) and
               all(result['unresolved']==0 for result in covers.values()),'all six physical domains covered')
    model.need(contacts['pair_counts_by_class']==[1,16,8,8,8,4] and
               len(contacts['local_boxes'])==15 and contacts['positive_gram_blocks']==12 and
               contacts['positive_second_variation_blocks']==2,'exact contacts and equality blocks')
    counts={key:sum(result[piece]['visited'] for piece in ('normal','face')) if key=='inner'
            else sum(root['visited'] for root in result['roots']) if key in general else result['visited']
            for key,result in covers.items()}
    report={'passed':True,'n':10,'normalization_mass':9,'pair_classes':len(covers),
            'contact_classes':len(contacts['local_boxes']),'coverage_nodes_by_class':counts,
            'total_coverage_nodes':sum(counts.values()),'compensation':'0','unresolved':0,
            'sharp_global_value_proved':True,'equality':'unique up to rotation and point permutation'}
    files={'verification_report.json':contacts,
           **{name+'_cover_report.json':result for name,result in covers.items()}}
    return report,files


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    report,files=verification()
    files['global_verification_report.json']=report
    for name,result in files.items():
        path=model.HERE/name
        if args.record_report:
            path.write_bytes((json.dumps(result,indent=2)+'\n').encode())
        else:
            model.need(json.loads(path.read_bytes())==result,'matches frozen report '+name)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
