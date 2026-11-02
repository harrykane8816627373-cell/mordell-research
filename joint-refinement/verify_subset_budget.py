"""Independent exact subset identities, strict cuts and preservation checks."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import gzip,json,tempfile
import subset_budget as D
from verify import C
from verify_force_enclosure import need,inside

E=D.E
HERE=Path(__file__).resolve().parent


def physical_point():
    base=[C(1),C(-1),C(0,1),C(0,-1),C(1,Q(1,2)),C(-1,-Q(1,2)),C(Q(1,2)),C(-Q(1,2))]
    raw=json.loads((HERE/'reference_centers.json').read_bytes())['centers'][0]
    x=[Q(round(Q(v,2**40)*1024),1024) for v in raw]
    target=[C(x[0])]+[C(x[2*i-1],x[2*i]) for i in range(1,7)]
    target.append(-sum(target,C()))
    delta=[a-b for a,b in zip(base,target)]
    alpha=2*sum((z.conjugate()*a).re for a,z in zip(base,delta))/sum(z.norm() for z in delta)
    q=[a-alpha*z for a,z in zip(base,delta)]
    need(sum(q,C())==C() and sum(z.norm() for z in q)==7,'independent exact centred mass-seven reflection')
    return q


def arithmetic():
    E.setup(40);settings=E.model('global_maximum');q=physical_point()
    variables=[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]]
    need(max(z.norm() for z in q)<Q(49,32),'retained point obeys the declared radius')
    for members,outside,coefficients in D.RECORDS:
        m=len(members)
        value=sum((q[i]-q[j]).norm() for index,i in enumerate(members) for j in members[index+1:])
        need(value==m*sum(q[i].norm() for i in members)-sum((q[i] for i in members),C()).norm(),
            'independent exact subset diameter identity')
        need(value>Q(m*m*(m-1),8),'the known point obeys every strict global subset inequality')
    retained=0
    for width in [0,1,E.S//256,E.S//32,E.S//8]:
        box=[E.rat(v.numerator,v.denominator) for v in variables]
        box=[(a-width-int(k==0),b+width) for k,(a,b) in enumerate(box)]
        result,reason,evidence=D.propagate(box,E.positions(box),settings)
        need(result is not None and all(inside(p,v) for p,v in zip(result,variables)),
            'all conditional norm and coordinate contractions retain a known feasible point')
        for p,z in zip(evidence['norms'],q):need(inside(p,z.norm()),'every propagated point norm retains the same actual point')
        retained+=1
    box=[E.rat(v.numerator,v.denominator) for v in variables]
    for k in [0,3,10]:box[k]=(box[k][0]-E.S//64,box[k][1]+E.S//64)
    identities=0;centroids=0
    for choices in product([0,1],repeat=3):
        values=[Q(a,E.S) for a,b in box]
        for k,side in zip([0,3,10],choices):values[k]=Q(box[k][side],E.S)
        actual=[C(values[0])]+[C(values[2*i-1],values[2*i]) for i in range(1,7)]
        actual.append(-sum(actual,C()))
        for members,outside,coefficients in D.RECORDS:
            total=sum((actual[i] for i in members),C());m=len(members)
            diameter=sum((actual[i]-actual[j]).norm() for index,i in enumerate(members) for j in members[index+1:])
            need(diameter==m*sum(actual[i].norm() for i in members)-total.norm(),'independent actual-corner quadratic identity')
            bound=D.centroid(box,E.positions(box),members,outside,coefficients)
            need(inside(bound[0],total.re) and inside(bound[1],total.im),'exact dependent-centroid interval includes its true corner')
            identities+=1;centroids+=1
    # Four small boxes: independent pair upper sums lose a factor compared with the shared norm budget.
    a=E.rat(2,5)[1];box=[E.rat(11,10)]+[(-a,a) for _ in range(8)]+[(-E.S,E.S) for _ in range(4)]
    qb=E.positions(box);members=(1,2,3,4);m=4
    naive=sum(E.cnorm(E.cs(qb[i],qb[j]))[1] for index,i in enumerate(members) for j in members[index+1:])
    shared=m*sum(E.cnorm(qb[i])[1] for i in members)
    threshold=6*E.S
    need(naive>threshold and shared<threshold,'strict shared-budget rejection where separate pair bounds are inconclusive')
    result,reason,_=D.propagate(box,qb,settings)
    need(result is None,'the entire continuous four-inner-particle box is rejected')
    # The strict threshold includes equality: a rational right triangle has diameter sum 9/4.
    triangle=[C(),C(Q(3,4)),C(Q(3,4),Q(3,4))]
    need(sum((triangle[i]-triangle[j]).norm() for i in range(3) for j in range(i+1,3))==Q(9,4),
        'exact strict-threshold equality example')
    r=settings['radius']
    boundary=[E.rat(11,10),E.rat(-3,4),E.rat(0),E.rat(0),E.rat(0),E.rat(0),E.rat(3,4)]+[(-r,r)]*6
    result,reason,evidence=D.propagate(boundary,E.positions(boundary),settings)
    need(result is None and reason=='subset_centroid_budget' and evidence['members']==[1,2,3]
        and evidence['upper']==evidence['strict_threshold_lower'],
        'the actual whole-box contractor rejects equality at the required strict subset threshold')
    try:D.propagate(box,qb,E.model('user_physical_system'))
    except ValueError:pass
    else:raise ArithmeticError('global-only subset bounds escaped their scope')
    return {'subset_classes':len(D.RECORDS),'independent_corner_identities':identities,
        'independent_dependent_centroid_enclosures':centroids,'known_mass7_boxes_retained':retained,
        'all_246_strict_subset_inequalities_verified_for_retained_point':True,
        'continuous_shared_budget_strictly_improves_naive_pair_bounds':True,
        'strict_threshold_rational_equality_checked':True,'global_scope_guard_checked':True,
        'floating_point_acceptance':False}


def prefix_guards():
    import recent_force as P
    from check_force import check
    base,accepted,meta=P.load_recent()
    baseline=json.loads(json.dumps(base));baseline['accepted_prefix_sha256']=P.SHA256
    rejected=0
    with tempfile.TemporaryDirectory(prefix='subset-prefix-guards-',dir=HERE.parent) as temp:
        path=Path(temp)/'test.json.gz'
        def save(data):path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        save(baseline);report=check(path,reuse_recent=True)
        need(report['counts']==accepted['counts'] and report['frontier']==accepted['frontier']
            and report['processed_since_accepted_prefix']==0 and report['closed_accepted_frontier_ids']==[]
            and report['remaining_accepted_frontier_regions']==30,'reused accepted prefix adds no fake progress')
        import subset_force
        valid=json.loads(json.dumps(baseline))
        i,box,root=next(row for row in valid['pending'] if row[0]==60819)
        result,reason,_=subset_force.contract(box,E.model('global_maximum'),passes=1)
        need(result is None,'independently computed one-step whole accepted-leaf closure')
        valid['nodes'][i]={'kind':'excluded','engine':'joint_force_subset_fast','reason':reason}
        valid['pending']=[row for row in valid['pending'] if row[0]!=i]
        save(valid);single=check(path,reuse_recent=True)
        need(single['processed_since_accepted_prefix']==1 and single['closed_accepted_frontier_ids']==[60819]
            and single['remaining_accepted_frontier_regions']==single['open_leaf_boxes']==29
            and single['newly_closed_since_accepted_prefix']==[],
            'an actual closed old leaf is counted separately from its still-open original parent')
        for kind in range(8):
            data=json.loads(json.dumps(baseline))
            if kind==0:data['accepted_prefix_sha256']='0'*64
            elif kind==1:
                i=next(i for i,n in enumerate(data['nodes']) if n['kind']=='excluded');data['nodes'][i]['reason']='forged'
            elif kind==2:data['roots'].pop()
            elif kind==3:data['pending'].pop()
            elif kind==4:data['pending'][0][1][0][0]+=1
            elif kind==5:data['nodes'].append({'kind':'unknown'})
            elif kind==6:data['full_eight_point_exclusion_proved']=True
            else:data['scope']='arbitrary_stable_configuration'
            save(data)
            try:check(path,reuse_recent=True)
            except ArithmeticError:rejected+=1
            else:raise ArithmeticError('forged recent-prefix proof accepted')
    need(rejected==8,'all forged new-prefix proofs rejected under Python -O')
    return {'accepted_processed_nodes':60791,'unchanged_arithmetic_inputs':61,
        'false_recent_prefix_variants_rejected':rejected,'old_arithmetic_replayed_by_this_test':False,
        'single_valid_leaf_closure_accounting_checked':True,
        'acceptance_uses_assert':False}


if __name__=='__main__':
    print(json.dumps({'passed':True,'subset_budget_arithmetic':arithmetic(),
        'accepted_prefix_guards':prefix_guards(),'full_eight_point_exclusion_proved':False},indent=2))
