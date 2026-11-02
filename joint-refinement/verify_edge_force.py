"""Independent exact tests of correlated edges and arbitrary-reference force bounds."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import gzip,json,tempfile
import edge_force as S
import force_enclosure as F
import phase_force as P
from verify import C
from verify_force_enclosure import direct_force,need,inside

E=F.E


def checks():
    E.setup(40)
    reference=[C(1),C(-1),C(0,1),C(0,-1),C(1,1),C(-1,-1),C(2),C(-2)]
    variables=[reference[0].re]+[v for z in reference[1:7] for v in [z.re,z.im]]
    fixed=[int(v*E.S) for v in variables]
    box=[(v,v) for v in fixed]
    for k in [0,2,5]:box[k]=(fixed[k]-E.S//64-int(k==0),fixed[k]+E.S//64)
    qbox=E.positions(box)
    old=F.enclose(box,qbox,endpoint_separation=E.S//2)
    sharp=S.enclose(box,qbox,E.S//2)
    need(old is not None and sharp is not None,'both continuous midpoint enclosures exist')
    need(old['value']==sharp['value'] and old['jacobian']==sharp['jacobian']
        and old['affine']==sharp['affine'],'only the edge perturbation bounds change at a fixed reference')
    need(all(a<=b for a,b in zip(sharp['remainder'],old['remainder'])),
        'each sharper quadratic remainder dominates the prior point-radius sum bound')
    strict=sum(r['upper']<r['previous_upper'] for r in sharp['edge_bounds'])
    need(strict>0 and any(a<b for a,b in zip(sharp['remainder'],old['remainder'])),
        'a nonempty correlated continuous box exhibits strict improvement')
    phase=P.enclose(box,qbox,E.S//2)
    need(phase is not None,'complex remainder exists on the same continuous box')
    phase_strict=sum(v[1]-v[0]<2*r for z,r in zip(phase['complex_remainder'],phase['remainder']) for v in z)
    need(phase_strict>0,'complex remainder strictly improves at least one real-component interval')
    corner_count=0;edge_count=0;newton_count=0;phase_count=0;phase_newton=0
    for choices in product([0,1],repeat=3):
        values=list(variables)
        for k,side in zip([0,2,5],choices):values[k]=Q(box[k][side],E.S)
        points=[C(values[0])]+[C(values[2*i-1],values[2*i]) for i in range(1,7)]
        points.append(-sum(points,C()))
        need(min((points[i]-points[j]).norm() for i,j in E.PAIRS)>Q(1,2),
            'independent corner input meets the declared physical endpoint separation')
        actual=direct_force(points)
        target=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in actual]
        for center in [None,[2*v for v in sharp['center2']]]:
            data=S.enclose(box,qbox,E.S//2,center)
            phase=P.enclose(box,qbox,E.S//2,center)
            need(data is not None,'arbitrary distinct reference exists')
            need(phase is not None,'arbitrary-reference complex remainder exists')
            cv=[Q(v,2*E.S) for v in data['center2']]
            cq=[C(cv[0])]+[C(cv[2*i-1],cv[2*i]) for i in range(1,7)]
            cq.append(-sum(cq,C()))
            fc=direct_force(cq)
            for record,(i,j) in zip(data['edge_bounds'],E.PAIRS):
                h=(points[i]-cq[i])-(points[j]-cq[j])
                need(h.norm()<=Q(record['upper'],2*E.S)**2,'independent exact edge perturbation lies inside its correlated bound')
                edge_count+=1
            for i in range(8):
                derivative=[]
                for k in range(13):
                    h=[C(*E.L[j][k]) for j in range(8)]
                    derivative.append(-sum(((h[i]-h[j])/(cq[i]-cq[j])**2 for j in range(8) if i!=j),C())-4*h[i].conjugate())
                remainder=actual[i]-fc[i]-sum((v*(a-b) for v,a,b in zip(derivative,values,cv)),C())
                need(remainder.norm()<=Q(data['remainder'][i],E.S)**2,'independent true force remainder at this reference')
                for v,value in zip(phase['complex_remainder'][i],[remainder.re,remainder.im]):
                    need(inside(v,value),'independent exact real and imaginary remainder components are enclosed')
                    phase_count+=1
                for c,value in enumerate([actual[i].re,actual[i].im]):
                    lo,hi=data['affine'][i][c];r=data['remainder'][i]
                    need(Q(lo-r,E.S)<=value<=Q(hi+r,E.S),'actual force is enclosed for the entire original box')
            image=F.newton_image(data,target)
            need(image is not None and all(inside(p,v) for p,v in zip(image,values)),
                'arbitrary-reference Newton necessary image retains a known shifted-force solution')
            newton_count+=1
            image=P.newton_image(phase,target)
            need(image is not None and all(inside(p,v) for p,v in zip(image,values)),
                'complex-remainder Newton necessary image retains an exact shifted-force solution')
            phase_newton+=1
        corner_count+=1
    pole=[0]*13
    need(S.enclose(box,qbox,E.S//2,pole) is None,'reference collisions produce no bound')
    return {'independent_actual_corners':corner_count,'independent_exact_edge_records':edge_count,
        'known_shifted_force_Newton_images':newton_count,'strictly_improved_edges_at_midpoint':strict,
        'independent_complex_remainder_components':phase_count,
        'known_complex_remainder_Newton_images':phase_newton,
        'strictly_improved_remainder_components_at_midpoint':phase_strict,
        'arbitrary_reference_outside_box_checked':True,'odd_width_half_grid_center_checked':True,
        'previous_midpoint_bound_dominated':True,'floating_point_acceptance':False}


def checkpoint_guards():
    import checkpoint_force as P
    from check_force import check
    base,accepted,meta=P.load_checkpoint()
    baseline=json.loads(json.dumps(base));baseline['checkpoint_sha256']=P.SHA256
    rejected=0
    with tempfile.TemporaryDirectory(prefix='edge-proof-guards-',dir=P.ROOT) as temp:
        path=Path(temp)/'test.json.gz'
        def save(data):path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        save(baseline)
        result=check(path,reuse_checkpoint=True)
        need(result['counts']==accepted['counts'] and result['frontier']==accepted['frontier']
            and result['processed_since_checkpoint']==0 and result['newly_closed_since_checkpoint']==[],
            'reused authentic checkpoint has precisely its accepted forest and no new progress')
        for kind in range(9):
            data=json.loads(json.dumps(baseline))
            if kind==0:data['checkpoint_sha256']='0'*64
            elif kind==1:
                i=next(i for i,n in enumerate(data['nodes']) if n['kind']=='excluded')
                data['nodes'][i]['reason']='forged'
            elif kind==2:data['roots'].pop()
            elif kind==3:data['nodes'].pop()
            elif kind==4:data['pending'][0][1][0][0]+=1
            elif kind==5:data['pending'].pop()
            elif kind==6:data['full_eight_point_exclusion_proved']=True
            elif kind==7:data['nodes'].append({'kind':'unknown'})
            else:data['scope']='arbitrary_stable_configuration'
            save(data)
            try:check(path,reuse_checkpoint=True)
            except ArithmeticError:rejected+=1
            else:raise ArithmeticError('forged reused-checkpoint forest accepted')
        isolated=Path(temp)/'isolated'
        for row in meta['arithmetic_inputs']:
            destination=isolated/row['path'];destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes((P.ROOT/row['path']).read_bytes())
        for name in ['force_checkpoint_768598a.json','force_checkpoint_768598a.json.gz','force_checkpoint_768598a_report.json']:
            destination=isolated/'joint-refinement/certificates'/name
            destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes((P.HERE/'certificates'/name).read_bytes())
        old_root,old_here=P.ROOT,P.HERE
        try:
            P.ROOT,P.HERE=isolated,isolated/'joint-refinement'
            P.load_checkpoint()
            for name in ['joint-refinement/force_enclosure.py','joint-refinement/certificates/force_checkpoint_768598a.json']:
                target=isolated/name;original=target.read_bytes();target.write_bytes(original+b'\n')
                try:P.load_checkpoint()
                except ArithmeticError:rejected+=1
                else:raise ArithmeticError('changed checkpoint arithmetic or metadata accepted')
                target.write_bytes(original)
        finally:P.ROOT,P.HERE=old_root,old_here
    need(rejected==11,'all checkpoint lineage, pending coverage, source and metadata forgeries rejected')
    return {'accepted_processed_nodes':25790,'unchanged_arithmetic_inputs':50,
        'unchanged_checkpoint_has_no_new_progress':True,'false_checkpoint_variants_rejected':rejected,
        'acceptance_uses_assert':False,'old_arithmetic_replayed_by_this_test':False}


def gram_checks():
    import pair_gram as G
    import force_contract as K
    base=[C(1),C(-1),C(0,1),C(0,-1),C(1,Q(1,2)),C(-1,-Q(1,2)),C(Q(1,2)),C(-Q(1,2))]
    pair_count=0;mode_count=0;strict_det=None
    for scale in [Q(1),Q(45,32),Q(100)]:
        q=[scale*z for z in base]
        box=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in q]
        inverse,_=E.invert_edges(box)
        evidence=G.exclusion(box,inverse)
        b=[[C() for _ in range(8)] for _ in range(8)]
        for i,j in E.PAIRS:
            w=1/(q[i]-q[j])**2
            b[i][i]+=w;b[j][j]+=w;b[i][j]=-w;b[j][i]=-w
        d=[[C(14 if i==j else -2)-sum((b[k][i].conjugate()*b[k][j] for k in range(8)),C())
            for j in range(8)] for i in range(8)]
        need(all(d[i][i].im==0 for i in range(8)),'independent exact Gram diagonals are real')
        determinants={(i,j):d[i][i].re*d[j][j].re-d[i][j].norm() for i,j in E.PAIRS}
        pair_count+=len(determinants)
        if evidence:
            if evidence['reason']=='Hessian_Gram_diagonal':need(d[evidence['point']][evidence['point']].re<0,'independent true negative Gram diagonal')
            else:need(determinants[tuple(evidence['pair'])]<0,'independent true negative complex principal determinant')
        if scale==Q(45,32):
            need(evidence is not None and evidence['reason']=='Hessian_Gram_pair'
                and all(d[i][i].re>0 for i in range(8)),'strictly stronger than the individual column-norm budgets')
            need(K.subset_exclusion(box,inverse) is None,'the existing interval indicator probes are inconclusive here')
            for mask in range(1,128):
                action=[sum((b[i][j] for j in range(8) if mask&(1<<j)),C()) for i in range(8)]
                k=mask.bit_count()
                need(sum(z.norm() for z in action)<=2*k*(8-k),'all 127 exact old indicator inequalities hold')
                mode_count+=1
            strict_det=str(determinants[tuple(evidence['pair'])])
        if scale==100:
            row_bound=2*max(sum((Q(1,(q[i]-q[j]).norm()) for j in range(8) if i!=j),Q()) for i in range(8))
            need(row_bound<4 and evidence is None,'a rigorous strictly positive full-Hessian configuration survives the pair tests')
    return {'independent_complex_principal_determinants':pair_count,
        'strict_improvement_example_scale':'45/32','strict_improvement_example_mass':'14175/1024',
        'strict_improvement_negative_determinant':strict_det,
        'old_exact_indicator_inequalities_satisfied':mode_count,
        'all_single_column_budgets_satisfied_in_strict_example':True,
        'wide_positive_Hessian_input_retained':True,'all_two_column_tests_are_only_necessary':True,
        'floating_point_acceptance':False}


if __name__=='__main__':
    print(json.dumps({'passed':True,'correlated_edge_arithmetic':checks(),
        'full_Hessian_pair_Gram':gram_checks(),
        'checkpoint_reuse_guards':checkpoint_guards(),
        'full_eight_point_exclusion_proved':False},indent=2))
