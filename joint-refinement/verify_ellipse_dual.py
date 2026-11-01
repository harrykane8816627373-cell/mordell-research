"""Independent finite arithmetic checks for joint ellipse/mass-split bounds."""
from fractions import Fraction as F
import json
import ellipse_dual as J
import mass_dual as D
from verify import C,exact_log_product

E=J.E


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def independent_dual(lx,ly,a,b,lower,upper):
    return max(g-a*m for m,g in lx.items())+max(g-b*m for m,g in ly.items()) \
        +7*E.S**2*b+max(lower*(a-b),upper*(a-b))


def minimization_checks():
    E.setup(3)
    families=[({0:0,2:8,4:5},{1:2,3:10,5:1},F(3),F(4)),
        ({3:7,4:5},{3:3,4:4},F(16,5),F(19,5)),
        ({4:10},{4:10},F(7,2),F(7,2)),
        ({3:7,4:5},{3:3,4:4},F(3),F(3))]
    cases,points=0,0
    strict=0
    for left,right,lo,hi in families:
        lx={m*E.S**2:g*E.S**2 for m,g in left.items()}
        ly={m*E.S**2:g*E.S**2 for m,g in right.items()}
        lower,upper=lo*E.S**2,hi*E.S**2
        for tx,ty in [(-8,-3),(-1,5),(0,0),(4,12),(8,3),(16,16)]:
            result=J.minimize_split(lx,ly,tx,ty,lower,upper)
            # All activations of these independent input families lie above -100.
            expected=min(independent_dual(lx,ly,a,b,lower,upper)
                for a in range(-100,tx+1) for b in range(-100,ty+1))
            a,b=result['parameters']
            need(a<=tx and b<=ty and result['value']==expected
                and independent_dual(lx,ly,a,b,lower,upper)==expected,
                'two-parameter grid minimum equals independent exhaustive rectangular scan')
            old=D.minimize_dual(lx,ly,min(tx,ty))
            need(result['value']<=old['value'],'joint split bound dominates the common-parameter bound')
            strict+=result['value']<old['value']
            cases+=1
            points+=(tx+101)*(ty+101)
    need(strict>0,'there are strict gains beyond the common mass parameter')
    return {'independent_two_parameter_cases':cases,'integer_parameter_pairs_checked':points,
        'strict_improvement_cases':strict,'floating_point_acceptance':False}


def physical_witness():
    base=[C(F(2,3)),C(-F(2,3)),C(0,F(5,6)),C(0,-F(5,6)),
        C(1),C(-1),C(0,F(7,6)),C(0,-F(7,6))]
    ints=[(1123,0),(112,-1073),(427,-428),(112,1073),
        (-754,727),(-694,0),(427,428),(-754,-727)]
    target=[C(F(a,1000),F(b,1000)) for a,b in ints]
    mean=sum(target,C())/8
    target=[q-mean for q in target]
    direction=[a-b for a,b in zip(base,target)]
    numerator=sum((a.conjugate()*b for a,b in zip(direction,base)),C()).re
    denominator=sum(v.norm() for v in direction)
    q=[z-2*numerator/denominator*v for z,v in zip(base,direction)]
    return q


def physical_checks():
    E.setup(40)
    q=physical_witness()
    need(q[0].im==0 and sum(q,C())==C() and sum(z.norm() for z in q)==7 and len(set(q))==8,
        'actual thirteen-coordinate Gaussian-rational physical witness')
    s2=sum((z**2 for z in q),C())
    s3=sum((z**3 for z in q),C())
    ellipse=s2.norm()+F(54,175)*s3.norm()
    need(ellipse<1 and s2 and s3,'same actual point lies in the old nonzero-moment ellipse')
    values=[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]]
    center=[round(v*E.S) for v in values]
    mass_x=sum(z.re*z.re for z in q)
    lower_log,upper_log=exact_log_product(q)
    results=[]
    for divisor in [256,64,16]:
        radius=E.S//divisor
        box=[(v-radius,v+radius) for v in center]
        qbox=E.positions(box)
        for value,z in zip(q,qbox):
            need(all(F(interval[0],E.S)<=coord<=F(interval[1],E.S)
                for interval,coord in zip(z,[value.re,value.im])),'actual point is inside every tested continuous box')
        split=J.mass_split(qbox)
        need(split is not None and not split['empty']
            and split['lower']<=mass_x*E.S**2<=split['upper'],
            'joint mass split encloses exact actual component mass')
        new=J.target_upper(box,qbox,E.model()['radius'],center,split)
        old=D.target_upper(box,qbox,E.model()['radius'],center)
        need(new is not None and old is not None
            and upper_log<=F(new['log_Delta_over_target_upper'],E.S)
            and new['log_Delta_over_target_upper']<=old['log_Delta_over_target_upper'],
            'actual-point logarithm is below the joint whole-box bound, which dominates the old one')
        results.append({'radius':'1/'+str(divisor),
            'new_log_upper':str(F(new['log_Delta_over_target_upper'],E.S)),
            'old_log_upper':str(F(old['log_Delta_over_target_upper'],E.S)),
            'mass_split_lower':str(F(split['lower'],E.S**2)),
            'mass_split_upper':str(F(split['upper'],E.S**2))})
    return {'actual_center_mass_and_old_ellipse_exact':True,'both_low_moments_nonzero':True,
        'whole_boxes_checked':results,'witness_is_claimed_stationary':False,
        'floating_point_acceptance':False}


def ellipse_dependency():
    E.setup(40)
    q=[C(v) for v in [F(6,5),F(-6,5),F(9,10),F(-9,10),F(1),F(-1),F(1,2),F(-1,2)]]
    need(sum(q,C())==C() and sum(z.norm() for z in q)==7,'outside-ellipse comparison still has exact center and mass')
    box=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in q]
    split=J.mass_split(box)
    need(split is not None and split['empty'],'new split correctly rejects the ellipse constraint at this point')
    need(sum((z**2 for z in q),C()).norm()==49,'this rejection cannot be claimed for mass-only inputs')
    return {'old_ellipse_is_a_required_hypothesis':True,'outside_ellipse_mass_seven_example':True}


def coverage_guards():
    from pathlib import Path
    import gzip,tempfile
    from check_joint import accepted_parent,check,PARENT_SHA256
    parent,report,settings=accepted_parent()
    toy={'parent_sha256':PARENT_SHA256,'precision_bits':40,'settings':parent['settings'],
        'roots':[[i,j] for j,(i,b) in enumerate(parent['pending'])],
        'nodes':[{'kind':'unknown'} for _ in parent['pending']],
        'pending':[[j,b] for j,(i,b) in enumerate(parent['pending'])],
        'complete_cover':False,'full_eight_point_exclusion_proved':False}
    rejected=0
    with tempfile.TemporaryDirectory(prefix='joint-guard-',dir=Path(__file__).resolve().parent.parent) as temp:
        path=Path(temp)/'toy.gz'
        def save(data):
            path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        save(toy)
        need(check(path)['unresolved_regions']==19,'unchanged parent frontier remains fully visible')
        for kind in range(6):
            data=json.loads(json.dumps(toy))
            if kind==0:data['roots'].pop()
            elif kind==1:data['full_eight_point_exclusion_proved']=True
            else:
                data['pending'].pop(0)
                data['nodes'][0]={'kind':'excluded','reason':[
                    'ellipse_split_log','ellipse_mass_split','Hessian_radial_budget','Hessian_inner_particle'][kind-2]}
            save(data)
            try:check(path)
            except ArithmeticError:rejected+=1
            else:raise ArithmeticError('false joint coverage accepted')
    need(rejected==6,'all false joint certificates are rejected')
    return {'parent_open_regions':19,'false_proof_variants_rejected':rejected}


if __name__=='__main__':
    print(json.dumps({'passed':True,'joint_parameter_minimization':minimization_checks(),
        'actual_physical_boxes':physical_checks(),'scope_guards':ellipse_dependency(),
        'coverage_guards':coverage_guards(),
        'full_eight_point_exclusion_proved':False},indent=2))
