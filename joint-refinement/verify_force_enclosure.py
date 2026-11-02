"""Independent exact force/Taylor/Newton, imported inner-radius and spectrum checks."""
from fractions import Fraction as Q
from pathlib import Path
import gzip,json
import force_enclosure as F
import force_contract as C7
from verify import C
from verify_ellipse_dual import physical_witness
from verify_radial_stability import add,multiply,scale

E=F.E


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def inside(interval,x):return Q(interval[0],E.S)<=x<=Q(interval[1],E.S)


def direct_force(points):
    return [sum((C(1)/(q-z) for j,z in enumerate(points) if i!=j),C())-4*q.conjugate() for i,q in enumerate(points)]


def force_checks():
    E.setup(40)
    points=physical_witness()
    variables=[points[0].re]+[v for z in points[1:7] for v in [z.re,z.im]]
    center=[round(v*E.S) for v in variables]
    count=0
    actual=direct_force(points)
    for divisor in [1024,256,128]:
        radius=E.S//divisor
        box=[(c-radius-int(k%2==0),c+radius) for k,c in enumerate(center)]
        enclosure=F.enclose(box)
        need(enclosure is not None,'strict pole clearance in independent actual-point boxes')
        cv=[Q(v,2*E.S) for v in enclosure['center2']]
        reference=[C(cv[0])]+[C(cv[2*i-1],cv[2*i]) for i in range(1,7)]
        reference.append(-sum(reference,C()))
        forces=direct_force(reference)
        delta=[a-b for a,b in zip(variables,cv)]
        for i in range(8):
            need(inside(enclosure['value'][i][0],forces[i].re) and inside(enclosure['value'][i][1],forces[i].im),
                'exact half-grid center residual is independently enclosed')
            jac=[]
            for k in range(13):
                direction=[C(*E.L[j][k]) for j in range(8)]
                derivative=-sum(((direction[i]-direction[j])/(reference[i]-reference[j])**2
                    for j in range(8) if i!=j),C())-4*direction[i].conjugate()
                need(inside(enclosure['jacobian'][i][k][0],derivative.re)
                    and inside(enclosure['jacobian'][i][k][1],derivative.im),'independent rational center derivative')
                jac.append(derivative)
                count+=1
            residual=actual[i]-forces[i]-sum((j*d for j,d in zip(jac,delta)),C())
            need(residual.norm()<=Q(enclosure['remainder'][i],E.S)**2,'exact actual quadratic force remainder lies in the disk')
            for c,value in enumerate([actual[i].re,actual[i].im]):
                a,b=enclosure['affine'][i][c]
                r=enclosure['remainder'][i]
                need(Q(a-r,E.S)<=value<=Q(b+r,E.S),'actual full force is inside whole-box affine enclosure')
        target=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in actual]
        image=F.newton_image(enclosure,target)
        need(image is not None and all(inside(interval,value) for interval,value in zip(image,variables)),
            'the exact actual point survives a shifted-force Newton map with known target force')
    p=F.preconditioner([[E.rat(2*int(i==j)) for j in range(13)] for i in range(13)])
    need(p==[[E.S//2 if i==j else 0 for j in range(13)] for i in range(13)],'Newton inverse scaling is independently exact')
    need(F.preconditioner([[E.ZERO]*13 for _ in range(13)]) is None,'singular midpoint returns unknown')
    need(F.enclose([(-E.S,E.S)]*13) is None,'a pole-crossing disk never becomes an exclusion')
    return {'actual_continuous_boxes':3,'independent_complex_Jacobian_entries':count,
        'odd_width_half_grid_centers':True,'exact_actual_remainders_checked':24,
        'shifted_force_Newton_images_checked':3,'floating_point_acceptance':False}


def inherited_threshold():
    a=Q(729,1000)
    slope=Q(1083,20)
    f=[Q(0),Q(0),Q(-162),Q(288),Q(-120)]
    fa=sum(v*a**k for k,v in enumerate(f))
    constant=fa-slope*a
    lower=add(add(scale(multiply([0,1],multiply([Q(-1073,1000),1],[Q(-1073,1000),1])),120),
        scale(multiply([Q(-63718,59375),1],[Q(-63718,59375),1]),57)),[Q(474140217,11875000000)])
    identity=add([constant,slope],scale(f,-1))
    need(identity==multiply([-a,1],lower),'imported stronger radius polynomial identity checked by all coefficients')
    gap=7*slope+8*constant+4
    need(gap==-Q(49638643,3125000000)<0,'strict imported mass/ellipse contradiction at radius 729/1000')
    eighth_upper=(7-7*a)**4+7*a**4
    spectral_lower=16-2*gap/eighth_upper
    need(eighth_upper==Q(1865874186281,125000000000)
        and spectral_lower==Q(29857958071936,1865874186281)
        and spectral_lower>Q(16001,4000)**2,
        'the imported annular identity also gives a strictly negative full-Hessian gap below -1/4000')
    return {'source':'user attachment, not a newly attributed discovery','inner_squared_radius_threshold':'729/1000',
        'exact_contradiction_margin':str(gap),'coefficient_identity_checked':True,
        'additional_uniform_Hessian_min_strict_upper':'-1/4000',
        'curvature_gap_source':'derived here from the uploaded polynomial identity'}


def endpoint_checks():
    E.setup(40)
    reference=[C(1),C(-1),C(0,1),C(0,-1),C(1,1),C(-1,-1),C(2),C(-2)]
    points=[Q(11,10)*z for z in reference]
    values=[z.re for z in reference[:1]]+[v for z in reference[1:7] for v in [z.re,z.im]]
    box=[(int(v*E.S)-E.S//2,int(v*E.S)+E.S//2) for v in values]
    need(F.enclose(box) is None,'the original local disks cross their reference pole-clearance limit')
    separation=min((points[i]-points[j]).norm() for i,j in E.PAIRS)
    need(separation>Q(1,2),'the actual target point independently satisfies the inherited endpoint separation')
    data=F.enclose(box,endpoint_separation=E.S//2)
    need(data is not None,'known endpoint separation supplies a valid global finite-identity enclosure')
    center_force=direct_force(reference)
    actual=direct_force(points)
    variables=[points[0].re]+[v for z in points[1:7] for v in [z.re,z.im]]
    for i in range(8):
        derivative=[]
        for k in range(13):
            h=[C(*E.L[j][k]) for j in range(8)]
            derivative.append(-sum(((h[i]-h[j])/(reference[i]-reference[j])**2
                for j in range(8) if i!=j),C())-4*h[i].conjugate())
        delta=[a-b for a,b in zip(variables,values)]
        remainder=actual[i]-center_force[i]-sum((v*d for v,d in zip(derivative,delta)),C())
        need(remainder.norm()<=Q(data['remainder'][i],E.S)**2,'independent actual global endpoint remainder')
        for c,value in enumerate([actual[i].re,actual[i].im]):
            lo,hi=data['affine'][i][c];r=data['remainder'][i]
            need(Q(lo-r,E.S)<=value<=Q(hi+r,E.S),'the actual separated point survives the conditional global force enclosure')
    target=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in actual]
    image=F.newton_image(data,target)
    need(image is not None and all(inside(interval,value) for interval,value in zip(image,variables)),
        'the conditional endpoint Newton image retains the exact shifted-force solution')
    return {'local_disk_method_is_inconclusive':True,'actual_endpoint_separation_verified':True,
        'global_actual_remainders_checked':8,'global_shifted_force_Newton_checked':True,
        'mass_seven_is_not_required_for_this_force_lemma':True,'floating_point_acceptance':False}


def spectral_checks():
    E.setup(40)
    q=[C(1),C(-1),C(0,1),C(0,-1),C(1,Q(1,2)),C(-1,-Q(1,2)),C(Q(1,2)),C(-Q(1,2))]
    records=0
    for scale_factor in [1,100]:
        points=[scale_factor*z for z in q]
        box=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in points]
        inv,_=E.invert_edges(box)
        evidence=C7.subset_exclusion(box,inv)
        violated=[]
        for mask in range(1,128):
            direction=[C(int(bool(mask&(1<<i)))) for i in range(8)]
            action=[sum(((direction[i]-direction[j])/(points[i]-points[j])**2 for j in range(8) if i!=j),C()) for i in range(8)]
            lhs=sum(z.norm() for z in action)
            k=mask.bit_count()
            if lhs>2*k*(8-k):violated.append(mask)
            records+=1
        need((evidence is not None)==bool(violated),'interval subset rejection agrees with exact actual operator actions')
        if evidence:need(evidence['mask'] in violated,'the returned integer subset is an independently checked witness')
        if scale_factor==100:
            row_bound=2*max(sum((Q(1,(points[i]-points[j]).norm()) for j in range(8) if i!=j),Q()) for i in range(8))
            need(row_bound<4 and not violated,'a rigorously positive-Hessian widely spaced configuration is not rejected')
    return {'independent_exact_indicator_modes':records,'all_complement_classes_checked':127,
        'wide_positive_Hessian_input_retained':True,'floating_point_acceptance':False}


def proof_guards():
    import tempfile
    from propose_force import initialize
    from check_force import check
    base=initialize()
    rejected=0
    with tempfile.TemporaryDirectory(prefix='force-guards-',dir=Path(__file__).resolve().parent.parent) as temp:
        path=Path(temp)/'test.json.gz'
        def save(data):path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        save(base)
        result=check(path)
        need(result['total_closed_parent_regions']==6 and result['remaining_parent_regions']==15,
            'the initial composed upload proves exactly its six original closures')
        for kind in range(8):
            data=json.loads(json.dumps(base))
            if kind==0:data['full_eight_point_exclusion_proved']=True
            elif kind==1:data['roots'].pop()
            elif kind==2:data['pending'].pop()
            elif kind==3:
                node=next(n for n in data['nodes'] if n.get('engine')=='uploaded' and n['kind']=='excluded')
                node['reason']='forged'
            elif kind in [4,5]:
                data['pending']=[r for r in data['pending'] if r[0]!=0]
                data['nodes'][0]={'kind':'excluded','engine':'joint_force','reason':
                    'force_affine_exact' if kind==4 else 'force_quadratic_Newton'}
            elif kind==6:data['nodes'].append({'kind':'unknown'})
            else:
                data['pending']=[r for r in data['pending'] if r[0]!=0]
                data['nodes'][0]={'kind':'excluded','engine':'joint_force_global','reason':'force_affine_exact'}
            save(data)
            try:check(path)
            except ArithmeticError:rejected+=1
            else:raise ArithmeticError('false force proof accepted')
    need(rejected==8,'all forged force proof variants rejected even with Python -O')
    return {'original_six_closures_preserved':True,'false_proof_variants_rejected':rejected,
        'acceptance_uses_assert':False}


if __name__=='__main__':
    print(json.dumps({'passed':True,'affine_and_quadratic_Newton':force_checks(),
        'global_endpoint_force':endpoint_checks(),
        'uploaded_radial_threshold':inherited_threshold(),'full_Hessian_subset_modes':spectral_checks(),
        'proof_guards':proof_guards(),
        'full_eight_point_exclusion_proved':False},indent=2))
