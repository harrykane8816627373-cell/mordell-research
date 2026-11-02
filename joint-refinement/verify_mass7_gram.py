"""Independently check strict probe improvement within mass seven and the old ellipse.

This illustrative witness is not used as a leaf in the global cover. Its
stationarity and target-value failures remain explicit.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json
from verify import C
from verify_force_enclosure import direct_force,need
import pair_gram as G
import force_contract as K

HERE=Path(__file__).resolve().parent
E=G.E


def check():
    data=json.loads((HERE/'certificates/mass7_gram_witness.json').read_bytes())
    need(data['construction']=='householder_mass7_v1' and data['target_grid_denominator']==1024,
        'declared exact centred rational reflection construction')
    coords=data['target_coordinates']
    need(len(coords)==13 and all(type(v) is int for v in coords),'thirteen integer target coordinates')
    x=[Q(v,1024) for v in coords]
    target=[C(x[0])]+[C(x[2*i-1],x[2*i]) for i in range(1,7)]
    target.append(-sum(target,C()))
    base=[C(1),C(-1),C(0,1),C(0,-1),C(1,Q(1,2)),C(-1,-Q(1,2)),C(Q(1,2)),C(-Q(1,2))]
    delta=[a-b for a,b in zip(base,target)]
    norm=sum(z.norm() for z in delta)
    need(norm>0,'nonzero reflection direction')
    alpha=2*sum((z.conjugate()*a).re for a,z in zip(base,delta))/norm
    need(alpha==Q(data['reflection_coefficient']),'exact declared reflection coefficient')
    q=[a-alpha*z for a,z in zip(base,delta)]
    need(sum(q,C())==C() and sum(z.norm() for z in q)==7,'actual exact centred mass-seven points')
    distances=[(q[i]-q[j]).norm() for i,j in E.PAIRS]
    need(min(distances)>Q(1,2) and max(z.norm() for z in q)<Q(49,32),
        'the inherited global endpoint and radius conditions also hold')
    s2=sum((z*z for z in q),C());s3=sum((z*z*z for z in q),C())
    ellipse=s2.norm()+Q(54,175)*s3.norm()
    need(s2.norm()>0 and s3.norm()>0 and ellipse<Q(1,5),'strict old ellipse and both moments nonzero')
    deficit=sum(z.norm()*(z.norm()-1)**2 for z in q)
    need(deficit>Q(1,1200),'outside the already excluded radial region')
    b=[[C() for _ in range(8)] for _ in range(8)]
    for i,j in E.PAIRS:
        w=1/(q[i]-q[j])**2
        b[i][i]+=w;b[j][j]+=w;b[i][j]=-w;b[j][i]=-w
    d=[[C(14 if i==j else -2)-sum((b[k][i].conjugate()*b[k][j] for k in range(8)),C())
        for j in range(8)] for i in range(8)]
    need(all(d[i][i].im==0 and d[i][i].re>0 for i in range(8)),'all eight exact single-column budgets hold strictly')
    for mask in range(1,128):
        k=mask.bit_count()
        action=[sum((b[i][j] for j in range(8) if mask&(1<<j)),C()) for i in range(8)]
        need(sum(z.norm() for z in action)<2*k*(8-k),'every exact old indicator inequality holds strictly')
    i,j=data['negative_principal_pair']
    determinant=d[i][i].re*d[j][j].re-d[i][j].norm()
    need(Q(-7,500)<determinant<Q(-1,100),'independent strictly negative complex Gram determinant')
    E.setup(40)
    box=[(E.rat(z.re.numerator,z.re.denominator),E.rat(z.im.numerator,z.im.denominator)) for z in q]
    inverse,_=E.invert_edges(box)
    need(K.subset_exclusion(box,inverse) is None,'old interval probes retain this actual input')
    evidence=G.exclusion(box,inverse)
    need(evidence is not None and evidence['reason']=='Hessian_Gram_pair'
        and evidence['pair']==[i,j],'new integer interval implementation rejects the same true principal minor')
    product=Q(1)
    for distance in distances:product*=distance
    need(product<7**7 and sum(z.norm() for z in direct_force(q))>0,
        'this is explicitly not a stationary or target-threshold counterexample')
    return {'passed':True,'mass':7,'centred':True,'reflection_coefficient':str(alpha),
        'old_ellipse_strict_upper':'1/5','both_old_moments_nonzero':True,
        'squared_separation_strict_lower':'1/2','squared_radius_strict_upper':'49/32',
        'radial_deficit_strict_lower':'1/1200','single_column_budgets_strictly_satisfied':8,
        'indicator_inequalities_strictly_satisfied':127,'negative_principal_pair':[i,j],
        'negative_determinant_strict_interval':['-7/500','-1/100'],
        'old_interval_probes_retain_input':True,'new_interval_pair_Gram_rejects_input':True,
        'stationarity_satisfied':False,'target_threshold_satisfied':False,
        'illustrative_probe_separation_only':True,'floating_point_acceptance':False,
        'full_eight_point_exclusion_proved':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args();result=check();path=HERE/'reports/mass7_gram.json'
    if args.record_report:path.write_bytes((json.dumps(result,indent=2)+'\n').encode())
    else:need(result==json.loads(path.read_bytes()),'deterministic exact mass-seven witness report')
    print(json.dumps(result,indent=2))
