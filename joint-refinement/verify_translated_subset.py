"""Independent rational checks of the translated continuous-subset bound."""
from fractions import Fraction as Q
from itertools import product
import json
import translated_subset as T
from verify import C
from verify_force_enclosure import need,inside
from verify_subset_budget import physical_point

E=T.E


def main():
    E.setup(40);settings=E.model('global_maximum');q=physical_point()
    variables=[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]]
    retained=0
    for width in [0,1,E.S//256,E.S//32,E.S//8]:
        box=[E.rat(v.numerator,v.denominator) for v in variables]
        box=[(a-width-int(k==0),b+width) for k,(a,b) in enumerate(box)]
        qb=E.positions(box)
        need(T.exclusion(box,qb,settings) is None,'translated cuts retain the same exact mass-seven point')
        for members,outside,coefficients in T.D.RECORDS:
            m=len(members)
            actual=sum((q[i]-q[j]).norm() for k,i in enumerate(members) for j in members[k+1:])
            need(actual>Q(m*m*(m-1),8),'retained point obeys every strict subset inequality')
            evidence=T.bound(box,qb,members,outside,coefficients)
            need(actual<=Q(evidence['upper'],E.S),'outward bound contains the actual retained point')
        retained+=1
    need(sum(q,C())==C() and sum(z.norm() for z in q)==7,'known point is genuinely centred and mass seven')
    box=[E.rat(v.numerator,v.denominator) for v in variables]
    for k in [0,3,10]:box[k]=(box[k][0]-E.S//64,box[k][1]+E.S//64)
    qb=E.positions(box);identities=0;enclosures=0
    for choices in product([0,1],repeat=3):
        values=[Q(a,E.S) for a,b in box]
        for k,side in zip([0,3,10],choices):values[k]=Q(box[k][side],E.S)
        actual=[C(values[0])]+[C(values[2*i-1],values[2*i]) for i in range(1,7)]
        actual.append(-sum(actual,C()))
        for members,outside,coefficients in T.D.RECORDS:
            evidence=T.bound(box,qb,members,outside,coefficients);m=len(members)
            center=C(*(Q(v,E.S) for v in evidence['center']))
            total=sum((actual[i] for i in members),C())
            diameter=sum((actual[i]-actual[j]).norm() for k,i in enumerate(members) for j in members[k+1:])
            translated=m*sum((actual[i]-center).norm() for i in members)-(total-m*center).norm()
            need(diameter==translated,'independent Gaussian-rational translated variance identity')
            need(translated<=Q(evidence['upper'],E.S),'every actual corner is below the outward upper bound')
            sum_box=T.D.centroid(box,qb,members,outside,coefficients)
            need(inside(sum_box[0],total.re) and inside(sum_box[1],total.im),'dependent centroid enclosure retains its actual sum')
            identities+=1;enclosures+=1
    # Three continuous boxes with centre -1/2-i/2 and half-width 1/3.
    # Independent pair sums allow 8/3; this fixed-centre variance bound gives 2.
    radius=settings['radius'];box=[E.rat(11,10)]+[(-radius,radius)]*12
    for i in [2,3,4]:
        box[2*i-1]=(E.rat(-5,6)[0],E.rat(-1,6)[1])
        box[2*i]=(E.rat(-5,6)[0],E.rat(-1,6)[1])
    qb=E.positions(box)
    record=next(row for row in T.D.RECORDS if row[0]==(2,3,4))
    evidence=T.bound(box,qb,*record);threshold=evidence['strict_threshold_lower'];members=record[0]
    naive=sum(E.cnorm(E.cs(qb[i],qb[j]))[1] for k,i in enumerate(members) for j in members[k+1:])
    center=T.D.centroid(box,qb,*record)
    zero=3*sum(E.cnorm(qb[i])[1] for i in members)-E.cnorm(center)[0]
    need(naive>threshold and zero>threshold and evidence['upper']<threshold,
        'entire continuous box improves both separate-pair and zero-centre variance bounds')
    rejection=T.exclusion(box,qb,settings)
    need(rejection is not None and rejection['reason']=='translated_subset_budget','new whole-box criterion rejects the continuous example')
    # A dyadic right triangle sits exactly at the strict 9/4 threshold.
    boundary=[E.rat(11,10),E.rat(-3,4),E.rat(0),E.rat(0),E.rat(0),E.rat(0),E.rat(3,4)]+[(-radius,radius)]*6
    qb=E.positions(boundary);record=next(row for row in T.D.RECORDS if row[0]==(1,2,3))
    evidence=T.bound(boundary,qb,*record)
    need(evidence['upper']==evidence['strict_threshold_lower']==9*E.S//4,'exact boundary equality is preserved by all outward operations')
    rejection=T.exclusion(boundary,qb,settings)
    need(rejection is not None and rejection['reason']=='translated_subset_budget','strict threshold also rejects equality')
    inconsistent=[E.rat(1)]+[E.ZERO]*12
    qb_bad=E.positions(inconsistent);qb_bad[7]=(E.rat(1),E.ZERO)
    rejection=T.exclusion(inconsistent,qb_bad,settings)
    need(rejection is not None and rejection['reason']=='translated_subset_geometry',
        'incompatible direct and exactly centred subset sums are rejected as an empty joint enclosure')
    rejected=0
    for call in [lambda:T.exclusion(boundary,qb,E.model('user_physical_system')),
        lambda:T.contract(boundary,E.model('user_physical_system'))]:
        try:call()
        except ValueError:rejected+=1
        else:raise ArithmeticError('global subset theorem escaped its scope')
    print(json.dumps({'passed':True,'subset_classes':len(T.D.RECORDS),
        'independent_translated_corner_identities':identities,'dependent_centroid_enclosures':enclosures,
        'exact_mass7_boxes_retained':retained,'continuous_strict_strengthening_checked':True,
        'strict_threshold_equality_rejected':True,'empty_centred_intersection_rejected':True,
        'scope_guard_rejections':rejected,
        'floating_point_acceptance':False,'acceptance_uses_assert':False,
        'full_eight_point_exclusion_proved':False},indent=2))


if __name__=='__main__':main()
