"""Global-maximizer subset inequalities coupled to one exact mass and centroid."""
from itertools import combinations
import force_enclosure as F

E=F.E
RECORDS=[]
for size in range(2,8):
    for members in combinations(range(8),size):
        outside=tuple(i for i in range(8) if i not in members)
        coefficients=[(sum(E.L[i][k][0] for i in members),sum(E.L[i][k][1] for i in members)) for k in range(13)]
        RECORDS.append((members,outside,coefficients))


def centroid(box,qbox,members,outside,coefficients):
    linear=E.ci(0)
    for value,(a,b) in zip(box,coefficients):
        linear=E.ca(linear,(E.scale(value,a),E.scale(value,b)))
    direct=E.ci(0);complement=E.ci(0)
    for i in members:direct=E.ca(direct,qbox[i])
    for i in outside:complement=E.cs(complement,qbox[i])
    return tuple(E.meet(E.meet(a,b),c) for a,b,c in zip(linear,direct,complement))


def propagate(box,qbox,settings,passes=2):
    if settings['scope']!='global_maximum':raise ValueError('subset inequalities require the global-maximizer reduction')
    q=[tuple(z) for z in qbox]
    norms=[(max(0,E.cnorm(z)[0]),min(settings['radius2'],E.cnorm(z)[1])) for z in q]
    original=[tuple(p) for p in box]
    changed=0;tested=0
    try:
        for step in range(passes):
            before=list(norms)
            for i in range(8):
                lower=max(norms[i][0],7*E.S-sum(norms[j][1] for j in range(8) if j!=i))
                upper=min(norms[i][1],7*E.S-sum(norms[j][0] for j in range(8) if j!=i))
                if lower>upper:return None,'subset_mass_budget',{'point':i}
                norms[i]=lower,upper
            for members,outside,coefficients in RECORDS:
                m=len(members)
                center=centroid(original,q,members,outside,coefficients)
                center_lower=E.cnorm(center)[0]
                target=m*m*(m-1)*E.S//8
                upper_mass=min(sum(norms[i][1] for i in members),7*E.S-sum(norms[i][0] for i in outside))
                upper=m*upper_mass-center_lower
                tested+=1
                if upper<=target:
                    return None,'subset_centroid_budget',{'members':list(members),'upper':upper,'strict_threshold_lower':target}
                required=(target+center_lower)//m
                total_upper=sum(norms[i][1] for i in members)
                for i in members:
                    lower=max(norms[i][0],required-total_upper+norms[i][1])
                    if lower>norms[i][1]:return None,'subset_mass_budget',{'point':i,'members':list(members)}
                    norms[i]=lower,norms[i][1]
                total_lower=sum(norms[i][0] for i in outside)
                for i in outside:
                    upper=min(norms[i][1],7*E.S-required-total_lower+norms[i][0])
                    if upper<norms[i][0]:return None,'subset_mass_budget',{'point':i,'members':list(members)}
                    norms[i]=norms[i][0],upper
            if norms==before:break
        for i in range(8):
            lo,hi=norms[i]
            for component in [0,1]:
                other=E.square(q[i][1-component])
                if hi<other[0]:return None,'subset_radius_budget',{'point':i}
                bound=E.sqrt_hi(hi-other[0])
                value=E.meet(q[i][component],(-bound,bound))
                if lo>other[1]:
                    minimum=E.sqrt_lo(lo-other[1])
                    if value[0]>=0:value=E.meet(value,(minimum,value[1]))
                    elif value[1]<=0:value=E.meet(value,(value[0],-minimum))
                if value!=q[i][component]:changed+=1
                q[i]=(value,q[i][1]) if component==0 else (q[i][0],value)
        E.close_center(q)
        result=[E.meet(a,b) for a,b in zip(original,E.variables(q))]
    except ValueError:return None,'subset_centroid_geometry',{}
    return result,None,{'point_boxes':q,'norms':norms,'tested':tested,'changed_components':changed}
