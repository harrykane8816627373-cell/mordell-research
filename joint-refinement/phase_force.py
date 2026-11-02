"""Complex interval remainders in the exact finite reciprocal identity."""
import edge_force as S
import force_enclosure as F

E=F.E


def enclose(box,qbox,separation,center2=None):
    data=S.enclose(box,qbox,separation,center2)
    if data is None:return None
    if any(E.cnorm(E.cs(qbox[i],qbox[j]))[1]<separation for i,j in E.PAIRS):return None
    inverse,_=E.invert_edges(qbox,separation)
    cq=E.positions([E.point(v) for v in data['center2']])
    remainder=[E.ci(0) for _ in range(8)]
    for record,zinv,(i,j) in zip(data['edge_bounds'],inverse,E.PAIRS):
        dx=cq[i][0][0]-cq[j][0][0];dy=cq[i][1][0]-cq[j][1][0]
        squared=dx*dx+dy*dy
        inv2=(E.rat(4*E.S**2*(dx*dx-dy*dy),squared*squared),
            E.rat(-8*E.S**2*dx*dy,squared*squared))
        x,y=E.divide(record['real2'],2),E.divide(record['imag2'],2)
        h2=E.sub(E.square(x),E.square(y)),E.scale(E.mul(x,y),2)
        term=E.cm(E.cm(h2,inv2),zinv)
        remainder[i]=E.ca(remainder[i],term)
        remainder[j]=E.cs(remainder[j],term)
    try:
        data['complex_remainder']=[tuple(E.meet(v,(-r,r)) for v in z)
            for z,r in zip(remainder,data['remainder'])]
    except ValueError:return None
    return data


def exclusion(data):
    if data is None:return None
    for i,(affine,remainder) in enumerate(zip(data['affine'],data['complex_remainder'])):
        for c,interval in enumerate(E.ca(affine,remainder)):
            if interval[0]>0 or interval[1]<0:
                return {'point':i,'component':c,'force_interval':list(interval)}
    return None


def newton_image(data,target=None):
    if data is None:return None
    values,matrix=F.reduced(data['value'],data['jacobian'])
    if target is not None:
        tv,_=F.reduced(target)
        values=[E.sub(a,b) for a,b in zip(values,tv)]
    p=F.preconditioner(matrix)
    if p is None:return None
    errors,_=F.reduced(data['complex_remainder'])
    image=[]
    for i in range(13):
        result=E.ZERO
        for j in range(13):
            result=E.sub(result,E.mul(E.point(p[i][j]),E.add(values[j],errors[j])))
        for k in range(13):
            coefficient=E.rat(int(i==k))
            for j in range(13):coefficient=E.sub(coefficient,E.mul(E.point(p[i][j]),matrix[j][k]))
            result=E.add(result,F.half_product(coefficient,data['delta2'][k]))
        c=data['center2'][i]
        image.append((c//2+result[0],-(-c//2)+result[1]))
    return image


def contract(original,settings,passes=2):
    box=[tuple(p) for p in original]
    for step in range(passes):
        box,reason,aux=E.contract(box,settings,passes=5)
        if box is None:return None,reason,aux
        if all(1000*E.cnorm(q)[0]>=729*E.S for q in aux[0]):
            return None,'Hessian_inner_729',{}
        before=list(box)
        for source,center in S.references(before,aux[0],settings):
            data=enclose(before,aux[0],settings['separation'],center)
            evidence=exclusion(data)
            if evidence is not None:return None,'phase_force_affine',dict(evidence,reference=source)
            evidence=F.exclusion(data,'taylor')
            if evidence is not None:return None,'phase_force_taylor',dict(evidence,reference=source)
            if source.startswith('reference_'):continue
            image=newton_image(data)
            if image is None:continue
            try:box=[E.meet(a,b) for a,b in zip(box,image)]
            except ValueError:return None,'phase_force_Newton',{'reference':source}
        if box==before:break
    box,reason,aux=S.K.J.contract(box,settings)
    if box is None:return box,reason,aux
    evidence=S.K.subset_exclusion(aux[0],aux[1])
    if evidence is not None:return None,'Hessian_subset_norm',evidence
    return box,None,aux
