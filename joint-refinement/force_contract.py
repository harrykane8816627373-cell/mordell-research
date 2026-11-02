"""Apply inherited physics, force affine exclusion and quadratic Newton contraction."""
import ellipse_dual as J
import force_enclosure as F

E=F.E


def subset_exclusion(qbox,inverse=None):
    if inverse is None:
        try:inverse,_=E.invert_edges(qbox)
        except ZeroDivisionError:return None
    weights=[E.cm(z,z) for z in inverse]
    for mask in range(1,128):
        count=mask.bit_count()
        values=[E.ci(0) for _ in range(8)]
        for (i,j),w in zip(E.PAIRS,weights):
            sign=bool(mask&(1<<i))-bool(mask&(1<<j))
            if sign:
                w=E.cscale(w,sign)
                values[i]=E.ca(values[i],w)
                values[j]=E.cs(values[j],w)
        lower=sum(E.cnorm(z)[0] for z in values)
        upper=2*count*(8-count)*E.S
        if lower>upper:return {'mask':mask,'lower':lower,'upper':upper}
    return None


def contract(original,settings,passes=3,endpoint_separation=None):
    box=[tuple(p) for p in original]
    for step in range(passes):
        box,reason,aux=E.contract(box,settings,passes=5)
        if box is None:return None,reason,aux
        if all(1000*E.cnorm(q)[0]>=729*E.S for q in aux[0]):
            return None,'Hessian_inner_729',{}
        data=F.enclose(box,aux[0],endpoint_separation)
        evidence=F.exclusion(data,'affine')
        if evidence is not None:return None,'force_affine_exact',evidence
        evidence=F.exclusion(data,'taylor')
        if evidence is not None:return None,'force_taylor_exact',evidence
        image=F.newton_image(data)
        if image is None:break
        try:new=[E.meet(a,b) for a,b in zip(box,image)]
        except ValueError:return None,'force_quadratic_Newton',{'empty_image':True}
        if new==box:break
        box=new
    box,reason,aux=J.contract(box,settings)
    if box is None:return box,reason,aux
    evidence=subset_exclusion(aux[0],aux[1])
    if evidence is not None:return None,'Hessian_subset_norm',evidence
    return box,None,aux
