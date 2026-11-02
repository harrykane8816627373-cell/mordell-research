"""Force enclosures retaining exact edge correlations and arbitrary rational references."""
from math import isqrt
import force_enclosure as F
import force_contract as K

E=F.E


def scaled_interval(interval,n):
    a,b=interval
    return (a*n,b*n) if n>=0 else (b*n,a*n)


def enclose(box,qbox,separation,center2=None):
    if center2 is None:center2=[a+b for a,b in box]
    if len(center2)!=13 or any(type(v) is not int for v in center2):
        raise ValueError('thirteen exact doubled-grid reference coordinates')
    reference_box=[(v//2,-(-v//2)) for v in center2]
    data=F.enclose(reference_box,qbox,endpoint_separation=separation)
    if data is None:return None
    if data['center2']!=list(center2):raise ArithmeticError('exact half-grid reference')
    cq=E.positions([E.point(v) for v in center2])
    mid=[(z[0][0],z[1][0]) for z in cq]
    radii=[]
    for z,(x,y) in zip(qbox,mid):
        rx=max(abs(2*z[0][0]-x),abs(2*z[0][1]-x))
        ry=max(abs(2*z[1][0]-y),abs(2*z[1][1]-y))
        radii.append(F.ceil_sqrt(rx*rx+ry*ry))
    delta2=[(2*a-c,2*b-c) for (a,b),c in zip(box,center2)]
    remainder=[0]*8
    first=[2*r for r in radii]
    bounds=[]
    for edge,(i,j) in enumerate(E.PAIRS):
        dx,dy=mid[i][0]-mid[j][0],mid[i][1]-mid[j][1]
        distance=isqrt(dx*dx+dy*dy)
        real,imag=(0,0),(0,0)
        for (a,b),d in zip(E.EDGE_L[edge],delta2):
            real=E.add(real,scaled_interval(d,a))
            imag=E.add(imag,scaled_interval(d,b))
        actual=E.cs(qbox[i],qbox[j])
        xr=(2*actual[0][0]-dx,2*actual[0][1]-dx)
        yr=(2*actual[1][0]-dy,2*actual[1][1]-dy)
        real=max(real[0],xr[0]),min(real[1],xr[1])
        imag=max(imag[0],yr[0]),min(imag[1],yr[1])
        if real[0]>real[1] or imag[0]>imag[1]:return None
        rectangle=F.ceil_sqrt(max(abs(v) for v in real)**2+max(abs(v) for v in imag)**2)
        perturbation=min(radii[i]+radii[j],rectangle)
        endpoint=max(distance-perturbation,2*E.sqrt_lo(max(separation,E.cnorm(actual)[0])))
        if distance<=0 or endpoint<=0:return None
        error1=-(-2*E.S**2*perturbation//(distance*endpoint))
        error2=-(-2*E.S**2*perturbation**2//(distance**2*endpoint))
        first[i]+=error1;first[j]+=error1
        remainder[i]+=error2;remainder[j]+=error2
        bounds.append({'pair':[i,j],'real2':real,'imag2':imag,'upper':perturbation,
            'previous_upper':radii[i]+radii[j]})
    affine=[]
    for i in range(8):
        value=data['value'][i]
        for k in range(13):
            value=E.ca(value,(F.half_product(data['jacobian'][i][k][0],delta2[k]),
                F.half_product(data['jacobian'][i][k][1],delta2[k])))
        affine.append(value)
    data.update({'box':[tuple(p) for p in box],'delta2':delta2,'affine':affine,
        'remainder':remainder,'first_error':first,'edge_bounds':bounds})
    return data


def references(box,qbox,settings):
    result=[('midpoint',[a+b for a,b in box])]
    try:
        midpoint,_=E.slab_center(box,settings['radius'],qbox[7])
        candidate=[2*v for v in midpoint]
        if candidate!=result[0][1]:result.append(('centroid_slab',candidate))
    except (ValueError,ArithmeticError):pass
    midpoint=[(a+b)//2 for a,b in box]
    nearest=min(enumerate(K.J.R.REFERENCES['centers']),
        key=lambda item:sum((a-b)**2 for a,b in zip(midpoint,item[1])))
    candidate=[2*v for v in nearest[1]]
    if all(candidate!=v for name,v in result):result.append(('reference_'+str(nearest[0]),candidate))
    return result


def contract(original,settings,passes=2):
    box=[tuple(p) for p in original]
    for step in range(passes):
        box,reason,aux=E.contract(box,settings,passes=5)
        if box is None:return None,reason,aux
        if all(1000*E.cnorm(q)[0]>=729*E.S for q in aux[0]):
            return None,'Hessian_inner_729',{}
        before=list(box)
        for source,center in references(before,aux[0],settings):
            data=enclose(before,aux[0],settings['separation'],center)
            evidence=F.exclusion(data,'affine')
            if evidence is not None:return None,'edge_force_affine',dict(evidence,reference=source)
            evidence=F.exclusion(data,'taylor')
            if evidence is not None:return None,'edge_force_taylor',dict(evidence,reference=source)
            if source.startswith('reference_'):continue
            image=F.newton_image(data)
            if image is None:continue
            try:box=[E.meet(a,b) for a,b in zip(box,image)]
            except ValueError:return None,'edge_force_Newton',{'reference':source}
        if box==before:break
    box,reason,aux=K.J.contract(box,settings)
    if box is None:return box,reason,aux
    evidence=K.subset_exclusion(aux[0],aux[1])
    if evidence is not None:return None,'Hessian_subset_norm',evidence
    return box,None,aux
