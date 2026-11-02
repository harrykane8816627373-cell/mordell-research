"""Exact half-grid force affine enclosures and quadratic-remainder Newton maps."""
from fractions import Fraction as F
from math import isqrt
import refinement as R

E=R.E


def ceil_sqrt(n):
    a=isqrt(n)
    return a+(a*a!=n)


def gaussian(z,a,b):
    return E.sub(E.scale(z[0],a),E.scale(z[1],b)),E.add(E.scale(z[0],b),E.scale(z[1],a))


def half_product(a,b):
    values=[x*y for x in a for y in b]
    return min(values)//(2*E.S),-(-max(values)//(2*E.S))


def enclose(box,qbox=None,endpoint_separation=None):
    if len(box)!=13 or any(type(v) is not int for p in box for v in p) or any(a>b for a,b in box):
        raise ValueError('thirteen ordered exact intervals')
    if qbox is None:qbox=E.positions(box)
    center2=[a+b for a,b in box]
    cq=E.positions([E.point(v) for v in center2])
    mid=[(q[0][0],q[1][0]) for q in cq]
    radius=[]
    for q,(x,y) in zip(qbox,mid):
        rx=max(abs(2*q[0][0]-x),abs(2*q[0][1]-x))
        ry=max(abs(2*q[1][0]-y),abs(2*q[1][1]-y))
        radius.append(ceil_sqrt(rx*rx+ry*ry))
    value=[E.ci(0) for _ in range(8)]
    jacobian=[[E.ci(0) for _ in range(13)] for _ in range(8)]
    remainder=[0]*8
    first_error=[2*r for r in radius]
    for i,j in E.PAIRS:
        dx,dy=mid[i][0]-mid[j][0],mid[i][1]-mid[j][1]
        squared=dx*dx+dy*dy
        distance=isqrt(squared)
        perturbation=radius[i]+radius[j]
        if endpoint_separation is None:
            if distance<=perturbation:return None
            endpoint=distance-perturbation
        else:
            if type(endpoint_separation) is not int or endpoint_separation<0:
                raise ValueError('exact nonnegative endpoint squared-distance lower bound')
            measured=E.cnorm(E.cs(qbox[i],qbox[j]))[0]
            endpoint=max(distance-perturbation,2*E.sqrt_lo(max(measured,endpoint_separation)))
            if distance<=0 or endpoint<=0:return None
        linear_error=-(-2*E.S**2*perturbation//(distance*endpoint))
        quadratic_error=-(-2*E.S**2*perturbation**2//(distance**2*endpoint))
        first_error[i]+=linear_error
        first_error[j]+=linear_error
        remainder[i]+=quadratic_error
        remainder[j]+=quadratic_error
        inverse=E.rat(2*E.S*dx,squared),E.rat(-2*E.S*dy,squared)
        inverse2=E.rat(4*E.S**2*(dx*dx-dy*dy),squared*squared),E.rat(-8*E.S**2*dx*dy,squared*squared)
        value[i]=E.ca(value[i],inverse)
        value[j]=E.cs(value[j],inverse)
        for k in range(13):
            a,b=E.L[i][k][0]-E.L[j][k][0],E.L[i][k][1]-E.L[j][k][1]
            if a or b:
                v=gaussian(inverse2,a,b)
                jacobian[i][k]=E.cs(jacobian[i][k],v)
                jacobian[j][k]=E.ca(jacobian[j][k],v)
    for i in range(8):
        value[i]=E.cs(value[i],(E.point(2*mid[i][0]),E.point(-2*mid[i][1])))
        for k,(a,b) in enumerate(E.L[i]):
            jacobian[i][k]=E.ca(jacobian[i][k],(E.rat(-4*a),E.rat(4*b)))
    delta2=[(2*a-c,2*b-c) for (a,b),c in zip(box,center2)]
    affine=[]
    for i in range(8):
        z=value[i]
        for k in range(13):
            z=E.ca(z,(half_product(jacobian[i][k][0],delta2[k]),half_product(jacobian[i][k][1],delta2[k])))
        affine.append(z)
    return {'box':[tuple(p) for p in box],'center2':center2,'value':value,'jacobian':jacobian,
        'remainder':remainder,'first_error':first_error,'affine':affine,'delta2':delta2}


def exclusion(data,method='affine'):
    if data is None:return None
    if method=='taylor':
        for i,(value,error) in enumerate(zip(data['value'],data['first_error'])):
            if E.cnorm(value)[0]*E.S>error*error:
                return {'point':i,'method':'taylor','error':error}
        return None
    if method!='affine':raise ValueError('known force enclosure')
    for i,(value,error) in enumerate(zip(data['affine'],data['remainder'])):
        for component in [0,1]:
            lo,hi=value[component]
            if lo>error or hi<-error:
                return {'point':i,'component':component,'interval':list(value[component]),'error':error,
                    'signed_margin':lo-error if lo>error else hi+error,'method':'affine'}
    return None


def reduced(value,jacobian=None):
    rows=[(0,0)]+[(i,c) for i in range(1,7) for c in [0,1]]
    out=[]
    matrix=[]
    for i,c in rows:
        sign=-1 if c==0 else 1
        out.append(E.scale(E.sub(value[i][c],value[7][c]),sign))
        if jacobian is not None:
            matrix.append([E.scale(E.sub(jacobian[i][k][c],jacobian[7][k][c]),sign) for k in range(13)])
    return out,matrix


def preconditioner(matrix):
    # Invert a rational midpoint approximation. Correctness never assumes it is an exact inverse of the true Jacobian.
    a=[[F((x+y)//2) for x,y in row]+[F(int(i==j)) for j in range(13)] for i,row in enumerate(matrix)]
    for k in range(13):
        pivot=max(range(k,13),key=lambda i:abs(a[i][k]))
        if not a[pivot][k]:return None
        a[k],a[pivot]=a[pivot],a[k]
        p=a[k][k]
        a[k]=[v/p for v in a[k]]
        for i in range(13):
            if i!=k:
                factor=a[i][k]
                if factor:a[i]=[v-factor*w for v,w in zip(a[i],a[k])]
    return [[round(v*E.S**2) for v in row[13:]] for row in a]


def newton_image(data,target=None):
    if data is None:return None
    values,matrix=reduced(data['value'],data['jacobian'])
    if target is not None:
        tv,_=reduced(target)
        values=[E.sub(a,b) for a,b in zip(values,tv)]
    p=preconditioner(matrix)
    if p is None:return None
    rows=[(0,0)]+[(i,c) for i in range(1,7) for c in [0,1]]
    errors=[data['remainder'][i]+data['remainder'][7] for i,c in rows]
    image=[]
    for i in range(13):
        result=E.ZERO
        for j in range(13):result=E.sub(result,E.mul(E.point(p[i][j]),values[j]))
        for k in range(13):
            coefficient=E.rat(int(i==k))
            for j in range(13):coefficient=E.sub(coefficient,E.mul(E.point(p[i][j]),matrix[j][k]))
            result=E.add(result,half_product(coefficient,data['delta2'][k]))
        error=-(-sum(abs(p[i][j])*errors[j] for j in range(13))//E.S)
        c=data['center2'][i]
        image.append((c//2+result[0]-error,-(-c//2)+result[1]+error))
    return image
