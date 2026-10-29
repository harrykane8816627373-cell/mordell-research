"""Exact multipoint compressions of the original saturated stability operator.

A tuple (a,b) means a+i*sqrt(d)*b. All arithmetic is Fraction.
Passing a block test is necessary, never sufficient, for local maximality.
See BLOCK_STABILITY.md for the universal proof and scope of the examples.
"""
from fractions import Fraction as F
from itertools import combinations

ZERO = (F(0), F(0))
ONE = (F(1), F(0))


def add(z,w): return (z[0]+w[0], z[1]+w[1])
def neg(z): return (-z[0],-z[1])
def sub(z,w): return add(z,neg(w))
def scale(z,a): return (a*z[0],a*z[1])
def conj(z): return (z[0],-z[1])
def mul(z,w,d=1): return (z[0]*w[0]-d*z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def norm(z,d=1): return z[0]*z[0]+d*z[1]*z[1]
def divide(z,w,d=1): return scale(mul(z,conj(w),d),1/norm(w,d))


def geometry(points,d=1):
    if type(d) is not int or d<1:
        raise ValueError('positive integer d required')
    raw=tuple((F(a),F(b)) for a,b in points)
    n=len(raw)
    if n<2 or len(set(raw))!=n:
        raise ValueError('at least two distinct points required')
    mean=scale((sum(p[0] for p in raw),sum(p[1] for p in raw)),F(1,n))
    v=tuple(sub(p,mean) for p in raw)
    Q=sum(norm(p,d) for p in v)
    return {'n':n,'d':d,'points':v,'Q':Q,'lambda':F(n-1,2)}


def projection(data,i,j):
    n,Q,d=data['n'],data['Q'],data['d']
    v=data['points']
    return sub((F(int(i==j))-F(1,n),F(0)),scale(mul(v[i],conj(v[j]),d),1/Q))


def residual_cross(data,i,j):
    if i==j:
        raise ValueError('cross entries require distinct indices')
    n,Q,d=data['n'],data['Q'],data['d']
    v,l=data['points'],data['lambda']
    delta=sub(v[i],v[j])
    inverse=divide(ONE,mul(delta,delta,d),d)
    return sub(scale(inverse,-Q/n),scale(mul(conj(v[i]),conj(v[j]),d),l/Q))


def indices(data,rows,columns):
    rows,columns=tuple(rows),tuple(columns)
    all_indices=rows+columns
    if not rows or not columns or len(set(all_indices))!=len(all_indices):
        raise ValueError('nonempty disjoint row and column index sets required')
    if any(type(i) is not int or not 0<=i<data['n'] for i in all_indices):
        raise ValueError('point index outside range')
    return rows,columns


def cross_block_matrix(data,rows,columns):
    rows,columns=indices(data,rows,columns)
    l=data['lambda']
    matrix=[]
    for i in rows:
        matrix.append([scale(conj(projection(data,i,j)),l) for j in rows]
                      +[residual_cross(data,i,j) for j in columns])
    for i in columns:
        matrix.append([conj(residual_cross(data,j,i)) for j in rows]
                      +[scale(projection(data,i,j),l) for j in columns])
    return matrix


def row_schur_matrix(data,row,columns):
    rows,columns=indices(data,(row,),columns)
    l,d=data['lambda'],data['d']
    p=projection(data,row,row)[0]
    c=[residual_cross(data,row,j) for j in columns]
    return [[sub(scale(projection(data,j,k),l*l*p),mul(conj(c[a]),c[b],d))
             for b,k in enumerate(columns)] for a,j in enumerate(columns)]


def determinant(matrix,d=1):
    a=[list(row) for row in matrix]
    n=len(a)
    if any(len(row)!=n for row in a):
        raise ValueError('square matrix required')
    result=ONE
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]!=ZERO),None)
        if pivot is None:
            return ZERO
        if pivot!=j:
            a[j],a[pivot]=a[pivot],a[j]
            result=neg(result)
        result=mul(result,a[j][j],d)
        for i in range(j+1,n):
            ratio=divide(a[i][j],a[j][j],d)
            for k in range(j+1,n):
                a[i][k]=sub(a[i][k],mul(ratio,a[j][k],d))
            a[i][j]=ZERO
    return result


def quadratic(matrix,vector,d=1):
    n=len(matrix)
    if len(vector)!=n:
        raise ValueError('matching matrix and vector required')
    result=ZERO
    for i in range(n):
        for j in range(n):
            result=add(result,mul(conj(vector[i]),mul(matrix[i][j],vector[j],d),d))
    if result[1]:
        raise ArithmeticError('Hermitian quadratic form must be real')
    return result[0]


def negative_vector(matrix,d=1):
    """Return a rational negative quadratic witness if positive LDL pivots fail.

    Zero-pivot degeneracies return None; complete principal-minor acceptance
    below does not rely on this optional witness extraction.
    """
    n=len(matrix)
    L=[[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
    D=[]
    for j in range(n):
        pivot=matrix[j][j][0]-sum(norm(L[j][k],d)*D[k] for k in range(j))
        if pivot<0:
            x=[ZERO]*n
            x[j]=ONE
            for i in range(j-1,-1,-1):
                x[i]=neg(sum_complex((mul(conj(L[k][i]),x[k],d) for k in range(i+1,j+1))))
            q=quadratic(matrix,x,d)
            if q!=pivot:
                raise ArithmeticError('LDL negative witness mismatch')
            return {'vector':x,'quadratic_value':q,'first_negative_pivot':j}
        if pivot==0:
            return None
        D.append(pivot)
        for i in range(j+1,n):
            t=matrix[i][j]
            for k in range(j):
                t=sub(t,scale(mul(L[i][k],conj(L[j][k]),d),D[k]))
            L[i][j]=scale(t,1/pivot)
    return None


def sum_complex(values):
    result=ZERO
    for value in values:
        result=add(result,value)
    return result


def psd_certificate(matrix,d=1,max_minors=4095):
    n=len(matrix)
    if not n or any(len(row)!=n for row in matrix):
        raise ValueError('nonempty square matrix required')
    if any(matrix[i][j]!=conj(matrix[j][i]) for i in range(n) for j in range(n)):
        raise ValueError('Hermitian matrix required')
    planned=(1<<n)-1
    if type(max_minors) is not int or max_minors<1 or planned>max_minors:
        raise ValueError('complete PSD check exceeds principal-minor budget')
    checked=0
    for size in range(1,n+1):
        for ids in combinations(range(n),size):
            q=determinant([[matrix[i][j] for j in ids] for i in ids],d)
            if q[1]:
                raise ArithmeticError('Hermitian principal determinant must be real')
            checked+=1
            if q[0]<0:
                return {'passes':False,'negative_principal_indices':ids,
                        'negative_principal_determinant':q[0],
                        'principal_minors_checked':checked,
                        'negative_quadratic_witness':negative_vector(matrix,d)}
    return {'passes':True,'principal_minors_checked':checked,
            'scope':'Necessary compressed stability only, not local or global maximality.'}


def cross_block_certificate(points,rows,columns,d=1,max_minors=4095):
    data=geometry(points,d)
    matrix=cross_block_matrix(data,rows,columns)
    return {'n':data['n'],'rows':tuple(rows),'columns':tuple(columns),
            'matrix':matrix,**psd_certificate(matrix,d,max_minors)}


def triangle_certificate(points,row,j,k,d=1):
    data=geometry(points,d)
    matrix=row_schur_matrix(data,row,(j,k))
    det=determinant(matrix,d)
    return {'n':data['n'],'row':row,'columns':(j,k),'schur_matrix':matrix,
            'schur_determinant':det[0],**psd_certificate(matrix,d)}


# Real and complex rectangular interval operations for entire classes.
def iv(q): return (F(q),F(q))
def iadd(a,b): return (a[0]+b[0],a[1]+b[1])
def ineg(a): return (-a[1],-a[0])
def isub(a,b): return iadd(a,ineg(b))
def imul(a,b):
    values=[x*y for x in a for y in b]
    return min(values),max(values)
def iscale(a,c): return imul(a,iv(c))
def isquare(a):
    return (F(0) if a[0]<=0<=a[1] else min(a[0]*a[0],a[1]*a[1]),max(a[0]*a[0],a[1]*a[1]))
def iinverse(a):
    if a[0]<=0<=a[1]:
        raise ValueError('interval reciprocal crosses zero')
    return F(1)/a[1],F(1)/a[0]
def izadd(z,w): return (iadd(z[0],w[0]),iadd(z[1],w[1]))
def izneg(z): return (ineg(z[0]),ineg(z[1]))
def izsub(z,w): return izadd(z,izneg(w))
def izconj(z): return (z[0],ineg(z[1]))
def izscale(z,a): return (imul(z[0],a),imul(z[1],a))
def izmul(z,w): return (isub(imul(z[0],w[0]),imul(z[1],w[1])),iadd(imul(z[0],w[1]),imul(z[1],w[0])))
def iznorm(z): return iadd(isquare(z[0]),isquare(z[1]))


def box_geometry(point_boxes):
    raw=tuple(((F(a),F(b)),(F(c),F(e))) for (a,b),(c,e) in point_boxes)
    n=len(raw)
    if n<2 or any(a>b or c>e for (a,b),(c,e) in raw):
        raise ValueError('at least two ordered complex coordinate boxes required')
    mean=(iv(0),iv(0))
    for p in raw:
        mean=izadd(mean,p)
    mean=izscale(mean,iv(F(1,n)))
    v=tuple(izsub(p,mean) for p in raw)
    Q=iv(0)
    for p in v:
        Q=iadd(Q,iznorm(p))
    if Q[0]<=0:
        raise ValueError('box moment must be strictly positive')
    return {'n':n,'raw':raw,'points':v,'Q':Q,'lambda':F(n-1,2)}


def box_projection(data,i,j):
    n,Q,v=data['n'],data['Q'],data['points']
    return izsub((iv(F(int(i==j))-F(1,n)),iv(0)),
                 izscale(izmul(v[i],izconj(v[j])),iinverse(Q)))


def box_cross(data,i,j):
    n,Q,l=data['n'],data['Q'],data['lambda']
    # Translation cancels in a difference; avoid artificial mean dependence.
    delta=izsub(data['raw'][i],data['raw'][j])
    inverse=izscale(izconj(delta),iinverse(iznorm(delta)))
    first=izscale(izmul(inverse,inverse),iscale(Q,-F(1,n)))
    second=izscale(izmul(izconj(data['points'][i]),izconj(data['points'][j])),iscale(iinverse(Q),l))
    return izsub(first,second)


def triangle_box_certificate(point_boxes,row,j,k):
    data=box_geometry(point_boxes)
    indices({'n':data['n']},(row,),(j,k))
    l=data['lambda']
    pi=box_projection(data,row,row)[0]
    pj=box_projection(data,j,j)[0]
    pk=box_projection(data,k,k)[0]
    cj,ck=box_cross(data,row,j),box_cross(data,row,k)
    hj=isub(iscale(imul(pi,pj),l*l),iznorm(cj))
    hk=isub(iscale(imul(pi,pk),l*l),iznorm(ck))
    cross=izsub(izscale(box_projection(data,j,k),iscale(pi,l*l)),izmul(izconj(cj),ck))
    det=isub(imul(hj,hk),iznorm(cross))
    return {'n':data['n'],'row':row,'columns':(j,k),'diagonal_intervals':(hj,hk),
            'determinant_interval':det,'excludes_entire_box':det[1]<0}


def all_pair_box_certificate(point_boxes):
    """Lower bounds for ALL old scalar pair slacks on a shape-coordinate box."""
    data=box_geometry(point_boxes)
    n,Q=data['n'],data['Q']
    factor=iscale(iinverse(Q),n)
    radii=[imul(factor,iznorm(v)) for v in data['points']]
    s0=F(2*n,(n-1)**2)
    rows=[]
    for i,j in combinations(range(n),2):
        s=imul(factor,iznorm(izsub(data['raw'][i],data['raw'][j])))
        a,b=radii[i],radii[j]
        bracket=iadd(iv(1),iscale(imul(isub(iv(n-1),iadd(a,b)),s),F(n-1,2*n)))
        slack=isub(iscale(imul(isub(s,iv(s0)),bracket),n-1),isquare(isub(a,b)))
        rows.append({'pair':(i,j),'distance_lower':s[0],'slack_lower':slack[0]})
    return {'all_old_local_pairs_pass':all(r['distance_lower']>s0 and r['slack_lower']>0 for r in rows),
            'minimum_distance_lower':min(r['distance_lower'] for r in rows),
            'maximum_radius_upper':max(a[1] for a in radii),
            'all_old_global_geometry_pass':n>=3 and min(r['distance_lower'] for r in rows)>F(4,n-1)
                and max(a[1] for a in radii)<F(2*(n-1),n),
            'minimum_slack_lower':min(r['slack_lower'] for r in rows),'rows':rows}


def cross_block_box_quadratic(point_boxes,rows,columns,vector):
    """Upper-bound a FIXED rational quadratic witness on an entire box."""
    data=box_geometry(point_boxes)
    rows,columns=indices({'n':data['n']},rows,columns)
    l=data['lambda']
    matrix=[]
    for i in rows:
        matrix.append([izscale(izconj(box_projection(data,i,j)),iv(l)) for j in rows]
                      +[box_cross(data,i,j) for j in columns])
    for i in columns:
        matrix.append([izconj(box_cross(data,j,i)) for j in rows]
                      +[izscale(box_projection(data,i,j),iv(l)) for j in columns])
    if len(vector)!=len(matrix):
        raise ValueError('matching fixed complex witness vector required')
    v=[(iv(F(a)),iv(F(b))) for a,b in vector]
    result=(iv(0),iv(0))
    for i in range(len(v)):
        for j in range(len(v)):
            result=izadd(result,izmul(izconj(v[i]),izmul(matrix[i][j],v[j])))
    return {'quadratic_interval':result[0],'excludes_entire_box':result[0][1]<0,
            'rows':rows,'columns':columns,'fixed_vector':vector}
