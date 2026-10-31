"""Exact central-heptagon spectrum and a quantitative power-sum neighborhood."""
from fractions import Fraction as F
from math import comb
from acceptance import require,write_report


def reduce(v):
    a=list(v)+[F(0)]*(7-len(v));last=a[6]
    return tuple(F(x)-last for x in a)
Z=reduce([0]);ONE=reduce([1])
def add(a,b):return reduce([x+y for x,y in zip(a,b)])
def scale(a,c):return tuple(c*x for x in a)
def mul(a,b):
    out=[F(0)]*7
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[(i+j)%7]+=x*y
    return reduce(out)
def omega(k):
    out=[F(0)]*7;out[k%7]=F(1);return reduce(out)
def inv_one_minus(k):return scale(sum_fields([scale(omega(k*j),F(j)) for j in range(1,7)]),F(-1,7))
def sum_fields(rows):
    result=Z
    for row in rows:result=add(result,row)
    return result
def conjugate(a):
    out=[F(0)]*7
    for k,x in enumerate(a):out[-k%7]+=x
    return reduce(out)
def mm(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),F(0)) for j in range(len(B[0]))] for i in range(len(A))]
def poly_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def run():
    # Cyclotomic field Q[zeta_7], with 1+...+zeta_7^6=0.
    for k in range(1,7):require(mul(add(ONE,scale(omega(k),-1)),inv_one_minus(k))==ONE,'exact inverse in cyclotomic field')
    points=[Z]+[scale(omega(j),2) for j in range(7)]
    inverse=[[Z]*8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            if i==j:continue
            if i==0:inverse[i][j]=scale(omega(-(j-1)),F(-1,2))
            elif j==0:inverse[i][j]=scale(omega(-(i-1)),F(1,2))
            else:inverse[i][j]=scale(mul(omega(-(i-1)),inv_one_minus(j-i)),F(1,2))
    for i in range(8):
        bar_point=Z if i==0 else scale(omega(-(i-1)),2)
        require(sum_fields(inverse[i])==bar_point,'actual eight-point force, exactly')
    B=[[scale(mul(inverse[i][j],inverse[i][j]),-1) if i!=j else Z for j in range(8)] for i in range(8)]
    for i in range(8):B[i][i]=scale(sum_fields(B[i]),-1)
    basis=[[ONE]+[Z]*7]+[[Z]+[omega(k*j) for j in range(7)] for k in range(7)]
    for i,u in enumerate(basis):
        for j,v in enumerate(basis):
            actual=sum_fields([mul(conjugate(x),y) for x,y in zip(u,v)])
            expected=scale(ONE,1 if i==0 else 7) if i==j else Z
            require(actual==expected,'exact Fourier-basis Hermitian metric')
    M=[[F(0)]*8 for _ in range(8)]
    for row,col,value in [(6,0,F(-1,4)),(6,1,F(1,4)),(7,2,F(1)),
            (0,3,F(-7,4)),(1,3,F(1,4)),(2,4,F(-1,4)),
            (3,5,F(-1,2)),(4,6,F(-1,2)),(5,7,F(-1,4))]:M[row][col]=value
    for k,v in enumerate(basis):
        actual=[sum_fields([mul(B[i][j],v[j]) for j in range(8)]) for i in range(8)]
        expected=[sum_fields([scale(basis[j][i],M[j][k]) for j in range(8)]) for i in range(8)]
        require(actual==expected,'independent cyclotomic action versus rational Fourier matrix')
    metric=[1]+[7]*7
    T=[[sum((M[k][i]*metric[k]*M[k][j] for k in range(8)),F(0))/metric[i] for j in range(8)] for i in range(8)]
    power=[[F(int(i==j)) for j in range(8)] for i in range(8)];traces=[F(8)]
    for j in range(1,9):power=mm(power,T);traces.append(sum((power[i][i] for i in range(8)),F(0)))
    coeff=[F(1)]
    for k in range(1,9):coeff.append(-sum((coeff[j]*traces[k-j] for j in range(k)),F(0))/k)
    expected=[F(1)]
    for eigen in [0,1,F(1,2),F(1,2),F(1,4),F(1,4),F(1,16),F(1,16)]:expected=poly_mul(expected,[1,-F(eigen)])
    require(coeff==expected,'exact complete singular-square characteristic polynomial')
    # Rational pi enclosure from Machin, with opposite alternating tails.
    def atan_bounds(x,N):
        a=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
        b=a+(-1)**N*x**(2*N+1)/F(2*N+1)
        return min(a,b),max(a,b)
    a=atan_bounds(F(1,5),12);b=atan_bounds(F(1,239),6)
    pl,pu=16*a[0]-4*b[1],16*a[1]-4*b[0]
    require(pl/F(7)-(pu/F(7))**3/6>F(3,7),'minimum root separation greater than 12/7')
    rho=F(1,50);d0=F(12,7);r=F(1,150);eps=F(1,10)
    perturb=56*rho/(d0-2*rho)**3
    require(F(25,49)>F(1,2) and perturb<F(2,7),'strict Hessian gap throughout the rotation slice ball')
    bounds=[F(0)]*9
    for k in range(2,9):bounds[k]=(eps+eps*sum(bounds[2:k-1]))/k
    inner_boundary=r*(128-r**7)
    outer_boundary=(2-r)*(448*r-sum((F(comb(7,k))*2**(7-k)*r**k for k in range(2,8)),F(0)))
    inner_error=sum((bounds[k]*r**(8-k) for k in range(2,9)),F(0))
    outer_error=sum((bounds[k]*(2+r)**(8-k) for k in range(2,9)),F(0))
    require(inner_error<inner_boundary and outer_error<outer_boundary and 8*r*r<rho*rho,'Rouche disks imply the Euclidean isolation ball')
    require(0<eps<=1,'Newton coefficient induction domain')
    return {'passed':True,'n':8,'force_constant':1,'reference_polynomial':'X^8-128*X',
        'reference_powers':'p1=...=p6=p8=0, p7=896',
        'singular_squares':[0,1,F(1,2),F(1,2),F(1,4),F(1,4),F(1,16),F(1,16)],
        'characteristic_coefficients_descending':coeff,
        'strict_slice_hessian_gap_lower':'2/7','Euclidean_isolation_radius':rho,
        'root_disk_radius':r,'power_coordinate_radius':eps,
        'Hessian_perturbation_norm_upper':perturb,
        'Newton_coefficient_error_upper':bounds,
        'Rouche_inner_boundary_lower':inner_boundary,'Rouche_inner_error_upper':inner_error,
        'Rouche_outer_boundary_lower':outer_boundary,'Rouche_outer_error_upper':outer_error,
        'scope':'Any actual force-one eight-point stationary configuration in this full power-sum neighborhood is a rotated/relabelled central regular heptagon. This excludes its punctured nonzero-p2/p3 stationary neighborhood; it does not cover every remaining power vector.'}


if __name__=='__main__':write_report('isolation_verification.json',run())
