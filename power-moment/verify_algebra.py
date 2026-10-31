"""Independent direct-root and division-free determinant checks."""
from fractions import Fraction as F
from itertools import product
from acceptance import require,write_report
from complex_interval import C,I,as_i,root_bounds,inverse_det
from power_coordinates import coefficients,extend,hankel,objective_box


def det_dp(A):
    n=len(A);states={0:C(1)}
    for mask in range(1,1<<n):
        row=mask.bit_count()-1;value=C(0)
        for col in range(n):
            if mask>>col&1:
                before=mask^(1<<col)
                sign=(-1)**((mask>>(col+1)).bit_count())
                value+=sign*states[before]*A[row][col]
        states[mask]=value
    return states[(1<<n)-1]


def run():
    cases=0
    for n in range(2,11):
        roots=[C(j,j*j%7) for j in range(n)]
        average=sum(roots,C(0))/n;roots=[z-average for z in roots]
        p=[sum((z**k for z in roots),C(0)) for k in range(n+1)]
        c=coefficients(p);direct=[C(1)]
        for z in roots:
            new=[C(0)]*(len(direct)+1)
            for j,a in enumerate(direct):new[j]+=a;new[j+1]-=z*a
            direct=new
        require(c==direct,'Newton coefficients versus direct root multiplication')
        tail=extend(p,2*n-2 if n>2 else n)
        require(all(tail[k]==sum((z**k for z in roots),C(0)) for k in range(len(tail))),'all higher powers constrained by the same recurrence')
        determinant=det_dp(hankel(p));vand=C(1)
        for i in range(n):
            for j in range(i):vand*=((roots[i]-roots[j])**2)
        require(determinant==vand,'holomorphic Hankel determinant versus independent distance product')
        if n<=6:
            inv,d=inverse_det(hankel(p));require(d==determinant,'Gaussian elimination versus division-free determinant')
        cases+=1
    corner_checks=0
    for n in (2,3,4):
        center=[C(n),C(0)]+[C(F(k,3),F(1,k+1)) for k in range(2,n+1)]
        radius=F(1,100000);bound=objective_box(center,radius)['objective_modulus_upper']
        for signs in product((-1,1),repeat=2*(n-1)):
            point=center[:2]+[center[k]+C(signs[2*(k-2)]*radius,signs[2*(k-2)+1]*radius) for k in range(2,n+1)]
            d=det_dp(hankel(point));require(d.norm()<=bound*bound,'entire-box bound checked against independent corner determinants')
            corner_checks+=1
    interval_checks=0
    A=I((F(-2,3),F(4,5)),(F(-3,4),F(5,6)))
    B=I((F(-1,5),F(7,8)),(F(-2,7),F(3,8)))
    for ar,ai,br,bi in product(A.re,A.im,B.re,B.im):
        require((A*B).contains(C(ar,ai)*C(br,bi)),'independent rectangular product corner')
        interval_checks+=1
    # Exact stationary examples verify Ward separately from the coordinate engine.
    ward=0
    for roots,kappa in [([C(-1),C(1)],F(1,2)),
                        ([C(1),C(0,1),C(-1),C(0,-1)],F(3,2)),
                        ([C(0),C(1),C(0,1),C(-1),C(0,-1)],F(5,2))]:
        n=len(roots);powers=[sum((z**k for z in roots),C(0)) for k in range(2*n+1)]
        for z in roots:require(sum((1/(z-y) for y in roots if y!=z),C(0))==kappa*z.conjugate(),'actual rational stationary force')
        for k in range(1,n+2):
            lhs=2*kappa*sum((z.conjugate()*z**k for z in roots),C(0))
            rhs=sum((powers[j]*powers[k-1-j] for j in range(k)),C(0))-k*powers[k-1]
            require(lhs==rhs,'general n,k mixed Ward identity');ward+=1
    rejected=0
    for action in [lambda:objective_box([C(2),C(1),C(3)],F(1,100)),
                   lambda:objective_box([C(2),C(0),C(3)],0),
                   lambda:I((2,1)),lambda:extend([C(2),C(0),C(1)],1)]:
        try:action()
        except ValueError:rejected+=1
    require(rejected==4,'illegal input rejection')
    return {'passed':True,'independent_root_cases':cases,'independent_box_corner_determinants':corner_checks,
            'complex_interval_corner_checks':interval_checks,'stationary_Ward_checks':ward,
            'invalid_inputs_rejected':rejected,
            'scope':'General identities are analytic; finite checks test implementations. No enumeration of all n8 stationary solutions.'}


if __name__=='__main__':write_report('algebra_verification.json',run())
