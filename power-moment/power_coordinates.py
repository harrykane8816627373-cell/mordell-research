"""Finite generating power sums and certified whole-box objective bounds."""
from fractions import Fraction as F
from complex_interval import C,I,as_i,inverse_det,root_bounds


def checked(p):
    n=len(p)-1
    if not 2<=n<=12:raise ValueError('implemented particle count 2..12')
    return n


def coefficients(p):
    """P(X)=sum_(j=0)^n c_j X^(n-j), c_0=1."""
    n=checked(p);c=[p[0]*0+1]
    for k in range(1,n+1):c.append(-sum((c[k-j]*p[j] for j in range(1,k+1)),0)/k)
    return c


def extend(p,last):
    n=checked(p)
    if type(last)is not int or last<n or last>2*n:raise ValueError('implemented extension n..2n')
    c=coefficients(p);q=list(p)
    for k in range(n+1,last+1):q.append(-sum((c[j]*q[k-j] for j in range(1,n+1)),0))
    return q


def hankel(p):
    n=checked(p);q=extend(p,max(n,2*n-2))
    return [[q[a+b] for b in range(n)] for a in range(n)]


def objective_box(center,radius):
    """Every actual configuration whose first n powers lie in this box.

    p_0=n and p_1=0 are fixed. Each remaining real/imag part has the same
    radius; the box is not assumed to consist of physical stationary data.
    """
    center=[C(x.re,x.im) if isinstance(x,C) else C(x) for x in center]
    n=checked(center);r=F(radius)
    if center[0]!=C(n) or center[1]!=C(0) or r<=0:raise ValueError('centered first n powers and positive box radius')
    P=[as_i(center[0]),as_i(center[1])]+[I((z.re-r,z.re+r),(z.im-r,z.im+r)) for z in center[2:]]
    H0=hankel(center);inv,det=inverse_det(H0);H=hankel(P)
    E=[[H[i][j]-as_i(H0[i][j]) for j in range(n)] for i in range(n)]
    rows=[]
    for i in range(n):
        total=F(0)
        for j in range(n):
            value=sum((as_i(inv[i][k])*E[k][j] for k in range(n)),as_i(0))
            total+=value.abs1_upper()
        rows.append(total)
    bound=root_bounds(det.norm())[1]
    for rsum in rows:bound*=1+rsum
    return {'objective_modulus_upper':bound,'center_determinant':det,
            'preconditioned_row_error_sums':rows,'radius':r,
            'scope':'Entire independent rectangular box in p2,...,pn, including all actual configurations it contains. This is an objective bound, not a stationarity or realizability certificate.'}
