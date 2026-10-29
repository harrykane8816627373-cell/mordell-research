"""Overlapping annular angular tests and certified crossing-disk averages.

All accepted bounds use Fraction arithmetic. Proof: OVERLAPPING_ENERGY.md.
The full Gram matrix is retained. Disk averages remain valid when a disk
crosses any number of test interfaces; no center-value substitution is made.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'joint-bounds'))
from subset_moment import log_interval,outward_decimal


def require_mode(k):
    if type(k)is not int or k<1:
        raise ValueError('integer angular mode k>=1 required')


def validate_hat(hat):
    a,b,c=map(F,hat)
    if not 0<=a<b<c:
        raise ValueError('squared radii 0<=a<b<c required')
    return a,b,c


def sqrt_upper(q,scale=10**6):
    q=F(q)
    if q<0 or type(scale)is not int or scale<1:
        raise ValueError('nonnegative input and positive root scale required')
    value=isqrt(q.numerator*scale*scale//q.denominator)
    if value*value*q.denominator<q.numerator*scale*scale:
        value+=1
    return F(value,scale)


def circle_coefficients(k,hat):
    """Phi_hat=-K_a+(1+t)K_b-tK_c; K_0 is identically zero."""
    require_mode(k)
    a,b,c=validate_hat(hat)
    t=(b**k-a**k)/(c**k-b**k)
    rows=[(b,1+t),(c,-t)]
    if a:
        rows.insert(0,(a,F(-1)))
    return rows


def gram_matrix(k,hats):
    """Exact Gram of arbitrary overlapping hats, for real OR imaginary part."""
    require_mode(k)
    hats=[validate_hat(hat) for hat in hats]
    if not hats:
        raise ValueError('nonempty hat family required')
    coefficients=[circle_coefficients(k,hat) for hat in hats]
    return [[k*sum((v*w*min(s,t)**k for s,v in left for t,w in right),F(0))
             for right in coefficients] for left in coefficients]


def solve_positive(matrix,rhs):
    """Exact positive-definite Gaussian elimination; singular tests rejected."""
    n=len(matrix)
    if n==0 or len(rhs)!=n or any(len(row)!=n for row in matrix):
        raise ValueError('nonempty square Gram and matching right side required')
    a=[list(map(F,row))+[F(value)] for row,value in zip(matrix,rhs)]
    for j in range(n):
        if a[j][j]<=0:
            raise ValueError('dependent tests or nonpositive Gram pivot')
        pivot=a[j][j]
        for i in range(j+1,n):
            ratio=a[i][j]/pivot
            for ell in range(j,n+1):
                a[i][ell]-=ratio*a[j][ell]
    result=[F(0)]*n
    for j in range(n-1,-1,-1):
        result[j]=(a[j][n]-sum((a[j][ell]*result[ell]
                              for ell in range(j+1,n)),F(0)))/a[j][j]
    return result


def phi_value(k,hat,point):
    """Exact complex center value, used only inside validated quadrature."""
    require_mode(k)
    a,b,c=validate_hat(hat)
    x,y=map(F,point)
    q=x*x+y*y
    if q<=a or q>=c:
        return F(0),F(0)
    re,im=F(1),F(0)
    for _ in range(k):
        re,im=re*x-im*y,re*y+im*x
    if q<=b:
        weight=1-(a/q)**k
    else:
        weight=(b**k-a**k)/(c**k-b**k)*((c/q)**k-1)
    return re*weight,im*weight


def lipschitz_upper(k,hat):
    """A global Euclidean Lipschitz constant for the complex-valued hat."""
    require_mode(k)
    a,b,c=validate_hat(hat)
    inner=k*sqrt_upper(b**(k-1))
    if a:
        inner+=k*sqrt_upper(a**(k-1))
    t=(b**k-a**k)/(c**k-b**k)
    outer=k*t*(c**k*sqrt_upper(b**(-k-1))+sqrt_upper(c**(k-1)))
    return max(inner,outer)


@lru_cache(maxsize=16)
def pi_interval(terms=14):
    """Machin's identity with alternating rational arctangent remainders."""
    if type(terms)is not int or terms<1:
        raise ValueError('positive integer arctangent length required')
    def atan(q):
        total=sum(((-1)**j*q**(2*j+1)/(2*j+1) for j in range(terms)),F(0))
        next_term=(-1)**terms*q**(2*terms+1)/(2*terms+1)
        return min(total,total+next_term),max(total,total+next_term)
    a,b=atan(F(1,5))
    c,d=atan(F(1,239))
    return 16*a-4*d,16*b-4*c


@lru_cache(maxsize=24)
def disk_grid(r2,grid,box_scale=64):
    """Square cells wholly inside, outside, or crossing a rational disk.

    Only retained cells are returned. Integer corner comparisons prove the
    classification; boundary-cell area is never silently treated as full.
    """
    r2=F(r2)
    if r2<=0 or type(grid)is not int or grid<2:
        raise ValueError('positive squared radius and integer grid>=2 required')
    radius=sqrt_upper(r2,box_scale)
    step=2*radius/grid
    area=step*step
    cell_radius_upper=sqrt_upper(F(2))*radius/grid
    threshold=r2*grid*grid/(radius*radius)
    cells=[]
    for i in range(grid):
        left,right=2*i-grid,2*i+2-grid
        minx=0 if left<=0<=right else min(left*left,right*right)
        maxx=max(left*left,right*right)
        for j in range(grid):
            low,high=2*j-grid,2*j+2-grid
            miny=0 if low<=0<=high else min(low*low,high*high)
            maxy=max(low*low,high*high)
            if minx+miny>=threshold:
                continue
            inside=maxx+maxy<=threshold
            cells.append((F(left+right,2)*radius/grid,
                          F(low+high,2)*radius/grid,inside))
    return radius,area,cell_radius_upper,tuple(cells)


def divide_interval(low,high,denlow,denhigh):
    if denlow<=0 or denhigh<denlow:
        raise ValueError('ordered positive denominator interval required')
    possibilities=(low/denlow,low/denhigh,high/denlow,high/denhigh)
    return min(possibilities),max(possibilities)


def disk_average_rectangle(k,hat,center,r2,grid=64):
    """Enclose mean Phi(center+u) over |u|²<=r2, including interface crossings.

    Returns Re/Im lower and upper endpoints, all rational. Refining the grid
    converges because the test is globally Lipschitz and the boundary cells
    have area tending to zero. This is not a harmonic center-value formula.
    """
    require_mode(k)
    hat=validate_hat(hat)
    center=tuple(map(F,center))
    r2=F(r2)
    if len(center)!=2 or r2<=0:
        raise ValueError('rational center pair and positive squared radius required')
    if center==(F(0),F(0)):
        return F(0),F(0),F(0),F(0)
    radius,area,cell_radius,cells=disk_grid(r2,grid)
    error=lipschitz_upper(k,hat)*cell_radius
    amplitude=sqrt_upper(hat[1]**k)
    re_sum_low=re_sum_high=im_sum_low=im_sum_high=0
    value_scale=10**12*grid
    for dx,dy,inside in cells:
        re,im=phi_value(k,hat,(center[0]+dx,center[1]+dy))
        rlo,rhi=max(-amplitude,re-error),min(amplitude,re+error)
        ilo,ihi=max(-amplitude,im-error),min(amplitude,im+error)
        if not inside:
            rlo,rhi=min(F(0),rlo),max(F(0),rhi)
            ilo,ihi=min(F(0),ilo),max(F(0),ihi)
        # A common outward mesh prevents denominator growth when thousands
        # of unrelated inverse-power values are summed. The rounding is
        # part of the certificate, not an inexact arithmetic shortcut.
        bounds=[value*value_scale for value in (rlo,rhi,ilo,ihi)]
        re_sum_low+=bounds[0].numerator//bounds[0].denominator
        re_sum_high+=-((-bounds[1].numerator)//bounds[1].denominator)
        im_sum_low+=bounds[2].numerator//bounds[2].denominator
        im_sum_high+=-((-bounds[3].numerator)//bounds[3].denominator)
    re_low,re_high=area*F(re_sum_low,value_scale),area*F(re_sum_high,value_scale)
    im_low,im_high=area*F(im_sum_low,value_scale),area*F(im_sum_high,value_scale)
    plo,phi=pi_interval(max(14,2*grid.bit_length()))
    re_low,re_high=divide_interval(re_low,re_high,plo*r2,phi*r2)
    im_low,im_high=divide_interval(im_low,im_high,plo*r2,phi*r2)
    if center[1]==0:
        im_low=im_high=F(0)
    if center[0]==0:
        if k%2:
            re_low=re_high=F(0)
        else:
            im_low=im_high=F(0)
    return re_low,re_high,im_low,im_high


def dual_certificate(k,hats,moment_rectangles):
    """Select rational midpoint dual coefficients, then validate the full box."""
    require_mode(k)
    hats=[validate_hat(hat) for hat in hats]
    rects=[tuple(map(F,rect)) for rect in moment_rectangles]
    if len(hats)!=len(rects) or not hats or any(len(rect)!=4 for rect in rects):
        raise ValueError('matching nonempty hats and four-endpoint rectangles required')
    if any(lo>hi or ilo>ihi for lo,hi,ilo,ihi in rects):
        raise ValueError('ordered real and imaginary intervals required')
    H=gram_matrix(k,hats)
    midpoint_re=[(lo+hi)/2 for lo,hi,ilo,ihi in rects]
    midpoint_im=[(ilo+ihi)/2 for lo,hi,ilo,ihi in rects]
    ar,ai=solve_positive(H,midpoint_re),solve_positive(H,midpoint_im)
    pairing=F(0)
    for r,i,(lo,hi,ilo,ihi) in zip(ar,ai,rects):
        pairing+=min(r*lo,r*hi)+min(i*ilo,i*ihi)
    norm=sum((H[j][ell]*(ar[j]*ar[ell]+ai[j]*ai[ell])
              for j in range(len(hats)) for ell in range(len(hats))),F(0))
    return {'k':k,'gram':H,'moment_rectangles':rects,'dual_real':ar,'dual_imaginary':ai,
            'angular_delta_lower':max(F(0),2*pairing-norm),
            'scope':'One validated dual over the complete moment box; overlapping same-mode tests are not added independently.'}


def shape_certificate(points,hats_by_mode,h=1,grid=64,radial_delta_lower=0):
    """Center rational coordinates, normalize Q(x)=n-1, and integrate all disks.

    hats_by_mode maps each k to (a,b,c) squared radii IN THE NORMALIZED x
    SCALE. Raw-coordinate rescaling cancels the common f^(k/2) factor in
    each Gram projection. A nonzero radial lower bound must be separately
    proved using the same smear h.
    """
    pts=[tuple(map(F,point)) for point in points]
    n=len(pts)
    h,radial=F(h),F(radial_delta_lower)
    if n<2 or any(len(point)!=2 for point in pts) or len(set(pts))!=n:
        raise ValueError('at least two distinct rational coordinate pairs required')
    if h<=0 or radial<0 or not hats_by_mode:
        raise ValueError('positive smear, nonnegative radial lower, and modes required')
    meanx=sum(point[0] for point in pts)/n
    meany=sum(point[1] for point in pts)/n
    v=[(a-meanx,b-meany) for a,b in pts]
    rawQ=sum(a*a+b*b for a,b in v)
    factor=F(n-1)/rawQ
    raw_r2=2*h*h/(n*factor)
    rows=[]
    total=F(0)
    for k,hats in sorted(hats_by_mode.items()):
        require_mode(k)
        raw_hats=[tuple(x/factor for x in validate_hat(hat)) for hat in hats]
        moments=[]
        for hat in raw_hats:
            rectangle=[F(0)]*4
            for center in v:
                value=disk_average_rectangle(k,hat,center,raw_r2,grid)
                rectangle=[a+b for a,b in zip(rectangle,value)]
            moments.append(tuple(rectangle))
        row=dual_certificate(k,raw_hats,moments)
        row['raw_hats']=raw_hats
        rows.append(row)
        total+=row['angular_delta_lower']
    loghlo,loghhi=log_interval(h)
    penalty=n*max(F(0),(h*h-1)/2-loghlo)
    return {'n':n,'normalization_factor_squared':factor,'raw_smear_squared':raw_r2,
            'grid':grid,'h':h,'angular_delta_lower':total,
            'radial_delta_lower':radial,'smear_penalty_upper':penalty,
            'delta_lower':max(F(0),total+radial-penalty),'modes':rows,
            'scope':'Rigorous configuration energy certificate including all crossing disks. No maximum or complete domain covering is inferred.'}


if __name__=='__main__':
    import json
    f=F(3,10)
    result=shape_certificate([(1,0),(-1,0),(2,0),(-2,0)],
        {2:[(0,f,9*f),(0,2*f,16*f)]},grid=32)
    for key in ('angular_delta_lower','delta_lower'):
        result[key]=outward_decimal((result[key],result[key]),10)[0]
    print(json.dumps({key:value for key,value in result.items() if key!='modes'},
                     default=str,indent=2))
