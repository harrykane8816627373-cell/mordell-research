"""Translated local Gauss energy, empty disks, and spatial population bounds.

All radii/centers use Coulomb x scale Q(x)=n-1. See LOCAL_POPULATION.md.
No point configuration, optimum, or full configuration-space grid is inferred.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'configuration-bounds'))
sys.path.insert(0,str(ROOT/'radial-angular'))
sys.path.insert(0,str(ROOT/'collective-stability'))
from radial_distribution import root_upper,polar_gap_upper
from radial_layers import smearing_penalty_upper,log_bounds
from ring_population import pi_interval


def require_n(n):
    if type(n)is not int or n<2:
        raise ValueError('integer n>=2 required')


def sqrt_lower(q,scale=10**12):
    q=F(q)
    if q<0 or type(scale)is not int or scale<1:
        raise ValueError('nonnegative squared input and positive root scale required')
    return F(isqrt(q.numerator*scale*scale//q.denominator),scale)


def center_pair(center):
    result=tuple(map(F,center))
    if len(result)!=2:
        raise ValueError('rational two-coordinate center required')
    return result


def disk_mass_bounds(center,squared_radius,grid=64,max_cells=1000000):
    """Background mass of B(center,sqrt(q)); exact cases or complete grid.

    Background density is 1/(2*pi) on |w|^2<=2. Boundary-cell area is
    UNKNOWN and therefore counted only in the upper bound.
    """
    c=center_pair(center);q=F(squared_radius)
    if q<0 or type(grid)is not int or grid<2 or type(max_cells)is not int or max_cells<1:
        raise ValueError('q>=0, integer grid>=2 and positive cell budget required')
    if grid*grid>max_cells:
        raise ValueError('complete background grid exceeds the stated cell budget')
    d2=c[0]*c[0]+c[1]*c[1]
    def exact(value,method):
        return {'lower':value,'upper':value,'method':method}
    if q==0:
        return exact(F(0),'zero-area query')
    v=2-d2-q
    if v>=0 and v*v>=4*d2*q:
        return exact(q/2,'query disk contained in background')
    v=q-d2-2
    if v>=0 and v*v>=8*d2:
        return exact(F(1),'query disk contains background')
    v=d2-2-q
    if v>=0 and v*v>=8*q:
        return exact(F(0),'disjoint disks')
    inside=possible=0
    for left,right,bottom,top,background_inside in background_grid(grid):
        mnx,mxx=axis_distance_bounds(c[0],left,right)
        mny,mxy=axis_distance_bounds(c[1],bottom,top)
        if mnx+mny>=q:
            continue
        possible+=1
        if background_inside and mxx+mxy<=q:
            inside+=1
    area=F(16,grid*grid)
    pilo,pihi=pi_interval()
    return {'lower':inside*area/(2*pihi),
            'upper':min(F(1),possible*area/(2*pilo)),
            'method':'complete rational two-disk intersection grid',
            'grid':grid,'inside_cells':inside,'possible_cells':possible,'cell_area':area}


def axis_distance_bounds(center,left,right):
    low=F(0) if left<=center<=right else min((left-center)**2,(right-center)**2)
    return low,max((left-center)**2,(right-center)**2)


@lru_cache(maxsize=16)
def background_grid(grid):
    cells=[]
    edges=[F(-2)+F(4*j,grid) for j in range(grid+1)]
    axes=[axis_distance_bounds(F(0),a,b) for a,b in zip(edges,edges[1:])]
    for i,(left,right) in enumerate(zip(edges,edges[1:])):
        for j,(bottom,top) in enumerate(zip(edges,edges[1:])):
            if axes[i][0]+axes[j][0]>=2:
                continue
            cells.append((left,right,bottom,top,axes[i][1]+axes[j][1]<=2))
    return tuple(cells)


def ball_population(n,proven_gap_upper,center,squared_radius,grid=64,scale=10**12):
    """Necessary integer population in an arbitrary translated closed disk.

    The supplied gap upper must independently bound delta=log(U/Delta).
    No angular-energy subtraction about a different center is performed.
    """
    require_n(n)
    gap,q=F(proven_gap_upper),F(squared_radius)
    c=center_pair(center)
    if gap<0 or q<0:
        raise ValueError('nonnegative proven gap and query squared radius required')
    r2=F(2,n)
    cross=root_upper(q*r2,2,scale)
    qplus=q+r2+2*cross
    qminus=max(F(0),q+r2-2*cross) if q>=r2 else F(0)
    error=min(F(1),root_upper(3*(qplus+2)*gap/(n*n),3,scale))
    inner=disk_mass_bounds(c,qminus,grid)
    outer=disk_mass_bounds(c,qplus,grid)
    low=n*max(F(0),inner['lower']-error)
    high=n*min(F(1),outer['upper']+error)
    lo=-((-low.numerator)//low.denominator)
    hi=high.numerator//high.denominator
    if q==0:
        hi=min(hi,1)  # Finite discriminant implies distinct point centers.
    return {'n':n,'center':c,'squared_radius':q,'lower':lo,'upper':hi,
            'CDF_error_upper':error,'inner_squared_radius_lower':qminus,
            'outer_squared_radius_upper':qplus,'background_inner':inner,'background_outer':outer,
            'excluded_by_bounds':lo>hi,
            'scope':'Counts |x_i-center|^2<=q in Coulomb scale; requires the supplied independently proved energy-gap upper.'}


def annulus_population(n,proven_gap_upper,center,inner_squared,outer_squared,grid=64,scale=10**12):
    a,b=F(inner_squared),F(outer_squared)
    if not 0<=a<=b:
        raise ValueError('0<=inner squared radius<=outer squared radius required')
    if a==b:
        require_n(n)
        if F(proven_gap_upper)<0:
            raise ValueError('nonnegative independently proved gap required')
        return {'n':n,'center':center_pair(center),'inner_squared':a,'outer_squared':b,
                'lower':0,'upper':0,'scope':'The stated annulus is empty.'}
    inner=ball_population(n,proven_gap_upper,center,a,grid,scale)
    outer=ball_population(n,proven_gap_upper,center,b,grid,scale)
    return {'n':n,'center':center_pair(center),'inner_squared':a,'outer_squared':b,
            'lower':max(0,outer['lower']-inner['upper']),
            'upper':min(n,outer['upper']-inner['lower']),
            'scope':'Counts inner_squared<|x_i-center|^2<=outer_squared; separate disk bounds are combined safely.'}


def spatial_grid_populations(n,proven_gap_upper,x_values,y_values,squared_radius,grid=64,max_queries=256):
    xs,ys=tuple(map(F,x_values)),tuple(map(F,y_values))
    if not xs or not ys or type(max_queries)is not int or max_queries<1 or len(xs)*len(ys)>max_queries:
        raise ValueError('nonempty coordinate lists within the complete query budget required')
    return [ball_population(n,proven_gap_upper,(x,y),squared_radius,grid) for y in ys for x in xs]


def uniform_disk_discrepancy(n,proven_gap_upper,radius_limit,scale=10**12):
    """Uniform over EVERY center and every query radius <=radius_limit.

    The x error compares with the radius-sqrt(2) uniform background.
    Original-z error adds the exact 1/n dilation discrepancy.
    """
    require_n(n)
    gap,T=F(proven_gap_upper),F(radius_limit)
    if gap<0 or T<0:
        raise ValueError('nonnegative proved energy gap and query radius limit required')
    r=root_upper(F(2,n),2,scale)
    smooth_error=root_upper(3*((T+r)**2+2)*gap/(n*n),3,scale)
    xerror=min(F(1),smooth_error+T*r+F(1,n))
    return {'n':n,'radius_limit':T,'x_error_upper':xerror,
            'original_z_error_upper':min(F(1),xerror+F(1,n)),
            'scope':'Uniform mass discrepancy over all centers and radii <=radius_limit. Original z uses its natural radius-sqrt(2) background and the 1/n scale correction.'}


def multiple_hole_certificate(n,holes,h=1,scale=10**12):
    """Sum disjoint local spatial energy balls arising from empty point disks.

    holes=(center_x,center_y,empty_squared_radius). Each OPEN point disk
    must independently be empty. Returned shrunken balls must be disjoint.
    """
    require_n(n)
    holes=list(holes)
    h=F(h)
    if h<=0 or not holes:
        raise ValueError('positive smear and nonempty hole list required')
    r=root_upper(2*h*h/n,2,scale)
    R0=sqrt_lower(2,scale)
    rows=[]
    for cx,cy,q in holes:
        c=center_pair((cx,cy));q=F(q)
        if q<0:
            raise ValueError('nonnegative empty squared radius required')
        point_clearance=max(F(0),sqrt_lower(q,scale)-r)
        background_clearance=max(F(0),R0-root_upper(c[0]**2+c[1]**2,2,scale))
        b=min(point_clearance,background_clearance)
        rows.append({'center':c,'empty_squared_radius':q,'energy_ball_radius':b})
    for i,row in enumerate(rows):
        for other in rows[:i]:
            if row['energy_ball_radius']==0 or other['energy_ball_radius']==0:
                continue
            d2=sum((a-b)**2 for a,b in zip(row['center'],other['center']))
            if d2<(row['energy_ball_radius']+other['energy_ball_radius'])**2:
                raise ValueError('energy balls overlap; local energy cannot be added independently')
    energy=sum((row['energy_ball_radius']**4/16 for row in rows),F(0))
    penalty=n*smearing_penalty_upper(h)
    return {'n':n,'h':h,'holes':rows,'local_energy_lower':energy,
            'raw_delta_lower':n*n*energy-penalty,
            'delta_lower':max(F(0),n*n*energy-penalty),'smear_penalty_upper':penalty,
            'scope':'Conditional on all open point disks being empty; disjoint spatial Dirichlet regions are added, not overlapping whole-energy bounds.'}


def local_deficit_certificate(n,center,squared_radius,population_upper,h=1,scale=10**12):
    """Local energy if a translated point disk contains at most m centers.

    Uses the exact log integral after shrinking the disk and clipping it
    to the uniform-background interior. Empty disks are the case m=0.
    """
    require_n(n)
    c=center_pair(center);q,h=F(squared_radius),F(h)
    m=population_upper
    if q<0 or h<=0 or type(m)is not int or not 0<=m<=n:
        raise ValueError('q>=0, h>0 and integer population upper in [0,n] required')
    r=root_upper(2*h*h/n,2,scale)
    b=min(max(F(0),sqrt_lower(q,scale)-r),
          max(F(0),sqrt_lower(2,scale)-root_upper(c[0]**2+c[1]**2,2,scale)))
    mass=F(m,n)
    if b*b<=2*mass:
        energy=F(0)
    elif m==0:
        energy=b**4/16
    else:
        loglo,loghi=log_bounds(b*b/(2*mass))
        energy=max(F(0),b**4/16-mass*b*b/2+mass*mass*loglo/2+3*mass*mass/4)
    penalty=n*smearing_penalty_upper(h)
    return {'n':n,'center':c,'squared_radius':q,'population_upper':m,
            'energy_ball_radius':b,'local_energy_lower':energy,
            'delta_lower':max(F(0),n*n*energy-penalty),
            'scope':'Conditional on at most the supplied number of point centers in the disk. This is one local radial energy term; do not add the full empty-disk cost again.'}


def covering_radius_upper(n,proven_gap_upper,scale=10**12):
    """Covers the ENTIRE closed background disk |center|<=sqrt(2).

    Every such center has some x_i at distance <= returned radius. Uses
    the inscribed half-radius ball even when the hole touches the boundary.
    """
    require_n(n)
    gap=F(proven_gap_upper)
    if gap<0:
        raise ValueError('nonnegative independently proved gap required')
    energetic=2*root_upper(F(2,n),2,scale)+4*root_upper(gap/(n*n),4,scale)
    trivial=root_upper(2,2,scale)+root_upper(F(n-1,n),2,scale)
    return {'n':n,'radius_upper':min(energetic,trivial),
            'energy_radius_upper':energetic,'moment_radius_upper':trivial,
            'scope':'Uniform nearest-point radius for every center in |center|<=sqrt(2), Coulomb x scale; not a separation bound.'}


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('n',type=int)
    args=p.parse_args()
    print(json.dumps(covering_radius_upper(args.n,polar_gap_upper(args.n)),default=str,indent=2))
