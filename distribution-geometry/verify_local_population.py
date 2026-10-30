"""Independent exact checks for translated energy and spatial populations."""
from fractions import Fraction as F
from pathlib import Path
import json

from local_population import (disk_mass_bounds,ball_population,annulus_population,
    spatial_grid_populations,multiple_hole_certificate,covering_radius_upper,
    background_grid,axis_distance_bounds,pi_interval,polar_gap_upper,
    uniform_disk_discrepancy,local_deficit_certificate)


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def run():
    exact_mass_checks=0
    for center,q,value in (((0,0),F(1),F(1,2)),((F(1,2),0),F(1,4),F(1,8)),
        ((0,0),F(2),F(1)),((0,0),F(9),F(1)),((3,0),F(1),F(0)),
        ((F(3,2),F(3,2)),F(0),F(0))):
        row=disk_mass_bounds(center,q)
        require(row['lower']==value==row['upper'],'exact contained/containing/disjoint mass')
        exact_mass_checks+=1
    # Equal radius sqrt(2), centers 0 and 2: lens area pi-2.
    pilo,pihi=pi_interval()
    true_low,true_high=F(1,2)-1/pilo,F(1,2)-1/pihi
    grid_rows=[]
    for grid in (16,32,64):
        row=disk_mass_bounds((2,0),2,grid)
        require(row['lower']<=true_low<=true_high<=row['upper'],'lens closed-form inside grid bounds')
        grid_rows.append(row)
    require(all(a['lower']<=b['lower'] and a['upper']>=b['upper']
                for a,b in zip(grid_rows,grid_rows[1:])),'nested complete grid improves area bounds')
    geometry_checks=0
    for center in ((0,0),(F(1,3),F(-2,5)),(2,0)):
        for left,right,bottom,top,bginside in background_grid(16):
            lx,hx=axis_distance_bounds(F(center[0]),left,right)
            ly,hy=axis_distance_bounds(F(center[1]),bottom,top)
            corners=[(x-center[0])**2+(y-center[1])**2
                     for x in (left,right) for y in (bottom,top)]
            require(max(corners)==hx+hy,'rectangle maximum at a corner')
            nearest_x=min(max(F(center[0]),left),right)
            nearest_y=min(max(F(center[1]),bottom),top)
            require((nearest_x-center[0])**2+(nearest_y-center[1])**2==lx+ly,
                    'rectangle closest-point projection')
            geometry_checks+=1
    # Triangle cost used for arbitrary translated CDF targets of slope <=1/2.
    triangle_checks=0
    for epsilon in (F(1,100),F(1,7),F(1,2),F(1)):
        end=2*epsilon
        integral=epsilon**2*end-epsilon*end**2/2+end**3/12
        require(integral==2*epsilon**3/3,'translated CDF triangle integral')
        for s in (F(0),F(1,2),F(2),F(10)):
            require(integral/(2*(s+2*epsilon))==epsilon**3/(3*(s+2*epsilon)),
                    'positive CDF error denominator')
            triangle_checks+=1
    hole_checks=0
    for b in (F(1,10),F(1,2),F(1)):
        require((b**4/4)/4==b**4/16,'uniform-background empty-ball fourth power')
        hole_checks+=1
    # Inscribed half-radius ball is valid even at a background boundary.
    for d in (F(0),F(1,2),F(1),F(7,5)):
        for rho in (F(1,10),F(1),F(5,2)):
            shifted=max(F(0),d-rho/2)
            require(d-shifted<=rho/2,'half-ball stays in original empty disk')
            require((shifted+rho/2)**2<=2,'half-ball stays in background disk')
            hole_checks+=1
    pair=multiple_hole_certificate(10000,
        [(F(-1,2),0,F(2401,10000)),(F(1,2),0,F(2401,10000))])
    require(pair['delta_lower']>F(4879681,8),'two arbitrary-center holes strengthen single-hole cost')
    require(pair['delta_lower']>polar_gap_upper(10000),'two-hole whole class excluded from global maxima')
    require(F(9999,10000)>F(99,100)**2,'normalized regular circle gives a nonempty two-hole class')
    offcenter=multiple_hole_certificate(10000,[(F(3,4),F(1,4),F(1,16))])
    require(offcenter['delta_lower']>0,'translated non-axis hole energy positive')
    outside=multiple_hole_certificate(10000,[(3,0,1)])
    require(outside['delta_lower']==0,'hole outside background has no automatic density penalty')
    deficits=[local_deficit_certificate(10000,(F(1,2),0),F(1,4),m)
              for m in (0,10,100,500,1000,10000)]
    require(all(a['delta_lower']>=b['delta_lower'] for a,b in zip(deficits,deficits[1:])),
            'less stringent population caps weaken the local deficit')
    require(deficits[0]['delta_lower']==multiple_hole_certificate(10000,[(F(1,2),0,F(1,4))])['delta_lower'],
            'zero-population deficit equals the empty-disk certificate')
    require(deficits[-1]['delta_lower']==0,'full population cap gives no artificial deficit')
    examples=[]
    discrepancy_examples=[]
    for n in (10000,1000000,100000000,1000000000000):
        row=covering_radius_upper(n,polar_gap_upper(n))
        require(row['radius_upper']>0,'positive covering-radius outer bound')
        examples.append(row)
        discrepancy=uniform_disk_discrepancy(n,polar_gap_upper(n),F(1))
        require(0<=discrepancy['x_error_upper']<=discrepancy['original_z_error_upper']<=1,
                'uniform translated-disk scale correction')
        discrepancy_examples.append(discrepancy)
    populations=[]
    n=1000000;gap=polar_gap_upper(n)
    for center in ((0,0),(1,0),(F(7,5),0)):
        row=ball_population(n,gap,center,F(1,16),64)
        require(0<=row['lower']<=row['upper']<=n,'translated integer population interval')
        populations.append(row)
    symmetric=ball_population(n,gap,(-1,0),F(1,16),64)
    require((symmetric['lower'],symmetric['upper'])==
            (populations[1]['lower'],populations[1]['upper']),'spatial reflection invariance')
    ring=annulus_population(n,gap,(F(1,2),0),F(1,16),F(1,4),32)
    require(0<=ring['lower']<=ring['upper']<=n,'arbitrary-center annulus interval')
    empty_ring=annulus_population(n,gap,(0,0),F(1,4),F(1,4))
    require(empty_ring['lower']==empty_ring['upper']==0,'zero-width annulus is empty')
    zero_ball=ball_population(n,gap,(0,0),0)
    require(zero_ball['upper']<=1,'zero-radius disk contains at most one distinct center')
    table=spatial_grid_populations(n,gap,[F(-1,2),F(0),F(1,2)],
                                  [F(-1,2),F(1,2)],F(1,16),32)
    require(len(table)==6,'complete stated spatial query grid')
    invalid=[lambda:disk_mass_bounds((0,0),-1),lambda:disk_mass_bounds((0,0),1,2,max_cells=1),
        lambda:ball_population(1,1,(0,0),1),lambda:ball_population(8,-1,(0,0),1),
        lambda:multiple_hole_certificate(10000,[(0,0,1),(0,0,1)]),
        lambda:annulus_population(8,1,(0,0),2,1),
        lambda:covering_radius_upper(8,-1),
        lambda:spatial_grid_populations(8,1,[0,1],[0,1],1,max_queries=3)]
    for call in invalid:
        try:
            call()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid input accepted')
    return {'passed':True,'arithmetic':'integer/Fraction, outward roots and certified pi',
        'exact_mass_checks':exact_mass_checks,'rectangle_geometry_checks':geometry_checks,
        'triangle_integral_checks':triangle_checks,'hole_and_halfball_checks':hole_checks,
        'rejected_invalid_cases':len(invalid),'lens_area_bounds':grid_rows,
        'two_hole_certificate':pair,'covering_radius_examples':examples,
        'local_population_deficit_examples':deficits,
        'uniform_disk_discrepancy_examples':discrepancy_examples,
        'translated_population_examples':populations,'annulus_example':ring,
        'listed_spatial_grid_query_count':len(table),
        'scope':'Finite arithmetic checks. LOCAL_POPULATION.md proves the local Gauss inequality, arbitrary-position populations, disjoint-hole additivity and whole-background covering theorem. No optimizer classification or global configuration grid is inferred.'}


if __name__=='__main__':
    report=run()
    Path(__file__).with_name('local_population_verification.json').write_text(
        json.dumps(report,default=str,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items()
        if k not in ('two_hole_certificate','translated_population_examples','annulus_example',
                     'local_population_deficit_examples')},
        default=str,indent=2))
