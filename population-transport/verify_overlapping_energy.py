"""Exact replay: overlapping Green Gram and all crossing-disk contributions."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sys

from overlapping_energy import (gram_matrix,solve_positive,phi_value,lipschitz_upper,
    pi_interval,disk_grid,disk_average_rectangle,dual_certificate,shape_certificate,
    outward_decimal,log_interval)


def require(ok,label):
    if not ok:
        raise ArithmeticError(label)


def harmonic_coefficients(k,hat,q):
    a,b,c=map(F,hat)
    if q<=a or q>=c:
        return F(0),F(0)
    if q<=b:
        return F(1),-a**k
    t=(b**k-a**k)/(c**k-b**k)
    return -t,t*c**k


def independent_cross_integral(k,left,right):
    """Integrate piecewise powers directly, independently of the min kernel."""
    stops=sorted(set([F(0)]+list(map(F,left))+list(map(F,right))))
    total=F(0)
    for a,b in zip(stops,stops[1:]):
        ar,br=harmonic_coefficients(k,left,(a+b)/2)
        al,bl=harmonic_coefficients(k,right,(a+b)/2)
        total+=F(k,2)*ar*al*(b**k-a**k)
        if br*bl:
            require(a>0,'inverse-power piece avoids zero')
            total+=F(k,2)*br*bl*(a**(-k)-b**(-k))
    return total


def run(write_report=False):
    hats=[(F(0),F(1),F(9)),(F(0),F(2),F(16)),
          (F(1,4),F(3,2),F(7)),(F(9),F(10),F(20))]
    gram_checks=0
    for k in range(1,5):
        H=gram_matrix(k,hats)
        for i,left in enumerate(hats):
            for j,right in enumerate(hats):
                require(H[i][j]==independent_cross_integral(k,left,right),
                        'full overlapping Gram equals direct polar integration')
                gram_checks+=1
        for i in range(len(hats)):
            a,b,c=hats[i]
            require(H[i][i]==k*(b**k-a**k)*(c**k-a**k)/(c**k-b**k),
                    'diagonal equals established single-hat norm')
        require(H[0][3]==0,'disjoint same-mode supports are orthogonal')
        solve_positive(H,[1,2,3,4])

    H=gram_matrix(2,hats[:2])
    require(H==[[F(81,40),F(88,45)],[F(88,45),F(512,63)]],
            'overlapping example exact Green Gram')
    # A potential equal to the first test has measurements H[:,0].  Adding
    # individual projection energies would double count its correlated part.
    exact=dual_certificate(2,hats[:2],[(H[0][0],H[0][0],0,0),
                                    (H[1][0],H[1][0],0,0)])
    naive=H[0][0]+H[1][0]**2/H[1][1]
    require(exact['angular_delta_lower']==H[0][0]<naive,
            'full Gram prevents same-mode overlap double counting')

    plo,phi=pi_interval()
    require(F(314159,100000)<plo<phi<F(355,113),'Machin rational pi enclosure')
    grid_checks=0
    for r2,g in product((F(1,100),F(1,2),F(5,3)),(8,16,32)):
        radius,area,error,cells=disk_grid(r2,g)
        inside=sum(1 for x,y,status in cells if status)
        boundary=len(cells)-inside
        require(inside*area<=plo*r2 and (inside+boundary)*area>=phi*r2,
                'grid brackets the entire disk area')
        for x,y,status in cells:
            half=radius/g
            if status:
                require((abs(x)+half)**2+(abs(y)+half)**2<=r2,
                        'inside cells include all their corners')
            grid_checks+=1

    # On a whole harmonic piece the new integration encloses the old exact
    # mean, without assuming that the disk misses interfaces in general.
    center=(F(1,4),F(1,3))
    expected=phi_value(2,hats[0],center)
    for grid in (16,32):
        rect=disk_average_rectangle(2,hats[0],center,F(1,100),grid)
        require(rect[0]<=expected[0]<=rect[1] and rect[2]<=expected[1]<=rect[3],
                'validated disk mean encloses exact harmonic benchmark')
    require(disk_average_rectangle(3,hats[0],(0,0),F(5,3),16)==(0,0,0,0),
            'centered radial disk has zero angular mean')

    # Every disk in this actual four-point example crosses at least one
    # relevant hat interface. Both compact tests also overlap each other.
    f=F(3,10)
    actual_hats={2:[(0,f,9*f),(0,2*f,16*f)]}
    points=[(1,0),(-1,0),(2,0),(-2,0)]
    result32=shape_certificate(points,actual_hats,grid=32)
    result64=shape_certificate(points,actual_hats,grid=64)
    require(result64['normalization_factor_squared']==f
            and result64['raw_smear_squared']==F(5,3),'normalized crossing-disk geometry')
    require(result32['delta_lower']>F(1,10),'coarse crossing-disk certificate positive')
    require(result64['delta_lower']>F(1,3),'refined crossing-disk certificate exceeds one third')
    for old,new in zip(result32['modes'][0]['moment_rectangles'],
                       result64['modes'][0]['moment_rectangles']):
        require(new[1]-new[0]<old[1]-old[0],'refinement narrows measured moment interval')
    # Verify that the former center-value interface was genuinely unavailable.
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'collective-stability'))
    from local_energy import harmonic_region
    crossing_rejections=0
    for q in (F(1),F(4)):
        for hat in hats[:2]:
            try:
                harmonic_region(q,q,F(5,3),*hat)
            except ValueError:
                crossing_rejections+=1
    require(crossing_rejections==4,'old whole-harmonic-piece test rejects every stated disk/hat type')

    # A nonempty entire neighborhood class: centered y_i with Q(y)=10 and
    # |y_i-v_i|<=1/1000. Same normalization and same smear remain fixed.
    eta=F(1,1000)
    row=result64['modes'][0]
    expanded=[]
    for hat,rect in zip(row['raw_hats'],row['moment_rectangles']):
        deviation=4*lipschitz_upper(2,hat)*eta
        expanded.append((rect[0]-deviation,rect[1]+deviation,
                         rect[2]-deviation,rect[3]+deviation))
    neighborhood=dual_certificate(2,row['raw_hats'],expanded)
    require(neighborhood['angular_delta_lower']>F(1,3),
            'nonempty interface-crossing neighborhood excluded by energy')
    # Known attainable lower bound M4>=4^4; the disk budget is <1/8.
    l2,h2=log_interval(F(2))
    l3,h3=log_interval(F(3))
    budget_lo,budget_hi=14*l2-6*h3-3,14*h2-6*l3-3
    require(0<budget_lo<=budget_hi<F(1,8)<F(1,3),
            'entire neighborhood contains no global maximum')

    invalid_calls=[
        lambda:gram_matrix(0,hats),lambda:gram_matrix(1,[(1,1,2)]),
        lambda:dual_certificate(2,[hats[0],hats[0]],[(1,1,0,0),(1,1,0,0)]),
        lambda:disk_average_rectangle(2,hats[0],(1,0),0,16),
        lambda:disk_average_rectangle(2,hats[0],(1,0),1,1),
        lambda:dual_certificate(2,[hats[0]],[(2,1,0,0)]),
        lambda:shape_certificate([(0,0),(0,0)],{1:[hats[0]]}),
    ]
    for callback in invalid_calls:
        try:
            callback()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid overlapping-energy input accepted')

    report={
        'status':'PASS','arithmetic':'Fraction, integer disk geometry, outward value rounding and Machin pi intervals',
        'overlapping_gram_integral_checks':gram_checks,'classified_grid_cells_checked':grid_checks,
        'rejected_invalid_cases':len(invalid_calls),'old_interface_rejections':crossing_rejections,
        'example_points':points,'example_actual_squared_hats':actual_hats[2],
        'example_raw_smear_squared':str(result64['raw_smear_squared']),
        'example_raw_gram':[[str(x) for x in line] for line in H],
        'grid32_delta_lower_outward':[str(x) for x in outward_decimal((result32['delta_lower'],result32['delta_lower']))],
        'grid64_delta_lower_outward':[str(x) for x in outward_decimal((result64['delta_lower'],result64['delta_lower']))],
        'grid64_raw_moment_rectangles':[[str(x) for x in outward_decimal(rect[:2])]+[str(x) for x in outward_decimal(rect[2:])]
                                        for rect in result64['modes'][0]['moment_rectangles']],
        'neighborhood_radius_raw':'1/1000','neighborhood_total_raw_moment':'10',
        'neighborhood_delta_lower_outward':[str(x) for x in outward_decimal((neighborhood['angular_delta_lower'],neighborhood['angular_delta_lower']))],
        'existing_global_optimum_gap_upper':'1/8',
        'scope':'Finite exact implementation tests and an interface-crossing nonempty neighborhood class. The example is not claimed to evade older geometric exclusions; no unknown sharp value or full-domain cover is established.',
    }
    if write_report:
        Path(__file__).with_name('overlapping_energy_verification.json').write_text(
            json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,ensure_ascii=False,indent=2))
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    run(parser.parse_args().write_report)
