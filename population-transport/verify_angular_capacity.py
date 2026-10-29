"""Exact acceptance for circle/band capacities and integer population costs."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import sys

from angular_capacity import *

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'collective-stability'))
import ring_population as old_population
from block_stability import geometry,norm,mul,divide,ONE,scale,sub
sys.path.insert(0,str(ROOT.parent/'symmetry-breaking'))
from pair_stability import pair_slack


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def compositions(n):
    if n==0:
        yield ()
    for first in range(1,n+1):
        for tail in compositions(n-first):yield (first,)+tail


def run():
    composition_checks=0
    for n in range(2,17):
        brute={}
        for p in compositions(n):
            k=len(p)
            c=sum(m*(2*(m//2))**2 for m in p)
            strict=any(m>1 and m%2 for m in p)
            composition_checks+=1
            if k not in brute or c<brute[k][0]:brute[k]=(c,strict)
            elif c==brute[k][0]:brute[k]=(c,brute[k][1] and strict)
        for k,(cost,strict) in brute.items():
            result=minimum_cost(n,k)
            require(result['cost']==cost,'independent complete ordered population minimum')
            require(result['all_minimizers_strict']==strict,'strict-boundary flag includes every minimizer')
            require(sum(m*c for m,c in result['population_groups'])==n,'compressed population witness sum')
            require(sum(c for m,c in result['population_groups'])==k,'compressed layer witness count')
            require(sum(c*population_cost(m) for m,c in result['population_groups'])==cost,'compressed cost witness')
        costs=[brute[k][0] for k in range(1,n+1)]
        require(all(a>b for a,b in zip(costs,costs[1:])),'strict decrease under layer splitting')
    pair_exchange_checks=0
    for a in range(2,22,2):
        for b in range(2,22,2):
            candidates=[population_cost(a-1)+population_cost(b+1),population_cost(a+1)+population_cost(b-1)]
            require(F(sum(candidates),2)==a**3+b**3-2*((a-1)**2+(b-1)**2),
                    'exact even-to-odd exchange identity')
            require(min(candidates)<a**3+b**3,'two even populations cannot minimize')
            pair_exchange_checks+=1
    # Independent closed rational odd-skip sums for M=2,4,6,8.
    chebyshev_checks=0
    for width,b in ((F(0),F(1)),(F(1,50),F(101,100)),(F(1,10),F(2))):
        c=width*width
        b2=b*b
        expected={2:1/(c+4*b2),4:2/(c+2*b2),
                  6:2/(c+b2)+1/(c+4*b2),
                  8:4*(c+2*b2)/((c+2*b2)**2-2*b2*b2)}
        for M,value in expected.items():
            t,dt=chebyshev_value_derivative(M//2,1+c/(2*b2))
            require(dt/(2*b2*(t+1))==value,'independent closed-form odd-node sum')
            chebyshev_checks+=1
    for M in range(2,52,2):
        t,dt=chebyshev_value_derivative(M//2,1)
        require(t==1 and dt==F(M*M,4),'Chebyshev endpoint exact derivative')
        require(dt/(2*(t+1))==F(M*M,16),'circle reciprocal csc-square sum')
        chebyshev_checks+=1
    # Direct complex identity on unevenly spaced unit-circle rational points.
    unit=[]
    for t in (F(-3),F(-2),F(-1),F(-1,2),F(0),F(1,3),F(1),F(3)):
        unit.append(((1-t*t)/(1+t*t),2*t/(1+t*t)))
    phase_checks=0
    for z,w in combinations(unit,2):
        delta=sub(z,w)
        B=scale(divide(ONE,mul(delta,delta)),-1)
        adjusted=mul(mul(z,w),B)
        require(adjusted==(1/norm(delta),F(0)),'circle phase adjustment is positive reciprocal distance')
        phase_checks+=1
    # The pairwise thick-band phase inequality is tested independently.
    thick_pair_checks=0
    sep=F(16,49)
    a,b=F(99,100),F(101,100)
    eta=1-2*(b-a)**2/sep
    for ri in (a,F(1),b):
        for rj in (a,F(1),b):
            z=(ri,F(0))
            for u in unit:
                w=scale(u,rj)
                delta=sub(z,w)
                distance=norm(delta)
                if distance<sep:continue
                B=scale(divide(ONE,mul(delta,delta)),-1)
                adjusted=mul(u,B)
                require(adjusted[0]>=eta/distance,'direct thick-band adjusted real part')
                thick_pair_checks+=1
    convexity_checks=0
    for b2 in (F(1,2),F(1),F(3)):
        for width2 in (F(0),F(1,100),F(1,10),F(1,2)):
            for j in range(1,33):
                u=F(j,32)
                if 4*b2*u<width2:continue
                numerator=-width2+(12*b2+2*width2)*u-8*b2*u*u
                require(numerator>=-width2+4*b2*u>=0,'uniform convexity numerator bound')
                convexity_checks+=1

    points=[(70,0),(49,49),(0,70),(-49,49),(-70,0),(-49,-49),(0,-70),(49,-49)]
    data=geometry(points)
    factor=F(8)/data['Q']
    radii=[factor*norm(p) for p in data['points']]
    slacks=[]
    distances=[]
    for i,j in combinations(range(8),2):
        s=factor*norm(sub(data['points'][i],data['points'][j]))
        distances.append(s)
        slacks.append(pair_slack(8,radii[i],radii[j],s))
    require(min(slacks)>0 and min(distances)>F(4,7),'nonempty annulus witness passes old local pairs and global distance')
    require(max(radii)<F(7,4),'nonempty annulus witness passes old global radius')
    require(min(radii)>a*a and max(radii)<b*b,'normalised witness is strictly inside thick annulus')
    require(not old_population.population_certificate((4,4))['excluded_from_global_maxima'],
            'old exact-layer population budget permits witness populations')
    for q in radii:
        require(not layer_capacity(8,4,q)['excluded_from_local_maxima'],
                'new separate single-layer capacity alone permits both four-point layers')
    band=annulus_capacity(8,8,a,b)
    require(band['operator_lower']==F(407120106,104142017),'exact whole-annulus angular lower')
    require(band['stability_slack']==F(-85246093,208284034)<0,'entire annulus excludes every local maximum')
    require((1-a*a)**2==F(39601,10**8) and 1-a*a<b*b-1,
            'eight-point outlier implies stated strict radial variance lower')
    require(layer_capacity(8,8,1)['excluded_from_local_maxima'],'even circle bound detects eight points on one circle')
    require(layer_capacity(8,3,F(1,14))['excluded_from_local_maxima'],'odd population equality is excluded strictly')
    require(not layer_capacity(10,9,F(10,9))['excluded_from_local_maxima'],
            'criterion preserves known ten-point central-nine local maximum')
    examples=[]
    for n in (8,100,1000,10000):
        new=minimum_layers(n)
        old=old_population.minimum_layers(n)
        if new['last_excluded_layers']:
            r=new['last_excluded_cost']
            require(r['cost']>new['cost_cap'] or
                    (r['cost']==new['cost_cap'] and r['all_minimizers_strict']),
                    'last excluded integer layer count')
        require(new['minimum_layers']>=old['minimum_layers'],'new stability integer bound is no weaker at recorded examples')
        examples.append({'n':n,'new_local_necessary':new,'old_global_packing_necessary':old['minimum_layers'],
                         'elementary_layer_lower_endpoint':elementary_layer_lower(n)})
    invalid=[lambda:layer_capacity(8,9,1),lambda:minimum_cost(8,0),
             lambda:minimum_cost(100,3,max_cases=1),
             lambda:annulus_capacity(8,8,0,1),lambda:annulus_capacity(8,8,1,F(1,2)),
             lambda:chebyshev_value_derivative(2,F(1,2))]
    for call in invalid:
        try:call()
        except ValueError:pass
        else:raise ArithmeticError('invalid or incomplete certificate accepted')
    return {'passed':True,'arithmetic':'integers/Fraction; no floating-point acceptance',
            'independent_ordered_population_compositions':composition_checks,
            'even_to_odd_exchange_checks':pair_exchange_checks,
            'independent_chebyshev_sum_checks':chebyshev_checks,
            'circle_phase_identity_checks':phase_checks,'thick_pair_identity_checks':thick_pair_checks,
            'convexity_numerator_checks':convexity_checks,'rejected_invalid_cases':len(invalid),
            'eight_point_annulus_certificate':band,
            'eight_point_necessary_radial_variance_strict_lower':'39601/100000000',
            'nonempty_witness':{'integer_points':points,'normalization_squared':factor,
                'minimum_squared_distance':min(distances),'maximum_squared_radius':max(radii),
                'minimum_old_pair_slack':min(slacks)},
            'integer_layer_examples':examples,
            'scope':'The universal circle/band conditions are proved analytically. Finite checks validate formulas, exact population minimization through n=16 and the stated entire-annulus certificate. These are necessary stability conditions, not energies or optimum certification.'}


if __name__=='__main__':
    report=run()
    ROOT.joinpath('angular_capacity_verification.json').write_text(json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('nonempty_witness','integer_layer_examples')},default=str,ensure_ascii=False,indent=2))
