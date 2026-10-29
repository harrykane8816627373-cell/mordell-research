"""Exact acceptance of saved integer candidates against the full prior envelope.

No NumPy, optimizer, float comparison, network, or optional library is used.
The floating search history is provenance only and is never an input to a
mathematical decision. Coordinate products are recomputed from scratch.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys

from exact_bounds import integer_value, bounded_log

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'radial-angular'))
from verify_population_search import (LEGACY, PRIOR, prior_ascent_values, contracted_polar,
    optimized_polar, verify_saved, gaussian_lower, central_polygon_lower, lattice_certificate,
    gaussian_norm, triangular_norm, improved_centered_polygon, refined_polar,
    twelve_two_ring_maximum, digest, population_search)
from subset_moment import log_interval, disk_log_interval, outward_decimal
from radial_layers import minimum_layers
from rational_ascent import value as independent_point_value


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def previous_values(point_counts):
    """Rebuild all earlier terms, retaining each previously certified gain."""
    saved=verify_saved()
    moments=json.loads((LEGACY/'moment_values.json').read_text(encoding='utf-8-sig'))
    rayleigh={r['n']:F(r['rayleigh_lower']) for r in moments['rows']}
    ascents, count=prior_ascent_values()
    prior_report=json.loads((ROOT/'radial-angular/population_search_verification.json').read_text(encoding='utf-8'))
    prior_rows={r['n']:r for r in prior_report['rows']}
    values={}
    for n in point_counts:
        base=max(F(n**n),central_polygon_lower(n),gaussian_lower(n),
            lattice_certificate(n,gaussian_norm)['value'], lattice_certificate(n,triangular_norm)['value'],
            refined_polar(n)['value'],saved.get(n,F(0)),rayleigh.get(n,F(0)),
            twelve_two_ring_maximum() if n==12 else F(0), optimized_polar(n)['value'],
            ascents.get(n,F(0)))
        perturbation=improved_centered_polygon(n)
        if perturbation is not None:
            base=max(base,perturbation)
        if n>=13:
            base=max(base,contracted_polar(n)['value'])
        require(digest(base)==prior_rows[n]['previous_complete'],'full pre-population envelope digest')
        population=population_search(n)
        require(digest(population['value'])==prior_rows[n]['actual'],'prior population witness digest')
        require(digest(population['surrogate'])==prior_rows[n]['surrogate'],'prior population proxy digest')
        values[n]=max(base,population['lower_bound'])
    return values,count


def run(path, partial=False):
    witness=json.loads(path.read_text(encoding='utf-8'))
    counts=[r['n'] for r in witness['rows']]
    require(counts==sorted(set(counts)) and counts and all(10<=n<=100 for n in counts),'unique valid point counts')
    if not partial:
        require(counts==list(range(10,101)),'complete claimed 10..100 witness range')
    # Independent Fraction coordinate arithmetic checks centering and scaling.
    examples=[[(0,0),(1,0)],[(2,3),(-1,7),(4,-2)],
              [(0,0),(1,0),(0,1),(-1,0),(0,-1)]]
    for points in examples:
        require(integer_value(points)==independent_point_value(points),'independent centered distance product')
        require(integer_value([(7*x+5,7*y-3) for x,y in points])==integer_value(points),'similarity invariance')
    log_checks=0
    for q in (F(1),F(2),F(1,2),F(3,7),F(19,11),F(2)**100,
              F(2)**100-F(1,2**12),F(2)**100+F(1,2**12),F(10)**400):
        lo,hi=bounded_log(q,bits=32,terms=20)
        a,b=log_interval(q,terms=28)
        require(lo<=b and a<=hi,'independent logarithm interval overlap')
        require(hi-lo<F(1,10**8),'documented small interval width')
        log_checks+=1
    old,old_count=previous_values(counts)
    rows=[]
    improvements=[]
    stronger_layers=[]
    for source in witness['rows']:
        n=source['n']
        points=source['integer_points']
        require(len(points)==n,'exact witness point count')
        value=integer_value(points)
        # For small witnesses also use the older coordinate implementation.
        if n<=14:
            require(value==independent_point_value(points),'independent actual coordinate product')
        ratio=value/old[n]
        if ratio>1:
            improvements.append(n)
        delta_log=bounded_log(ratio)
        log_value=bounded_log(value)
        disk_log=disk_log_interval(n)
        gap=outward_decimal((disk_log[0]-log_value[1],disk_log[1]-log_value[0]),12)
        require(gap[0]>=0,'candidate below proven disk upper bound')
        combined=max(value,old[n])
        combined_log=bounded_log(combined)
        combined_gap=outward_decimal((disk_log[0]-combined_log[1],disk_log[1]-combined_log[0]),12)
        old_log=bounded_log(old[n])
        old_gap_upper=outward_decimal((disk_log[0]-old_log[1],disk_log[1]-old_log[0]),12)[1]
        old_layers=minimum_layers(n,old_gap_upper)
        new_layers=minimum_layers(n,combined_gap[1])
        if new_layers['minimum_layers']>old_layers['minimum_layers']:
            stronger_layers.append(n)
        improvement=(ratio-1)*10000
        rows.append({'n':n,'actual':digest(value),'complete_previous':digest(old[n]),
            'strictly_improves_complete_previous':ratio>1,
            'relative_gain_basis_points_floor':improvement.numerator//improvement.denominator,
            'log_gain_interval':outward_decimal(delta_log,12),
            'log_U_over_candidate_interval':gap,
            'global_optimizer_gap_upper':combined_gap[1],
            'previous_minimum_layers':old_layers['minimum_layers'],
            'new_minimum_layers':new_layers['minimum_layers'],
            'last_excluded_layer_certificate':new_layers['last_excluded_certificate']})
    return {'passed':True,'point_counts':counts,'independent_small_coordinate_checks':len(examples),
        'independent_actual_coordinate_checks':sum(n<=14 for n in counts),
        'independent_log_checks':log_checks,'prior_ascent_witnesses_recomputed':old_count,
        'strict_complete_envelope_improvements':improvements,
        'stronger_necessary_layer_counts_at':stronger_layers,
        'scope':'Integer coordinates prove attainable lower bounds only. Gap upper bounds apply to every global optimizer; layer counts are necessary, not exact. No floating search stopping flag certifies stationarity or global optimality.',
        'rows':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--partial',action='store_true')
    p.add_argument('--input',type=Path,default=HERE/'integer_witnesses.json')
    p.add_argument('--output',type=Path,default=HERE/'witness_verification.json')
    args=p.parse_args()
    if hasattr(sys,'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    report=run(args.input,args.partial)
    args.output.write_text(json.dumps(report,ensure_ascii=False,default=str,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},ensure_ascii=False,indent=2))
