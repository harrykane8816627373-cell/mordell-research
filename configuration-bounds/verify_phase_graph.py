"""Independent sign enumeration and comparison to the published envelope.

Optional --baseline points to the source-review directory with its existing
standard-library baseline implementation/data. No network is used.
"""
from fractions import Fraction as F
from itertools import combinations,product
from math import prod
from pathlib import Path
import argparse,hashlib,json,sys
from phase_graph import (polar_data,ring_product,phase_graph,graph_gain,
                         optimize_signs,optimized_polar,v2)


def require(ok,name):
    if not ok:raise ArithmeticError(name)


def digest(q):
    def h(a):return hashlib.sha256(a.to_bytes(max(1,(a.bit_length()+7)//8),'big')).hexdigest()
    return {'numerator_sha256':h(q.numerator),'denominator_sha256':h(q.denominator)}


def audit(baseline=None,max_n=100):
    parity_checks=0;gain_checks=0;full_checks=0;states=0
    # Check the sign decomposition independently for arbitrary population
    # pairs, not merely the sixfold polar family.
    for mi in range(1,17):
        for mj in range(1,17):
            from math import gcd
            g=gcd(mi,mj);L=mi*mj//g
            if L%2:continue
            rings=((mi,2),(mj,5))
            groups,unary,edges=phase_graph(rings)
            base=ring_product(mi+mj,rings,(1,1))
            for signs in product((-1,1),repeat=2):
                direct=ring_product(mi+mj,rings,signs)/base
                require(direct==graph_gain(signs,unary,edges),'generic pair decomposition')
                parity_checks+=1
    rows=[];phase_wins=[];envelope_wins=[]
    if baseline is not None:
        baseline=Path(baseline).resolve();sys.path.insert(0,str(baseline))
        from compare_envelopes import verify_saved
        from baseline_finite_lower_bounds import gaussian_lower,central_polygon_lower
        from triangular_lattice import lattice_certificate,gaussian_norm,triangular_norm,improved_centered_polygon
        from polar_refinement import refined_polar,twelve_two_ring_maximum
        saved=verify_saved()
        moments=json.loads((baseline/'moment_values.json').read_text(encoding='utf-8-sig'))
        rayleigh={r['n']:F(r['rayleigh_lower']) for r in moments['rows']}
    for n in range(2,max_n+1):
        n0,r,rings=polar_data(n);result=optimized_polar(n)
        require(result['value']>=result['previous_two_sign_value'],'old two-sign comparison')
        if result['value']>result['previous_two_sign_value']:phase_wins.append(n)
        groups,unary,edges=phase_graph(rings)
        base=ring_product(n,rings,(1,)*len(rings),True)
        # Full independent enumeration is manageable for these actual
        # small ring counts. It proves each finite output is a full family
        # maximum independently of the block solver.
        best=F(0)
        for signs in product((-1,1),repeat=len(rings)):
            direct=ring_product(n,rings,signs,True)
            require(direct==base*graph_gain(signs,unary,edges),'all-sign direct product')
            best=max(best,direct);gain_checks+=1
        require(best==result['value'],'block versus independent complete sign enumeration')
        full_checks+=1;states+=result['states']
        row={'n':n,'rings':rings,'signs':result['signs'],'block_states':result['states'],
             'full_family_states':1<<len(rings),'class_sizes':result['class_sizes'],
             'strictly_above_previous_two_signs':result['value']>result['previous_two_sign_value'],
             'value':digest(result['value'])}
        if n<=30:row['value_exact']=str(result['value'])
        if baseline is not None:
            previous=max(F(n**n),central_polygon_lower(n),gaussian_lower(n),
                         lattice_certificate(n,gaussian_norm)['value'],
                         lattice_certificate(n,triangular_norm)['value'],
                         refined_polar(n)['value'],saved.get(n,F(0)),rayleigh.get(n,F(0)),
                         twelve_two_ring_maximum() if n==12 else F(0))
            old_perturbation=improved_centered_polygon(n)
            if old_perturbation is not None:previous=max(previous,old_perturbation)
            row['strictly_above_complete_published_envelope']=result['value']>previous
            row['previous_envelope']=digest(previous)
            if result['value']>previous:envelope_wins.append(n)
            if n<=30:row['previous_envelope_exact']=str(previous)
        rows.append(row)
    # Budget refusals must occur before any incomplete result is returned.
    thirty=optimized_polar(30)
    require(thirty['value']/thirty['previous_two_sign_value']==F(90155,84693)**12,'thirty-point exact factor')
    require(F(90155,84693)**12>2,'thirty-point gain strictly above two')
    try:optimize_signs(polar_data(100)[2],max_states=1)
    except ValueError as e:require('above budget' in str(e),'budget refusal detail')
    else:raise ArithmeticError('budget refusal missing')
    return {'passed':True,'n_min':2,'n_max':max_n,'generic_pair_sign_checks':parity_checks,
            'direct_sign_product_checks':gain_checks,'complete_family_maxima_checked':full_checks,
            'block_states_enumerated':states,'strict_two_sign_improvements':phase_wins,
            'strict_complete_envelope_improvements':envelope_wins if baseline is not None else None,
            'scope':'Exact finite +/- family optimization; not a global planar maximum or optimization over arbitrary phases/radii.',
            'rows':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline',type=Path)
    p.add_argument('--max-n',type=int,default=100);args=p.parse_args()
    report=audit(args.baseline,args.max_n)
    (Path(__file__).parent/'phase_graph_comparison.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},ensure_ascii=False,indent=2))
