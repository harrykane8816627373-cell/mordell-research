"""One explicit ascent step versus the COMPLETE restored envelope.

All old saved candidates are retained. The new exact phase-family result
is also combined, so a weaker new candidate is never called a new best.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import hashlib
import json
import sys

from rational_ascent import ascent_step, value
from phase_graph import optimized_polar

ROOT=Path(__file__).resolve().parent


def require(ok,name):
    if not ok:raise ArithmeticError(name)


def first_points(n,triangular):
    q=lambda p:p[0]*p[0]+(p[0]*p[1]if triangular else 0)+p[1]*p[1]
    chosen=sorted(((a,b)for a in range(-n,n+1)for b in range(-n,n+1)),key=lambda p:(q(p),p))[:n]
    if triangular:return tuple((F(a)+F(b,2),F(b,2))for a,b in chosen)
    return tuple((F(a),F(b))for a,b in chosen)


def digest(q):
    def ih(x):return hashlib.sha256(x.to_bytes(max(1,(x.bit_length()+7)//8),'big')).hexdigest()
    return {'numerator_sha256':ih(q.numerator),'denominator_sha256':ih(q.denominator)}


def integer_witness(points):
    denominator=lcm(*(v.denominator for p in points for v in p))
    # A common positive real scale does not change the normalized product.
    return tuple((int(a*denominator),int(b*denominator))for a,b in points)


def audit():
    old=json.loads((ROOT/'old_envelope_n2_30.json').read_text(encoding='utf-8'))
    old_values={r['n']:F(r['value'])for r in old['rows']}
    saved=json.loads((ROOT/'saved_rational_bounds.json').read_text(encoding='utf-8'))
    saved_points={r['n']:r['integer_points']for r in saved['rows']}
    rows=[];ascent_old_wins=[];ascent_phase_wins=[];final_old_wins=[]
    counts={'square':0,'triangular':0,'saved':0}
    for n in range(2,31):
        seeds=[('square',1,first_points(n,False)),('triangular',3,first_points(n,True))]
        if n in saved_points:seeds.append(('saved',1,saved_points[n]))
        certificates=[];values=[]
        for name,d,points in seeds:
            step=ascent_step(points,d)
            witness=integer_witness(step['new_points'])
            require(value(witness,d)==step['after'],'integer witness reproduces updated value')
            if not step['stationary']:counts[name]+=1
            certificates.append({'seed':name,'d':d,'stationary':step['stationary'],
                                 'normalized_residual':str(step['normalized_residual']),
                                 'step':str(step['step']),'log_gain_lower':str(step['log_gain_lower']),
                                 'integer_witness':witness,'before':digest(step['before']),
                                 'after':digest(step['after'])})
            values.append((step['after'],name))
        best,seed=max(values,key=lambda p:p[0])
        phase=optimized_polar(n)['value'];previous=old_values[n]
        combined=max(previous,phase,best)
        ratio_gain_basis_points=(best/max(previous,phase)-1)*10000
        if best>previous:ascent_old_wins.append(n)
        if best>max(previous,phase):ascent_phase_wins.append(n)
        if combined>previous:final_old_wins.append(n)
        winner='old' if combined==previous else ('phase' if combined==phase else 'ascent_'+seed)
        rows.append({'n':n,'old':digest(previous),'phase':digest(phase),
                     'best_ascent':digest(best),'combined':digest(combined),
                     'winner':winner,'best_ascent_seed':seed,
                     'ascent_improves_old':best>previous,
                     'ascent_improves_old_and_new_phase':best>max(previous,phase),
                     'ascent_gain_basis_points_floor_over_old_and_phase':ratio_gain_basis_points.numerator//ratio_gain_basis_points.denominator,
                     'certificates':certificates})
    report={'passed':True,'n_min':2,'n_max':30,'ascent_steps_per_seed':1,
            'strict_seed_improvements':counts,
            'ascent_improves_complete_old_envelope_at':ascent_old_wins,
            'ascent_improves_complete_old_plus_new_phase_at':ascent_phase_wins,
            'combined_improves_complete_old_at':final_old_wins,
            'scope':'Only one explicit rational step per supplied seed; all saved candidates and new exact phase-family values retained. No global sharpness or global convergence is claimed.',
            'rows':rows}
    return report


if __name__=='__main__':
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    result=audit()
    (ROOT/'ascent_envelope_comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items()if k!='rows'},indent=2))
