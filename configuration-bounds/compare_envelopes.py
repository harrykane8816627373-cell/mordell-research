"""Exact comparison after restoring all supplied attainable candidates.

The original polar family and saved integer candidates are source results.
The phase/centering refinement and the explicit twelve-point two-ring
family optimization are separate refinements. No floating-point decision.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
from baseline_finite_lower_bounds import gaussian_lower, central_polygon_lower
from triangular_lattice import (lattice_certificate, gaussian_norm,
                                triangular_norm, improved_centered_polygon)
from polar_refinement import ring_value, refined_polar, twelve_two_ring_maximum

ROOT = Path(__file__).resolve().parent


def require(ok, name):
    if not ok:
        raise ArithmeticError(name)


def verify_saved():
    payload = json.loads((ROOT/'saved_rational_bounds.json').read_text(encoding='utf-8'))
    result = {}
    for row in payload['rows']:
        n, points, D = row['n'], row['integer_points'], int(row['D'])
        require(len(points) == n and len(set(map(tuple,points))) == n, 'saved distinct coordinates')
        sx, sy = (sum(p[j] for p in points) for j in range(2))
        centered = [(n*x-sx,n*y-sy) for x,y in points]
        require(sum(x*x+y*y for x,y in centered) == n*D, 'saved exact prescribed moment')
        value = F(1)
        for (x,y),(u,v) in combinations(centered,2):
            value *= F((x-u)**2+(y-v)**2,D)
        require(value == F(row['value']), 'saved exact pair product')
        require(value > 0, 'saved positive value')
        result[n] = value
    require(sorted(result) == list(range(6,13)), 'source saved point-count range')
    require(result[10]/central_polygon_lower(10) > F(121,100), 'source ten-point stronger lower bound')
    return result


def digest(q):
    def ih(x):
        return hashlib.sha256(x.to_bytes(max(1,(x.bit_length()+7)//8),'big')).hexdigest()
    return {'numerator_sha256':ih(q.numerator),'denominator_sha256':ih(q.denominator)}


def audit(max_n=100):
    saved = verify_saved()
    moments = json.loads((ROOT/'moment_values.json').read_text(encoding='utf-8-sig'))
    rayleigh = {row['n']:F(row['rayleigh_lower']) for row in moments['rows']}
    rows=[]
    restored_polar_wins=[]
    saved_wins=[]
    phase_wins_previous_polar=[]
    final_wins=[]
    for n in range(2,max_n+1):
        candidates = {'regular_polygon':F(n**n),'central_polygon':central_polygon_lower(n),
                      'Gaussian_average':gaussian_lower(n),
                      'square_lattice':lattice_certificate(n,gaussian_norm)['value'],
                      'triangular_lattice':lattice_certificate(n,triangular_norm)['value']}
        old_perturbation=improved_centered_polygon(n)
        if old_perturbation is not None:
            candidates['old_explicit_perturbation']=old_perturbation
        if n in rayleigh:
            candidates['published_Rayleigh']=rayleigh[n]
        previous=max(candidates.values())
        original=ring_value(n)
        refined=refined_polar(n)['value']
        require(refined>=original and (n<7 or refined>original),'uniform polar refinement')
        restored=max(previous,original,saved.get(n,F(0)))
        best=max(restored,refined,twelve_two_ring_maximum() if n==12 else F(0))
        if original>previous:restored_polar_wins.append(n)
        if saved.get(n,F(0))>max(previous,original):saved_wins.append(n)
        if refined>max(previous,original):phase_wins_previous_polar.append(n)
        if best>restored:final_wins.append(n)
        row={'n':n,'original_polar_improves_previous_published':original>previous,
             'saved_candidate_improves_previous_and_polar':saved.get(n,F(0))>max(previous,original),
             'new_refinement_improves_complete_restored_envelope':best>restored,
             'previous':digest(previous),'restored':digest(restored),'final':digest(best)}
        if n<=30:
            row.update(previous_exact=str(previous),original_polar_exact=str(original),
                       refined_polar_exact=str(refined),saved_candidate_exact=str(saved[n]) if n in saved else None,
                       restored_exact=str(restored),final_exact=str(best))
        rows.append(row)
    return {'passed':True,'n_min':2,'n_max':max_n,'arithmetic':'integers and Fraction only',
            'restored_original_polar_improves_previous_at':restored_polar_wins,
            'restored_saved_candidate_improves_previous_and_polar_at':saved_wins,
            'phase_refinement_improves_previous_plus_polar_at':phase_wins_previous_polar,
            'new_refinement_improves_complete_restored_envelope_at':final_wins,
            'twelve_point_family_below_supplied_saved_candidate':twelve_two_ring_maximum()<saved[12],
            'scope':'Complete restored envelope includes source integer candidates n=6..12; no unknown global maximum is identified.',
            'rows':rows}


if __name__=='__main__':
    result=audit()
    (ROOT/'lower_bound_comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},ensure_ascii=False,indent=2))
