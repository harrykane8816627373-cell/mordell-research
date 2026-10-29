"""Independent exact checks of the polynomial stationarity certificates."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

from stationarity import (multiply, binomial_remainder, regular_component_possible,
                          ring_stationarity)
from exact_bounds import bounded_log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'radial-angular'))
sys.path.insert(0, str(ROOT/'configuration-bounds'))
from population_search import population_search, ring_value
from rational_ascent import residual_data, ascent_step, norm, sub
from subset_moment import outward_decimal


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def digest_vector(v):
    return hashlib.sha256(json.dumps([str(x) for x in v], separators=(',', ':')).encode()).hexdigest()


def run():
    # Polynomial long division, independent of the coefficient-folding code.
    remainders = 0
    for m in range(2, 12):
        for c in (F(-3, 2), F(2, 3), F(5)):
            poly = [F((k*k+3*k)%11-5, k+1) for k in range(3*m+4)]
            long_division = poly[:]
            for k in range(len(poly)-1, m-1, -1):
                long_division[k-m] += c*long_division[k]
                long_division[k] = 0
            require(long_division[:m] == binomial_remainder(poly, m, c), 'independent long division')
            remainders += 1
    # Rational-coordinate rings permit comparison with all pairwise forces.
    examples = [(((2, 1),), (1,), 0, [(1,0),(-1,0)]),
        (((2, 1), (2, 4)), (1, -1), 0, [(1,0),(-1,0),(0,2),(0,-2)]),
        (((4, 1), (4, 4)), (1, 1), 0,
         [(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)]),
        (((4, 1), (4, 4)), (1, 1), 1,
         [(0,0),(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)])]
    for rings, signs, center, points in examples:
        certificate = ring_stationarity(rings, signs, center)
        exact = residual_data(points)
        require(certificate['stationary'] == (exact['R'] == 0), 'direct force zero test')
        require(certificate['residual_norm_squared_lower'] <= exact['R'], 'direct residual lower bound')
        actual_separation = min(norm(sub(p,q), 1) for i,p in enumerate(exact['points'])
                                for q in exact['points'][i+1:])
        require(certificate['pair_separation_squared_lower'] <= actual_separation, 'direct separation')
        update = ascent_step(points)
        require(update['log_gain_lower'] >= certificate['log_global_gap_lower'], 'independent ascent gain')
    for m in range(2, 16):
        for M in range(0, 21):
            require(regular_component_possible(m, M) == (M == 0 or m <= M), 'population domain')
    rows = []
    stationary = []
    majority_excluded = []
    for n in range(2, 101):
        candidate = population_search(n)
        check = ring_stationarity(candidate['rings'], candidate['signs'], candidate['center'])
        if check['stationary']:
            stationary.append(n)
        majority = any(not regular_component_possible(m, n-candidate['center']-m, candidate['center'])
                       for m, b in candidate['rings'])
        if majority:
            require(not check['stationary'], 'general component obstruction agrees with remainder')
            majority_excluded.append(n)
        # Scale all squared radii by a rational factor: stationarity and gain unchanged.
        if n <= 30:
            scaled = ring_stationarity([(m,F(7,3)*b) for m,b in candidate['rings']],
                                       candidate['signs'], candidate['center'])
            require(scaled['scale_invariant_residual_lower'] == check['scale_invariant_residual_lower'],
                    'scale invariant residual bound')
            require(scaled['log_global_gap_lower'] == check['log_global_gap_lower'], 'scale invariant gain')
        rows.append({'n':n, 'rings': candidate['rings'], 'signs':candidate['signs'],
            'center':candidate['center'], 'stationary':check['stationary'],
            'all_phases_excluded_by_population':majority,
            'scale_invariant_residual_lower':str(check['scale_invariant_residual_lower']),
            'log_global_gap_lower':str(check['log_global_gap_lower']),
            'ring_remainders_sha256':[digest_vector(r['remainder']) for r in check['rings']]})
    require(stationary == [2,3,4,5,6,7,9,11], 'documented exact stationary candidates')
    # Eight-point DP uses two rings, unlike the retained stronger center+heptagon.
    bad = 0
    for action in (lambda:ring_stationarity([],()), lambda:ring_stationarity(((3,1),),(1,)),
                   lambda:ring_stationarity(((2,1),(2,1)),(1,1)),
                   lambda:ring_stationarity(((2,1),),(0,)),
                   lambda:regular_component_possible(1,3), lambda:binomial_remainder([1],2,0)):
        try:
            action()
        except ValueError:
            bad += 1
        else:
            raise ArithmeticError('invalid input accepted')
    return {'passed':True, 'candidate_range':[2,100], 'independent_remainders':remainders,
        'independent_direct_force_examples':len(examples), 'scale_invariance_checks':29,
        'invalid_inputs_rejected':bad, 'stationary_selected_candidates':stationary,
        'nonstationary_selected_count':99-len(stationary),
        'entire_continuous_phase_families_excluded_by_population':majority_excluded,
        'scope':'A nonzero remainder excludes exact stationarity of the specified candidate. The count obstruction excludes every radius and phase for the stated population pattern; passing either test proves no optimality.',
        'rows':rows}


if __name__ == '__main__':
    result=run()
    Path(__file__).with_name('stationarity_verification.json').write_text(
        json.dumps(result,ensure_ascii=False,default=str,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},ensure_ascii=False,indent=2))
