"""Replay exact population/radius search against the COMPLETE prior envelope."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import sys

from population_search import (LEGACY, optimized_radii, proxy_value,
                               population_dp, ring_value, population_search,
                               continuous_phase_upper)

PRIOR = Path(__file__).resolve().parents[1] / 'joint-bounds'
if str(PRIOR) not in sys.path:
    sys.path.insert(0, str(PRIOR))
from verify_phase_search import prior_ascent_values
from phase_search import contracted_polar
from phase_graph import optimized_polar
from rational_ascent import value as point_value
from compare_envelopes import verify_saved
from baseline_finite_lower_bounds import gaussian_lower, central_polygon_lower
from triangular_lattice import (lattice_certificate, gaussian_norm,
                                triangular_norm, improved_centered_polygon)
from polar_refinement import refined_polar, twelve_two_ring_maximum

ROOT = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def digest(q):
    def h(a):
        return hashlib.sha256(a.to_bytes(max(1, (a.bit_length()+7)//8), 'big')).hexdigest()
    return {'numerator_sha256': h(q.numerator), 'denominator_sha256': h(q.denominator)}


def even_compositions(n):
    if n == 0:
        yield ()
    for m in range(2, n+1, 2):
        for tail in even_compositions(n-m):
            yield (m,)+tail


def audit(max_n=100):
    if type(max_n) is not int or not 2 <= max_n <= 100:
        raise ValueError('comparison requires integer 2<=max_n<=100')
    # A ring plus the origin checks the added origin factors directly.
    square = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    require(ring_value(5, ((4, 1),), (1,), 1) == point_value(square+[(0, 0)]),
            'independent origin and four direct pair factors')
    double_square = square+[(2*x, 2*y) for x, y in square]
    require(ring_value(8, ((4, 1), (4, 4)), (1, 1)) == point_value(double_square),
            'independent double-square pair product')
    require(ring_value(9, ((4, 1), (4, 4)), (1, 1), 1) == point_value(double_square+[(0, 0)]),
            'independent two rings plus origin pair product')
    # Exhaust all compositions; do not reproduce the DP recurrence.
    composition_checks = 0
    for n in range(2, min(max_n, 22)+1):
        h = n % 2
        candidates = []
        for pops in even_compositions(n-h):
            rings = optimized_radii(pops, h)
            require(sum(m*b for m, b in rings) == n*(n-1), 'moment identity for arbitrary composition')
            candidates.append(proxy_value(n, rings))
            composition_checks += 1
        require(max(candidates) == population_dp(n)['surrogate'], 'independent complete composition maximum')
    try:
        population_dp(100, max_transitions=1)
    except ValueError as e:
        require('above budget' in str(e), 'DP budget refusal reason')
    else:
        raise ArithmeticError('DP incomplete enumeration refusal missing')

    saved = verify_saved()
    moments = json.loads((LEGACY/'moment_values.json').read_text(encoding='utf-8-sig'))
    rayleigh = {r['n']: F(r['rayleigh_lower']) for r in moments['rows']}
    phase_report = json.loads((LEGACY/'phase_graph_comparison.json').read_text(encoding='utf-8'))
    phase_rows = {r['n']: r for r in phase_report['rows']}
    contraction_report = json.loads((PRIOR/'phase_search_verification.json').read_text(encoding='utf-8'))
    contraction_rows = {r['n']: r for r in contraction_report['rows']}
    ascents, ascent_count = prior_ascent_values()
    rows = []
    explicit_wins = []
    existence_wins = []
    surrogate_stronger = []
    continuous_phase_exact = []
    sign_checks = 0
    for n in range(2, max_n+1):
        old_phase = optimized_polar(n)['value']
        prior = max(F(n**n), central_polygon_lower(n), gaussian_lower(n),
                    lattice_certificate(n, gaussian_norm)['value'],
                    lattice_certificate(n, triangular_norm)['value'],
                    refined_polar(n)['value'], saved.get(n, F(0)), rayleigh.get(n, F(0)),
                    twelve_two_ring_maximum() if n == 12 else F(0))
        perturbation = improved_centered_polygon(n)
        if perturbation is not None:
            prior = max(prior, perturbation)
        require(digest(prior) == phase_rows[n]['previous_envelope'], 'complete earlier envelope digest')
        require(digest(old_phase) == phase_rows[n]['value'], 'earlier phase digest')
        prior = max(prior, old_phase, ascents.get(n, F(0)))
        if n >= 13:
            require(digest(prior) == contraction_rows[n]['previous_complete'], 'complete pre-contraction envelope')
            contraction = contracted_polar(n)['value']
            require(digest(contraction) == contraction_rows[n]['value'], 'earlier contraction witness')
            prior = max(prior, contraction)
        result = population_search(n)
        require(sum(m*b for m, b in result['rings']) == n*(n-1), 'actual candidate exact moment')
        if n <= 30:
            independent = max(ring_value(n, result['rings'], signs, result['center'])
                              for signs in product((-1, 1), repeat=len(result['rings'])))
            sign_checks += 1 << len(result['rings'])
            require(independent == result['value'], 'independent full actual +/- phase maximum')
        require(result['continuous_phase_upper'] >= result['value'], 'all-continuous-phase upper covers explicit candidate')
        require(result['continuous_phase_upper'] >= result['surrogate'], 'all-continuous-phase upper covers angular-average existence bound')
        require(result['continuous_phase_upper'] <= (1+F(1, 1000000))*result['value'],
                'selected fixed-radius continuous-phase gap at most one ppm')
        if result['continuous_phase_upper'] == result['value']:
            continuous_phase_exact.append(n)
        if result['value'] > prior:
            explicit_wins.append(n)
        if result['lower_bound'] > prior:
            existence_wins.append(n)
        if result['surrogate'] > result['value']:
            surrogate_stronger.append(n)
        gain = (result['value']/prior-1)*10000
        angular_gap_ppm = (result['continuous_phase_upper']/result['value']-1)*1000000
        row = {'n': n, 'center': result['center'],
               'populations': result['populations'], 'rings': result['rings'],
               'signs': result['signs'], 'population_dp_transitions': result['transitions'],
               'phase_states': result['phase_states'],
               'actual': digest(result['value']), 'surrogate': digest(result['surrogate']),
               'continuous_phase_upper_at_selected_radii': digest(result['continuous_phase_upper']),
               'continuous_phase_relative_gap_ppm_ceiling': -((-angular_gap_ppm.numerator)//angular_gap_ppm.denominator),
               'previous_complete': digest(prior),
               'explicit_candidate_improves_complete_previous': result['value'] > prior,
               'existence_lower_improves_complete_previous': result['lower_bound'] > prior,
               'explicit_gain_basis_points_floor': gain.numerator//gain.denominator}
        if n <= 30:
            row.update(actual_exact=str(result['value']), surrogate_exact=str(result['surrogate']))
        rows.append(row)
    return {'passed': True, 'n_min': 2, 'n_max': max_n,
            'independent_direct_coordinate_examples': 3,
            'independent_ordered_even_compositions_checked': composition_checks,
            'independent_small_family_sign_values_checked': sign_checks,
            'prior_ascent_witnesses_recomputed': ascent_count,
            'strict_explicit_complete_previous_envelope_improvements': explicit_wins,
            'strict_existence_complete_previous_envelope_improvements': existence_wins,
            'surrogate_above_actual_sign_candidate_at': surrogate_stronger,
            'continuous_phase_exact_at_selected_radii': continuous_phase_exact,
            'all_selected_continuous_phase_relative_gaps_at_most': '1/1000000',
            'scope': 'DP is complete only for the angular-log-average surrogate over even populations with parity-fixed center. The actual candidate optimizes all +/- signs at the DP-selected radii. Prior envelope includes all published one-step ascent and contraction gains. No general optimal configuration, complete actual population optimum or new linear asymptotic coefficient is claimed.',
            'rows': rows}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--max-n', type=int, default=100)
    p.add_argument('--output', type=Path, default=ROOT/'population_search_verification.json')
    args = p.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    result = audit(args.max_n)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k != 'rows'}, ensure_ascii=False, indent=2))
