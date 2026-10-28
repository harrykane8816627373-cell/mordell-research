"""Replay the uniform radial theorem and compare complete prior bounds.

No floating-point decisions, network, or optional dependencies are used.
The all-n statement is proved in RADIAL_CONTRACTION.md; finite replay
checks the eleven exceptional n and the published comparison range.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import sys

from phase_search import (LEGACY, normalized_ring_value, worst_phase_gain,
                          analytic_log_gain_floor, contracted_polar)
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
        return hashlib.sha256(a.to_bytes(max(1, (a.bit_length() + 7) // 8), 'big')).hexdigest()
    return {'numerator_sha256': h(q.numerator), 'denominator_sha256': h(q.denominator)}


def prior_ascent_values():
    report = json.loads((LEGACY / 'ascent_envelope_comparison.json').read_text(encoding='utf-8'))
    result = {}
    count = 0
    for row in report['rows']:
        reconstructed = [point_value(c['integer_witness'], c['d']) for c in row['certificates']]
        for q, c in zip(reconstructed, row['certificates']):
            require(digest(q) == c['after'], 'prior exact ascent witness digest')
            count += 1
        result[row['n']] = max(reconstructed)
        require(digest(result[row['n']]) == row['best_ascent'], 'prior best ascent digest')
    require(sorted(result) == list(range(2, 31)), 'prior ascent comparison range')
    return result, count


def audit(max_n=100):
    if type(max_n) is not int or not 13 <= max_n <= 100:
        raise ValueError('comparison requires 13 <= max_n <= 100')
    # Direct rational-coordinate checks are independent of ring resultants.
    axis_square = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    double_square = axis_square + [(2*x, 2*y) for x, y in axis_square]
    require(normalized_ring_value(8, ((4, 1), (4, 4)), (1, 1)) == point_value(double_square),
            'independent eight rational point pairs')
    singleton = axis_square + [(3, 0)]
    require(normalized_ring_value(5, ((4, 1), (1, 9)), (1, 1)) == point_value(singleton),
            'independent singleton centering and normalization')
    # All-n proof constants. These checks are not a numerical extrapolation.
    require(F(69, 32) - F(225, 112) - F(1, 336) - F(12, 343) == F(3601, 32928),
            'exact uniform logarithmic floor')
    require(F(3601, 32928) > F(1, 10), 'uniform logarithmic floor exceeds one tenth')
    require(6 * F(6, 37)**3 < F(12, 343), 'small remainder interaction bound')
    finite_proof = []
    for n in range(13, 24):
        g = worst_phase_gain(n)
        require(g > F(21, 19), 'exceptional n exact worst-phase lower bound')
        finite_proof.append({'n': n, 'worst_phase_gain': str(g), 'strictly_above_21_over_19': True})
    require(analytic_log_gain_floor(24) == F(3601, 32928), 'all-n floor base case')
    require(analytic_log_gain_floor(10**6) > analytic_log_gain_floor(24), 'floor monotonicity exercise')

    saved = verify_saved()
    moments = json.loads((LEGACY / 'moment_values.json').read_text(encoding='utf-8-sig'))
    rayleigh = {r['n']: F(r['rayleigh_lower']) for r in moments['rows']}
    published_phase = json.loads((LEGACY / 'phase_graph_comparison.json').read_text(encoding='utf-8'))
    previous_rows = {r['n']: r for r in published_phase['rows']}
    ascents, ascent_count = prior_ascent_values()
    improved = []
    rows = []
    sign_checks = 0
    for n in range(13, max_n + 1):
        old_phase = optimized_polar(n)['value']
        prior = max(F(n**n), central_polygon_lower(n), gaussian_lower(n),
                    lattice_certificate(n, gaussian_norm)['value'],
                    lattice_certificate(n, triangular_norm)['value'],
                    refined_polar(n)['value'], saved.get(n, F(0)), rayleigh.get(n, F(0)),
                    twelve_two_ring_maximum() if n == 12 else F(0))
        perturbation = improved_centered_polygon(n)
        if perturbation is not None:
            prior = max(prior, perturbation)
        require(digest(prior) == previous_rows[n]['previous_envelope'], 'complete published envelope digest')
        require(digest(old_phase) == previous_rows[n]['value'], 'prior phase value digest')
        prior = max(prior, old_phase, ascents.get(n, F(0)))
        candidate = contracted_polar(n)
        gain = candidate['value'] / prior
        require(candidate['value'] > F(11, 10) * old_phase, 'all-n theorem exercised in range')
        if n <= 30:
            # Re-enumerate the entire new finite sign family independently
            # of the block graph solver used to produce the candidate.
            best = max(normalized_ring_value(n, candidate['rings'], s)
                       for s in product((-1, 1), repeat=len(candidate['rings'])))
            sign_checks += 1 << len(candidate['rings'])
            require(best == candidate['value'], 'independent full new sign enumeration')
        is_improved = gain > 1
        if is_improved:
            improved.append(n)
        gain_basis_points = (gain - 1) * 10000
        row = {'n': n,
               'rings': [[m, str(b)] for m, b in candidate['rings']],
               'signs': candidate['signs'], 'phase_states': candidate['phase_states'],
               'value': digest(candidate['value']), 'previous_complete': digest(prior),
               'strict_complete_envelope_improvement': is_improved,
               'gain_basis_points_floor_over_complete_previous': gain_basis_points.numerator // gain_basis_points.denominator,
               'gain_over_previous_phase': digest(candidate['value'] / old_phase)}
        if n <= 30:
            row['value_exact'] = str(candidate['value'])
        rows.append(row)
    return {'passed': True, 'n_min': 13, 'n_max': max_n,
            'all_n_theorem': 'For every n>=13 and every continuous old-ring rotation, log gain>1/10 and product gain>11/10.',
            'all_n_log_gain_floor_from_n24': '3601/32928',
            'exceptional_n_exact_proof': finite_proof,
            'independent_direct_coordinate_examples': 2,
            'independent_small_family_sign_values_checked': sign_checks,
            'prior_ascent_witnesses_recomputed': ascent_count,
            'strict_complete_previous_envelope_improvements': improved,
            'scope': 'New fixed-radius +/- family values compared with the entire published envelope, retaining all n<=30 published one-step ascents. No unrestricted maximum or new linear asymptotic coefficient is claimed.',
            'rows': rows}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--max-n', type=int, default=100)
    p.add_argument('--output', type=Path, default=ROOT / 'phase_search_verification.json')
    args = p.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    result = audit(args.max_n)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'exceptional_n_exact_proof')},
                     ensure_ascii=False, indent=2))
