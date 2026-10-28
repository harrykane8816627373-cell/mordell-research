"""Exact acceptance checks for population-aware radial empty-gap energies."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

from population_gap import (band_gap_certificate, band_gap_population,
                            free_gap_energy_interval, gap_energy_interval,
                            integer_gap_interval, lower_support_ge,
                            upper_support_le, outward_decimal, log_interval)


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def main(write_report=False):
    root_checks = 0
    for Rnum, rnum, boundary_num in product(range(11), range(1,7), range(16)):
        R, r, t = F(Rnum,5), F(rnum,7), F(boundary_num,5)
        require(upper_support_le(R*R, r*r, t*t) == ((R+r)**2 <= t*t),
                'outward support comparison against rational radii')
        require(lower_support_ge(R*R, r*r, t*t)
                == (max(R-r,F(0))**2 >= t*t),
                'inward support comparison against rational radii')
        root_checks += 2

    gaps = [(F(1,10), F(1,2)), (F(1,5), F(4,5)),
            (F(2,3), F(7,6)), (F(1), F(2))]
    integer_checks = count_tests = 0
    maximum_candidates = 0
    for n in range(2, 31):
        for a, b in gaps:
            for jlo, jhi in ((0,n), (0,n//2), (n//2,n), (n//3,2*n//3)):
                result = integer_gap_interval(n, a, b, jlo, jhi)
                candidates = result['candidate_counts']
                maximum_candidates = max(maximum_candidates, len(candidates))
                full = [(j, gap_energy_interval(a,b,F(j,n))) for j in range(jlo,jhi+1)]
                true_enclosure = (min(interval[0] for j, interval in full),
                                  min(interval[1] for j, interval in full))
                require(result['energy_interval'] == true_enclosure,
                        'convex integer candidates agree with all stated counts')
                for j, interval in full:
                    require(result['energy_interval'][0] <= interval[1],
                            'integer minimum is a lower bound for every count')
                    count_tests += 1
                # Integer restriction cannot improve the free minimum downward.
                free = free_gap_energy_interval(a,b)
                require(result['energy_interval'][0] > free[1],
                        'strict integer population penalty in specified tests')
                integer_checks += 1

    # Direct polynomial integration when the CDF is zero (including a=0).
    zero_mass_checks = 0
    for a, b in [(F(0),F(1,4)), (F(0),F(2))]+gaps:
        exact = (b*b-a*a)/16
        require(gap_energy_interval(a,b,0) == (exact,exact), 'zero-CDF exact integral')
        zero_mass_checks += 1

    bands = [(F(1,25),F(1,20),10), (F(109,100),F(11,10),90)]
    gap = (F(1,5),F(4,5))
    require(band_gap_population(100,bands,*gap) == 10, 'ten interior smeared disks')
    example = band_gap_certificate(100,bands,[gap])
    log2lo, log2hi = log_interval(F(2))
    expected = 75+100*log2lo, 75+100*log2hi
    energy = example['radial_energy_interval']
    require(10000*energy[0] <= expected[1]
            and expected[0] <= 10000*energy[1], 'exact example 75+100 log 2')
    # A one-term atanh lower estimate log2>2/3 suffices for exclusion.
    rational_delta_lower, global_gap_upper = F(425,3), F(1657,12)
    require(example['delta_lower'] > rational_delta_lower > global_gap_upper,
            'entire two-band radial population class excluded')
    # Nonempty centered example: two centered regular polygons.
    require(10*F(1,20)+90*F(197,180) == 99, 'centered two-polygon moment')
    require(F(109,100) < F(197,180) < F(11,10), 'outer example lies in band')
    # Old one-cell polynomial certificate on exactly the same interval.
    a,b,c = *gap,F(1,10)
    previous_cell = (c*c*(b-a)-c*(b*b-a*a)/2+(b**3-a**3)/12)/(2*b)
    require(10000*previous_cell == F(225,2) < global_gap_upper,
            'old single-cell certificate alone does not exclude this example')
    free_example = free_gap_energy_interval(*gap)
    require(10000*free_example[1] < global_gap_upper,
            'discarding the actual population loses this gap exclusion')

    # Disjoint gap energies add; overlapping gaps are rejected.
    pieces = [(F(1,5),F(2,5)), (F(2,5),F(4,5))]
    split = band_gap_certificate(100,bands,pieces)
    require(split['radial_energy_interval'][0] <= energy[1]
            and energy[0] <= split['radial_energy_interval'][1],
            'exact integration additive under gap subdivision')
    rejected = 0
    for callback in (
        lambda: gap_energy_interval(0,1,F(1,10)),
        lambda: integer_gap_interval(10,F(1,5),F(4,5),9,8),
        lambda: band_gap_population(100,bands,F(1,10),F(4,5)),
        lambda: band_gap_certificate(100,bands,[(F(1,5),F(3,5)),(F(2,5),F(4,5))]),
        lambda: band_gap_population(100,[(F(1),F(2),99)],*gap),
    ):
        try:
            callback()
        except ValueError:
            rejected += 1
        else:
            raise ArithmeticError('invalid population/gap input accepted')

    report = {
        'status': 'PASS', 'arithmetic': 'Fraction and outward atanh logarithm intervals',
        'radical_support_comparisons': root_checks,
        'integer_minimization_cases': integer_checks,
        'all_candidate_counts_cross_checked': count_tests,
        'maximum_retained_candidates': maximum_candidates,
        'zero_mass_integral_cases': zero_mass_checks,
        'invalid_inputs_rejected': rejected,
        'example_n': 100,
        'example_bands': [[str(a),str(b),j] for a,b,j in bands],
        'example_gap': [str(x) for x in gap], 'example_cumulative_count': 10,
        'example_exact_delta_formula': '75 + 100 log(2)',
        'example_delta_interval': [str(x) for x in outward_decimal((10000*energy[0],10000*energy[1]))],
        'example_simple_strict_delta_lower': str(rational_delta_lower),
        'existing_global_optimum_gap_upper': str(global_gap_upper),
        'previous_single_cell_delta_lower': str(10000*previous_cell),
        'free_constant_delta_interval': [str(x) for x in outward_decimal((10000*free_example[0],10000*free_example[1]))],
        'nonempty_example': 'centered regular 10-gon at squared radius 1/20 and regular 90-gon at 197/180',
        'scope': 'Conditional class certificates; integer counts on each gap are minimized independently. No full population-sequence or configuration enumeration is asserted.',
    }
    if write_report:
        Path(__file__).with_name('population_gap_verification.json').write_text(
            json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    main(parser.parse_args().write_report)
