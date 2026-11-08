"""Conservative necessary-condition exclusions; never existence certificates."""
from collections import Counter, defaultdict
from itertools import product


def require(ok, message):
    if not ok:
        raise ValueError(message)


def subset_sums(populations):
    values=[0]
    for a in populations:
        values += [e+a for e in values]
    return values


def isolated_force_exponents(populations, largest):
    """A single raw subset term has a nonzero monomial coefficient at any phase."""
    groups=defaultdict(list)
    for e in subset_sums(populations):
        groups[e % largest].append(e)
    return sorted(rows[0] for rows in groups.values() if len(rows)==1)


def distinct_subset_sum_obstruction(populations, largest):
    require(type(largest) is int and largest>=2,'largest population >=2 required')
    require(populations and all(type(a) is int and 2<=a<largest for a in populations),
            'all exterior populations must be smaller integers >=2')
    sums=subset_sums(populations)
    applies=sum(populations)<=2*largest and len(sums)==len(set(sums))
    exponents=isolated_force_exponents(populations,largest)
    return {'applies':applies,'exterior_total':sum(populations),
            'all_subset_sums_distinct':len(sums)==len(set(sums)),
            'isolated_force_exponents':exponents,
            'two_incompatible_force_constants':[-e for e in exponents[:2]] if len(exponents)>=2 else []}


def classify_four(populations):
    require(len(populations)==4 and all(type(a) is int and a>=2 for a in populations),
            'four populations >=2 required')
    a,b,c,m=sorted(populations)
    if a==m:
        return 'not_excluded_equal_populations'
    if c==m:
        return 'at_most_two_smaller_rings'
    if a==c and m!=2*a:
        return 'three_equal_outlier_not_doubled'
    exterior=(a,b,c);M=a+b+c
    if M<m:
        return 'degree_obstruction'
    if M==m:
        return 'balanced_four_ring_obstruction'
    if m % 2==0 and exterior.count(m//2)==2:
        return 'two_half_population_rings'
    if distinct_subset_sum_obstruction(exterior,m)['applies']:
        return 'distinct_subset_sums_low_degree'
    if len(isolated_force_exponents(exterior,m))>=2:
        return 'two_isolated_force_coefficients'
    return 'not_excluded_by_these_theorems'


def certificate():
    census=Counter();total=0
    for a in range(2,26):
        for b in range(a,(100-a)//3+1):
            for c in range(b,(100-a-b)//2+1):
                for m in range(c,100-a-b-c+1):
                    census[classify_four((a,b,c,m))]+=1;total+=1
    examples=[(2,2,2,2),(2,2,2,3),(2,2,2,4),(4,4,4,6),(2,2,4,6),(2,4,4,6),(2,4,6,8),
              (2,3,5,5),(2,2,5,5),(2,5,5,5),(2,3,4,9),(2,2,2,6),
              (2,3,3,6),(2,2,3,4),(2,4,8,10),(4,7,9,10)]
    sample=[]
    for row in examples:
        ordered=sorted(row);m=ordered[-1];outside=ordered[:-1]
        sample.append({'populations':list(row),'case':classify_four(row),
                       'isolated_force_exponents':isolated_force_exponents(outside,m)})
    return {'date':'2026-11-08',
            'hypotheses':'unit strengths, distinct positive regular-ring radii, arbitrary continuous phases, optional one center',
            'largest_population_defect_rule':'if populations differ, at least three rings have population strictly below the maximum',
            'four_ring_unequal_pattern_rule':'maximum population is unique and strictly smaller than the sum of the other three',
            'low_degree_resonance_rule':'if the other three sum to at most twice the unique maximum, two of them coincide or the largest of them equals the sum of the other two',
            'extra_excluded_family':'(a,s,s,2s), 2<=a<2s, a!=s',
            'three_equal_plus_outlier_rule':'(s,s,s,m), m!=s, requires m=2s',
            'census':{'maximum_noncentral_total':100,'profiles_checked':total,'counts':dict(sorted(census.items())),
                      'finite_census_replaces_general_proof':False},
            'examples':sample,
            'balanced_equal_smaller_branch':{
                'profile':'(s,s,s,3s), s>=2',
                'first_branch':'sum A_i=0 implies sum |A_i|^(2/s) A_i=0 by the full Ward identity; forces equal radii',
                'second_branch_mean_upper':'11/15',
                'second_branch_amgm_lower':'2^(-1/3)',
                'rational_cube_comparison':{'upper_cube_numerator':1331,'upper_cube_denominator':3375,
                                            'amgm_lower_cube':'1/2'}},
            'three_equal_half_degree_boundary':{
                'profile':'(2q,2q,2q,3q), q>=1',
                'second_branch_mean_upper':'3/4',
                'second_branch_amgm_lower':'2^(-1/3)',
                'rational_cube_comparison':{'upper_cube':'27/64','amgm_lower_cube':'1/2'}},
            'general_n_global_optima_determined':False,
            'all_four_ring_populations_equal_claimed':False,
            'remaining_patterns_certified_to_exist':False}


if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2))
