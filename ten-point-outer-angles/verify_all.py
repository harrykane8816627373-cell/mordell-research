"""Exact rational checks for ten-point outer angle geometry and curvature."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INNER = (F(27,500), F(217,1000))
OUTER = (F(37,50), F(139,100))
MIN_DISTANCE_SQUARED = F(2,5)
INNER_COSINE = -F(2,3)
MIDPOINT_SQUARED = F(1,40)
NARROW_INNER = (F(37,500),F(43,250))
NARROW_OUTER = (F(169,200),F(164,125))
NARROW_MIDPOINT = F(239,2000)


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def check_inputs():
    raw = (HERE/'source_manifest.json').read_bytes()
    manifest = json.loads(raw)
    need(set(manifest['files']) == {'.gitattributes', 'README.md', 'verify_all.py'},
         'complete list of local inputs')
    dependencies = manifest['dependencies']
    need(set(dependencies) == {
        'maximizer-domain/verification_report.json',
        'stationary-kernel/README.md',
        'ten-point-banded-kernel-refined/README.md',
        'ten-point-banded-kernel-refined/verification_report.json',
        'ten-point-counted-radii/verification_report.json',
        'ten-point-refined-radii/README.md',
        'ten-point-refined-radii/verification_report.json',
        'ten-point-tail-exclusion/verification_report.json',
    }, 'complete list of accepted prerequisites')
    for directory, section in ((HERE, 'files'), (ROOT, 'dependencies')):
        for path, sha in manifest[section].items():
            need(hashlib.sha256((directory/path).read_bytes()).hexdigest() == sha,
                 'fixed input '+path)

    refined = json.loads((ROOT/'ten-point-banded-kernel-refined/verification_report.json').read_bytes())
    inner = refined['inner_pair']
    need(refined['passed'] and refined['unresolved'] == 0
         and refined['exact_inner_count'] == 2 and refined['exact_outer_count'] == 8
         and inner['unresolved'] == 0 and inner['applies_to_every_global_maximizer']
         and tuple(map(F, inner['inner_squared_radius_open_interval'])) == INNER
         and F(inner['inner_midpoint_squared_strict_upper']) == MIDPOINT_SQUARED
         and F(inner['angle_cosine_strict_upper']) == INNER_COSINE,
         'refined inner count, radius, midpoint and angle')
    tail = json.loads((ROOT/'ten-point-tail-exclusion/verification_report.json').read_bytes())
    counted = json.loads((ROOT/'ten-point-counted-radii/verification_report.json').read_bytes())
    old_outer = tuple(map(F, tail['outer_squared_radius_open_interval']))
    expected_exclusions = {(old_outer[0], OUTER[0]), (OUTER[1], old_outer[1])}
    exclusions = {tuple(map(F, row['excluded_closed_interval'])) for row in counted['row_exclusions']
                  if row['unresolved'] == 0}
    need(tail['passed'] and tail['all_global_maximizers'] and tail['normalization_mass'] == 9
         and tail['outer_count'] == 8
         and counted['passed'] and counted['applies_to_every_global_maximizer'] and counted['mass'] == 9
         and expected_exclusions <= exclusions, 'certified outer radii and tail removal')
    domain = json.loads((ROOT/'maximizer-domain/verification_report.json').read_bytes())
    need(domain['passed'] and F(domain['pair_squared_distance_strict_lower']) == MIN_DISTANCE_SQUARED,
         'certified minimum squared distance')
    narrow = json.loads((ROOT/'ten-point-refined-radii/verification_report.json').read_bytes())
    need(narrow['passed'] and narrow['applies_to_every_global_maximizer'] and narrow['unresolved'] == 0
         and narrow['normalization_mass'] == 9
         and narrow['exact_inner_count'] == 2 and narrow['exact_outer_count'] == 8
         and tuple(map(F,narrow['inner_squared_radius_open_interval'])) == NARROW_INNER
         and tuple(map(F,narrow['outer_squared_radius_open_interval'])) == NARROW_OUTER
         and narrow['midpoint']['passed']
         and F(narrow['midpoint']['inner_midpoint_norm_strict_upper']) == NARROW_MIDPOINT,
         'refined radial bands and inner midpoint')
    return hashlib.sha256(raw).hexdigest()


def check_gap_cosine(corners, gap_cosine):
    for u,v in corners:
        need(u+v > MIN_DISTANCE_SQUARED
             and (u+v-MIN_DISTANCE_SQUARED)**2 < 4*u*v*gap_cosine**2,
             'acute gap cosine bound at all outer band corners')


def six_ray_sine_factor(gap_cosine):
    # Every sine below is a multiple of sin(alpha), where cos(alpha)=gap_cosine.
    sine = [F(0),F(1)]
    for k in range(2,6):
        sine.append(2*gap_cosine*sine[k-1]-sine[k-2])
    vertices = [sum(sine[:k])+sum(sine[:6-k]) for k in range(7)]
    vertex_minimum = 2*(sine[1]+sine[2])
    need(min(vertices) == vertex_minimum, 'all seven angular gap vertices')
    return vertex_minimum


def check_ward_weights(band, center, lower, upper):
    critical = center/3
    need(0 < band[0] < critical < band[1] < center, 'single interior maximum of squared weight')
    need(all(u*(center-u)**2 > lower**2 for u in band)
         and critical*(center-critical)**2 < upper**2, 'outer Ward weight range')


def semicircle_bounds(corners):
    gap_cosine = F(6,7)
    check_gap_cosine(corners,gap_cosine)
    vertex_minimum = six_ray_sine_factor(gap_cosine)
    need(vertex_minimum == F(38,7), 'original six-ray sine factor')
    sine_lower = F(279,100)
    need(vertex_minimum**2*(1-gap_cosine**2) > sine_lower**2,
         'strict sine projection lower bound')

    weight_center = F(13,4)
    weight_lower, weight_upper = F(43,20), F(23,10)
    check_ward_weights(OUTER,weight_center,weight_lower,weight_upper)
    midpoint_bound, inner_radius_bound = F(4,25), F(7,15)
    need(MIDPOINT_SQUARED < midpoint_bound**2 and INNER[1] < inner_radius_bound**2,
         'simple rational majorants of root lengths')
    inner_upper = 2*weight_center*midpoint_bound+2*inner_radius_bound*INNER[1]
    outer_lower = weight_lower*sine_lower-2*weight_upper
    need(outer_lower > inner_upper > 0, 'six outer rays contradict the first Ward moment')
    return gap_cosine, vertex_minimum, outer_lower, inner_upper


def check_outer_curvature(corners, lower):
    ratio_upper = 2*(1-lower)/lower
    for u,v in corners:
        need((u+v)*MIN_DISTANCE_SQUARED-(u-v)**2 > lower*MIN_DISTANCE_SQUARED**2,
             'outer negative curvature at minimum d on four corners')
        need((u+v)**2 < ratio_upper**2*u*v, 'outer radius imbalance at four corners')


def curvature_bounds(corners):
    outer = F(12,25)
    check_outer_curvature(corners,outer)
    max_ratio = F(7,15)/F(63,100)
    inner_near = F(13,5)
    need(INNER[1] < F(7,15)**2 and F(63,100)**2 < MIN_DISTANCE_SQUARED
         and 2*(max_ratio**2+max_ratio) < inner_near
         and 2*(OUTER[0]-INNER[1]) > MIN_DISTANCE_SQUARED,
         'near inner curvature and boundary monotonicity')

    far_cosine, ratio_low, ratio_high = F(5,12), F(19,100), F(11,20)
    need((1+INNER_COSINE)/2 < far_cosine**2, 'at least one far inner direction')
    need(ratio_low**2 < INNER[0]/OUTER[1] and INNER[1]/OUTER[0] < ratio_high**2 < 1,
         'inner over outer radius bounds')
    need(F(8,17) > far_cosine, 'nonpositive far curvature above radius ratio 1/4')
    inner_far = F(1,16)
    denominator_lower = (1-far_cosine/2)**2
    numerator_upper = ((1+F(1,4)**2)*far_cosine-2*ratio_low)/2
    need(0 < numerator_upper < denominator_lower*inner_far, 'positive far curvature upper bound')
    stiffness = 8*outer-inner_near-inner_far
    need(stiffness == F(471,400) and stiffness > 0, 'mean-zero angular Hessian gap')
    return outer, inner_near, inner_far, stiffness


def narrow_bounds():
    corners = tuple(product(NARROW_OUTER,repeat=2))
    gap_cosine, outer_curvature = F(17,20), F(49,100)
    check_gap_cosine(corners,gap_cosine)
    check_outer_curvature(corners,outer_curvature)
    radius_inner_upper = F(83,200)
    radius_outer_lower, radius_outer_upper = F(9,10), F(23,20)
    radius_gap, inner_curvature = F(63,100), F(219,100)
    near_ratio = radius_inner_upper/radius_gap
    need(NARROW_INNER[1] < radius_inner_upper**2 and radius_gap**2 < MIN_DISTANCE_SQUARED
         and 2*(NARROW_OUTER[0]-NARROW_INNER[1]) > MIN_DISTANCE_SQUARED
         and 2*(near_ratio+near_ratio**2) < inner_curvature, 'narrow near inner curvature')
    ratio_floor = F(23,100)
    need(ratio_floor**2 < NARROW_INNER[0]/NARROW_OUTER[1]
         and NARROW_INNER[1] < NARROW_OUTER[0]
         and 2*ratio_floor/(1+ratio_floor**2) > F(5,12), 'narrow far inner curvature is nonpositive')
    stiffness = 8*outer_curvature-inner_curvature
    need(stiffness == F(173,100), 'narrow mean-zero angular Hessian gap')

    # Compress an extended arc to a semicircle while preserving a lower angle gap.
    alpha_upper, shift, compressed_cosine = F(14,25), F(1,20), F(43,50)
    need(1-alpha_upper**2/2+alpha_upper**4/24 < gap_cosine, 'cosine upper series')
    need((compressed_cosine-gap_cosine)**2 >
         (2*shift*alpha_upper/3)**2*(1-gap_cosine**2), 'angle gap after linear compression')
    sine_lower, weight_center = F(111,40), F(13,4)
    sine_factor = six_ray_sine_factor(compressed_cosine)
    need(sine_factor**2*(1-compressed_cosine**2)>sine_lower**2, 'compressed six-ray sine lower')
    weight_lower, weight_upper = F(221,100), F(282,125)
    check_ward_weights(NARROW_OUTER,weight_center,weight_lower,weight_upper)
    outer_projection = weight_lower*sine_lower-(2+6*shift)*weight_upper
    inner_projection = 2*weight_center*NARROW_MIDPOINT+2*NARROW_INNER[1]*radius_inner_upper
    need(outer_projection-inner_projection == F(611,25000), 'extended arc Ward contradiction')

    # Bounds on the parameter perturbation of the angle gradient.
    cosines = [F(1),gap_cosine]
    for k in range(2,5):
        cosines.append(2*gap_cosine*cosines[-1]-cosines[-2])
    need(radius_outer_lower**2 < NARROW_OUTER[0]
         and NARROW_OUTER[1] < radius_outer_upper**2, 'outer radius conversion')
    lo, hi = NARROW_OUTER
    derivative = (hi-lo)/(4*lo**2)
    pair_coefficients = (F(4),F(3,5),F(1,4),F(1,8))
    need(all(2*derivative**2 < c**2*(1-s)**3
             for c,s in zip(pair_coefficients,cosines[1:])), 'cyclic angular gradient Lipschitz bounds')
    outer_force = 2*(2*sum(pair_coefficients[:3])+pair_coefficients[3])
    inner_force_radii = 4*radius_inner_upper/(MIN_DISTANCE_SQUARED*2*radius_outer_lower)
    dimension_factor = F(20,7)
    need(dimension_factor**2 > 8, 'strict square root of eight majorant')
    inner_force_midpoint = 8*radius_outer_upper*dimension_factor/MIN_DISTANCE_SQUARED
    return {
        'inner_squared_radius_open_interval': list(map(str,NARROW_INNER)),
        'outer_squared_radius_open_interval': list(map(str,NARROW_OUTER)),
        'outer_pair_angular_separation_cosine_strict_upper': str(gap_cosine),
        'fixed_mean_outer_angular_negative_hessian_lower': str(stiffness),
        'closed_arc_length_pi_plus': str(2*shift),
        'outer_count_per_closed_arc_upper': 5,
        'extended_arc_ward_contradiction_gap': str(outer_projection-inner_projection),
        'opposite_index_gradient_radius_coefficient': str(outer_force+inner_force_radii),
        'opposite_index_gradient_inner_midpoint_coefficient': str(inner_force_midpoint),
    }


def verify():
    sha = check_inputs()
    corners = tuple(product(OUTER, repeat=2))
    gap_cosine, sine_factor, outer_sum, inner_sum = semicircle_bounds(corners)
    outer_curvature, near, far, stiffness = curvature_bounds(corners)
    return {
        'passed': True,
        'source_manifest_sha256': sha,
        'normalization_mass': 9,
        'outer_count': 8,
        'outer_pair_angular_separation_cosine_strict_upper': str(gap_cosine),
        'outer_count_per_closed_semicircle_upper': 5,
        'outer_count_per_open_semicircle_lower': 3,
        'six_ray_gap_vertices_checked': 7,
        'six_ray_minimum_sine_sum_factor_times_sine_alpha': str(sine_factor),
        'six_ray_outer_directional_ward_sum_strict_lower': str(outer_sum),
        'inner_directional_ward_sum_modulus_strict_upper': str(inner_sum),
        'outer_pair_negative_curvature_lower': str(outer_curvature),
        'inner_outer_positive_curvature_upper': [str(near), str(far)],
        'fixed_mean_outer_angular_negative_hessian_lower': str(stiffness),
        'at_most_one_angular_stationary_vector_per_fixed_data_convex_branch': True,
        'refined_bands': narrow_bounds(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    out = (json.dumps(verify(),indent=2)+'\n').encode('utf-8')
    if args.record_report:
        (HERE/'verification_report.json').write_bytes(out)
    sys.stdout.buffer.write(out)


if __name__ == '__main__':
    main()
