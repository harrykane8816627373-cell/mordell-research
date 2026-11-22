"""Exact replay of a ninth-order odd-moment bound for the ten-point outer pairing."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CURVED = ROOT / 'ten-point-curved-kernel'
RADII = ROOT / 'ten-point-refined-radii'
MOMENTS = ROOT / 'ten-point-outer-moments'
HARMONICS = (1, 3, 5, 7, 9)
BUDGET = F(3, 625)
PROJECTION_CAP = F(21, 40)
ALPHA_COS = F(106, 125)
BETA_COS = F(1, 8)
OUTER_MIDPOINT_MASS_CAP = F(117, 4000)
TOTAL_MIDPOINT_NORM_CAP = F(363, 1000)


def read(path):
    return json.loads(path.read_bytes())


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


moments = load(MOMENTS / 'verify_all.py', 'tight_pairing_prior_moments')
need = moments.need
mul = moments.radial.mul
add = moments.add


def at(coeffs, x):
    value = F(0)
    for a in reversed(coeffs):
        value = a + x * value
    return value


def chebyshevs(n):
    polys = [[F(1)], [F(0), F(1)]]
    while len(polys) <= n:
        polys.append(add([F(0)]+[2*a for a in polys[-1]], polys[-2], -1))
    return polys


def check_inputs():
    manifest = read(HERE / 'source_manifest.json')
    for names, base in (('files', HERE), ('dependencies', ROOT)):
        for name, expected in manifest[names].items():
            need(sha256((base / name).read_bytes()).hexdigest() == expected, 'fixed input ' + name)
    radii = read(RADII / 'verification_report.json')
    need(radii['passed'] and radii['applies_to_every_global_maximizer'] and radii['unresolved']==0,
         'accepted radius and midpoint proof')
    bands = tuple(tuple(map(F, radii[k])) for k in ('inner_squared_radius_open_interval', 'outer_squared_radius_open_interval'))
    need(bands == moments.BANDS and (radii['exact_inner_count'], radii['exact_outer_count']) == (2,8),
         'two inner and eight outer points in the accepted bands')
    inner_midpoint = F(radii['midpoint']['inner_midpoint_squared_strict_upper'])
    need(inner_midpoint == F(57121, 4000000), 'accepted inner midpoint bound')
    prior_moments = read(MOMENTS / 'verification_report.json')
    need(prior_moments['passed'] and F(prior_moments['pairing']['angular_separation_cosine_strict_upper']) == ALPHA_COS,
         'accepted separation of all distinct outer directions')
    curved = read(CURVED / 'verification_report.json')
    kernel = read(CURVED / 'kernel.json')
    need(curved['passed'] and curved['unresolved'] == 0 and curved['PSD_blocks']==10
         and (curved['exact_inner_count'], curved['exact_outer_count']) == (2, 8),
         'accepted complete curved-kernel cover and ten positive definite Gram blocks')
    need(curved['kernel_sha256'] == sha256((CURVED / 'kernel.json').read_bytes()).hexdigest(),
         'covered kernel matches the Gram kernel used below')
    budget = moments.upper.hx(kernel['B']) - moments.upper.log.log_lower(moments.upper.witness())
    need(budget + F(curved['c_inner']) + 16*F(curved['c_mixed']) + 28*F(curved['c']) < BUDGET/2,
         'sharp enough total permitted Gram penalty')
    return kernel, bands, inner_midpoint


def projection(kernel, bands, data, cheb):
    weights = tuple(map(F, data['w']))
    need(len(weights)==len(HARMONICS) and weights[0]==1, 'fixed odd harmonic orders and normalization')
    rows = data['rows']
    es, slope = tuple(map(F, data['e'])), F(data['slope'])
    radius_bands = ((F(272,1000), F(415,1000)), (F(919,1000), F(1146,1000)))
    entries = [moments.coefficients(kernel, k, rows[str(k)]) for k in HARMONICS]
    covers = []
    for band, (rlo, rhi) in enumerate(radius_bands):
        ulo, uhi = bands[band]
        need(0 < rlo*rlo < ulo < uhi < rhi*rhi, 'verified radius enclosure')
        diff = [[F(0)] * 10 for _ in range(20)]
        for k, weight, (g, dual) in zip(HARMONICS, weights, entries):
            need(dual > 0 and len(g) == 6, 'positive exact dual norm and fixed radial feature order')
            terms = [(0, -band*weight)] + [(2*j+k, a) for j,a in enumerate(g)]
            for i, a in terms:
                for j, b in enumerate(cheb[k]):
                    diff[i][j] += a*b
        for sign in (-1, 1):
            polynomial = [[sign*x for x in row] for row in diff]
            polynomial[0][0] += es[band] - F(9,10)*slope
            polynomial[2][0] += slope
            covers.append(moments.check_bivariate(polynomial, (rlo,rhi)))
    common_error = 2*es[0] + 8*es[1]
    dual = sum(d for _,d in entries)
    cap = PROJECTION_CAP - common_error
    need(cap > 0 and cap**2 > BUDGET*dual, 'joint projection cap by the squared dual Gram inequality')
    return weights, {
        'global_projection_modulus_strict_upper': str(PROJECTION_CAP),
        'combined_exact_approximation_error': str(common_error),
        'dual_gram_norm_squared_decimal_display': float(dual),
        'two_variable_covers': covers,
        'two_variable_nodes': sum(c['visited'] for c in covers),
    }


def geometry(weights, cheb, outer, inner_midpoint):
    a,b = ALPHA_COS, BETA_COS
    shape, curv = [F(0)], [F(0)]
    for k, weight in zip(HARMONICS, weights):
        shape = add(shape, cheb[k], weight)
        curv = add(curv, cheb[k], k*k*weight)
    need(all(x==0 for x in curv[::2]), 'curvature polynomial odd in the cosine')
    concavity = moments.radial.strictly_above(curv[1::2], (F(0), F(1)), F(0))
    result = read(HERE / 'jacobian_concavity.json')
    sturm = result['output']
    source = sturm['source_polynomial']
    terms = source['polynomial']['terms']
    coefficients = [F(0)]*5
    need(result['operation_id']=='polynomial.root_count.compute'
         and source['domain']=='QQ' and source['variables']==['s']
         and [term['exponents'] for term in terms]==[[4],[3],[2],[1],[0]],
         'Jacobian checked the claimed univariate quartic')
    for term in terms:
        raw = term['coefficient']
        coefficients[term['exponents'][0]] = F(int(raw['num']),int(raw['den']))
    need(coefficients==curv[1::2] and sturm['root_count']==0
         and sturm['lower']=={'num':'0','den':'1'} and sturm['upper']=={'num':'1','den':'1'},
         'Jacobian root count matches the quartic verified independently by Bernstein')
    val = (1-b*b)*(1-a*a)
    scale = 1<<60
    x = isqrt((val*scale**2).numerator // (val*scale**2).denominator)
    lo,hi = F(x,scale), F(x+1,scale)
    need(lo*lo <= val < hi*hi, 'square root enclosure for the first cosine recurrence step')
    cb = [(b,b),(a*b+lo,a*b+hi)]
    for j in range(2,5):
        cb.append((2*a*cb[j-1][0]-cb[j-2][1], 2*a*cb[j-1][1]-cb[j-2][0]))
    need(all(0<c<=d<1 for c,d in cb), 'each five-point arc vertex projects to the forward semicircle')
    values = [at(shape, pair[0]) for pair in cb]
    vertices = [sum(values[:j])+sum(values[:5-j]) for j in range(6)]
    apex = at(shape,1)
    other_upper = apex+2*at(shape,a)
    need(other_upper > 2*apex, 'three-point projection upper bound also holds with a backward point')
    net = min(vertices)-other_upper
    need(net > PROJECTION_CAP, 'five points in a closed arc cause a forbidden odd projection')
    kappa = 1-2*b*b
    for u in outer:
        for v in outer:
            w = u+v-4*OUTER_MIDPOINT_MASS_CAP
            need(w>0 and w*w < 4*kappa*kappa*u*v, 'all endpoints certify the outer midpoint bound')
    need(4*OUTER_MIDPOINT_MASS_CAP+inner_midpoint < TOTAL_MIDPOINT_NORM_CAP**2 ,
         'four outer midpoint bounds and the inner midpoint bound imply the joint norm cap')
    return {
        'concavity_polynomial_ascending_coefficients': list(map(str, curv[1::2])),
        'concavity_replay': concavity,
        'jacobian_closed_unit_interval_root_count': sturm['root_count'],
        'closed_arc_half_angle_cosine': str(b),
        'closed_arc_occupancy_upper': 4,
        'five_point_arc_vertices_checked': len(vertices),
        'arc_odd_projection_lower_decimal_display': float(net),
        'cyclic_rank_pairing_cosine_strict_upper': str(-kappa),
        'paired_squared_distance_strict_lower': str(2*outer[0]*(1+kappa)),
        'paired_outer_midpoint_squared_strict_upper': str(OUTER_MIDPOINT_MASS_CAP),
        'joint_five_pair_midpoint_norm_strict_upper': str(TOTAL_MIDPOINT_NORM_CAP),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    kernel, bands, inner_midpoint = check_inputs()
    cheb = chebyshevs(max(HARMONICS))
    w, result = projection(kernel, bands, read(HERE / 'projection.json'), cheb)
    report = {'passed': True, 'applies_to_every_global_maximizer': True, 'normalization_mass': 9,
              'projection': result, 'pairing': geometry(w, cheb, bands[1], inner_midpoint)}
    saved = HERE / 'verification_report.json'
    text = json.dumps(report, indent=2) + '\n'
    if args.record_report:
        saved.write_text(text, encoding='utf-8', newline='\n')
    else:
        need(report == read(saved), 'default replay matches the saved report')
    print(text, end='')


if __name__ == '__main__':
    main()
