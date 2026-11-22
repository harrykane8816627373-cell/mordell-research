"""Check a signed quartic estimate for eliminating ten-point pair midpoints."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
AMPLITUDE = F(7, 250)
MIDPOINT = F(149, 500)
INNER_CAP = F(239, 2000)
OUTER_CAP = F(43, 250)
LIFT_BOUND = F(61, 10)


def read(path):
    return json.loads(path.read_bytes())


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def load(path):
    spec = importlib.util.spec_from_file_location('quartic_midpoint_basin', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_inputs(repo):
    manifest = read(HERE / 'source_manifest.json')
    for base, key in ((HERE, 'files'), (repo, 'dependencies')):
        for name, expected in manifest[key].items():
            require(sha256((base / name).read_bytes()).hexdigest() == expected,
                    'fixed input ' + name)
    pairing = read(repo / 'ten-point-tight-pairing/verification_report.json')
    joint = read(repo / 'ten-point-midpoint-budget/verification_report.json')
    radii = read(repo / 'ten-point-refined-radii/verification_report.json')
    for report in (pairing, joint, radii):
        require(report['passed'] and report['applies_to_every_global_maximizer']
                and report['normalization_mass'] == 9, 'certified global pairing premise')
    require(joint['unresolved'] == 0 and F(joint['strict_joint_midpoint_norm_upper']) == MIDPOINT,
            'certified joint midpoint cap')
    require(F(pairing['pairing']['paired_outer_midpoint_squared_strict_upper']) < OUTER_CAP**2
            and F(radii['midpoint']['inner_midpoint_squared_strict_upper']) == INNER_CAP**2,
            'certified midpoint caps for the prescribed inner and outer pairs')


def tensor_coefficient(quadratic, left, right):
    return sum((quadratic[i][j] for i, k in left for j, l in right if k == l),
               quadratic[0][0] * 0)


def check_jacobian(payload):
    source = payload['matrix']['entries']
    dim = len(source)
    integer_scale = 1 << 10
    rows = []
    for i, row in enumerate(source):
        rows.append([{'num': str((integer_scale*int(entry['num']))//int(entry['den'])
                                 - dim*(i == j)), 'den': '1'}
                     for j, entry in enumerate(row)])
    lower = {'matrix': {'domain': 'QQ', 'entries': rows}}
    serialized = json.dumps(lower, separators=(',', ':')).encode('ascii')
    result = read(HERE / 'jacobian_inertia.json')
    require(result['operation_id'] == 'matrix.inertia.compute'
            and result['source_sha256'] == sha256(serialized).hexdigest()
            and result['inertia'] == [dim,0,0], 'source-bound Jacobian integer lower matrix')
    return result['inertia'], lower


def bound_quartic(basin, representatives):
    zero = basin.ZERO
    exact = basin.exact
    scale = basin.SCALE
    coord = basin.eye(5)
    quadratic = [[10*s for s in row] for row in basin.eye(10)]
    nearest_displacement = F(10,7)*AMPLITUDE
    require(F(10,7)**2 > 2, 'safe uniform bound for amplitude displacement per edge')
    tensor_labels = list(combinations_with_replacement(range(10), 2))
    tensor_columns = [[(i,j)] if i == j else [(i,j),(j,i)] for i,j in tensor_labels]
    quartic = [[zero]*len(tensor_labels) for _ in tensor_labels]
    max_rho = F(0)
    for i, j in combinations(range(5), 2):
        eta_coeff = [x-y for x,y in zip(coord[i],coord[j])]
        midpoint_displacement = (INNER_CAP if i == 0 else OUTER_CAP) + OUTER_CAP
        for sign in (-1,1):
            d = tuple(x+sign*y for x,y in zip(representatives[i], representatives[j]))
            min_d = F(basin.sqrt(sum((x.square() for x in d), zero)).lo, scale)
            require(min_d > F(73,100), 'separated candidate cross edge')
            theta = 1-nearest_displacement/min_d
            rho = (midpoint_displacement/(min_d-nearest_displacement))**2
            require(theta > 0 and rho < 1, 'zero-free logarithm expansion on the full segment')
            max_rho = max(max_rho, rho)
            error = exact((1+2*rho/(3*(1-rho)))/theta**4-1)
            normalized = basin.normalized_functional(d, eta_coeff)
            quadratic = basin.add_quadratic(quadratic, normalized, 2)
            real, imag = [], []
            for (p, q), mult in zip(tensor_labels, map(len,tensor_columns)):
                px, py = normalized[p]
                qx, qy = normalized[q]
                real.append(mult*(px*qx-py*qy))
                imag.append(mult*(px*qy+py*qx))
            quartic = basin.add(quartic, basin.outer(real,real), beta=error-1)
            quartic = basin.add(quartic, basin.outer(imag,imag), beta=error+1)

    metric = [[tensor_coefficient(quadratic,l,r) for r in tensor_columns] for l in tensor_columns]
    matrix = basin.add(metric, quartic, alpha=exact(LIFT_BOUND), beta=-1)
    report, payload = basin.basin.rational_spd(matrix, scale, 'signed midpoint quartic')
    require(max_rho < F(1, 4), 'uniform interval logarithm ratio below one quarter')
    return report, payload


def angular_localization(repo, basin, representatives):
    gaps = read(repo / 'ten-point-outer-gap/verification_report.json')
    radii = read(repo / 'ten-point-angular-radii/verification_report.json')
    inner = read(repo / 'ten-point-inner-midpoint/verification_report.json')
    require(all(report['passed'] and report['applies_to_every_global_maximizer']
                and report['normalization_mass']==9 and report['unresolved']==0
                for report in (gaps, radii, inner)), 'certified angular and radial domain')
    lower,upper = map(F, gaps['outer_successive_angle_cosine_open_interval'])
    stiffness = F(inner['outer_mean_zero_negative_angle_hessian_lower_on_successive_gap_polytope'])
    radius = F(radii['outer_squared_radius_open_interval'][1])
    require(stiffness==F(161,50) and radius/2<F(4,5)**2, 'angular residual to half-difference coefficient')
    bx,by=representatives[1]
    cx,cy=representatives[3]
    require(0<bx.lo<bx.hi<by.lo and 0<cy.lo<cy.hi<cx.lo,
            'the larger-radius rectangle is nearer the real axis')
    b,c=bx.square()+by.square(),cx.square()+cy.square()
    successive=((bx*cx+by*cy)/basin.sqrt(b*c), (by.square()-bx.square())/b,
                (cx.square()-cy.square())/c)
    require(all(lower<F(v.lo,basin.SCALE)<=F(v.hi,basin.SCALE)<upper for v in successive),
            'all candidate successive gaps are in the convex angular domain')
    return {'candidate_successive_gap_types_verified':len(successive),
            'half_difference_error_coefficient_for_projected_angular_gradient':str(F(4,5)/stiffness),
            'free_shape_parameters_after_outer_angle_elimination':11}


def verify(repo):
    check_inputs(repo)
    basin = load(repo / 'ten-point-weighted-pair-basin/verify_all.py')
    require(basin.AMPLITUDE == AMPLITUDE, 'amplitudes in the known sharp basin')
    reps = basin.algebraic_data()
    basin.matrix_bounds(reps)
    basin.remainder_bounds()
    check, payload = bound_quartic(basin, reps)
    jacobian_inertia, jacobian_input = check_jacobian(payload)
    x = F(10,7)*AMPLITUDE/F(73,100)
    quadratic_error = (basin.CUBIC_MIDPOINT*AMPLITUDE
                       +4*AMPLITUDE**2*(3-2*x)*basin.LENGTH_MIDPOINT_RELATIVE_4/(1-x)**2)
    quarter_error = LIFT_BOUND*MIDPOINT**2
    margin = 1-quadratic_error-quarter_error
    require(margin > F(1,20), 'strict midpoint energy loss for the whole global midpoint range')
    return ({'passed':True, 'normalization_mass':9,
             'half_difference_radius':str(AMPLITUDE),
             'midpoint_norm_cap':str(MIDPOINT),
             'individual_midpoint_caps':[str(INNER_CAP),str(OUTER_CAP)],
             'signed_quartic_relative_upper':str(LIFT_BOUND),
             'relative_energy_loss_strict_lower':'1/20',
             'relative_energy_loss_certified_budget':str(margin),
             'signed_quartic_matrix':check,
             'jacobian_lower_matrix_inertia':jacobian_inertia,
             'angular_localization':angular_localization(repo, basin, reps)}, jacobian_input)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=HERE.parent)
    parser.add_argument('--record-report', action='store_true')
    parser.add_argument('--export-jacobian', type=Path)
    args = parser.parse_args()
    report, payload = verify(args.repo.resolve())
    saved = HERE / 'verification_report.json'
    content = json.dumps(report, indent=2)+'\n'
    if args.record_report:
        saved.write_text(content, encoding='utf-8', newline='\n')
    else:
        require(report == read(saved), 'exact replay matches the saved report')
    if args.export_jacobian:
        args.export_jacobian.write_text(json.dumps(payload,indent=2)+'\n', encoding='utf-8', newline='\n')
    print(content,end='')


if __name__ == '__main__':
    main()
