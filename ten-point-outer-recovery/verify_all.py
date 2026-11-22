"""Recover the inner pair from eight outer points and bound its perturbation."""

from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def inputs():
    spec = importlib.util.spec_from_file_location('recovery_anchors', ROOT/'ten-point-phase-sweep/verify_all.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.AngularBox()


def verify_vandermonde():
    record = json.loads((HERE/'jacobian_vandermonde.json').read_bytes())
    require(record['operation_id'] == 'matrix.symbolic.determinant.compute', 'Jacobian operation')
    determinant = record['output']['determinant']
    unit = (0, 0, 0, 0)
    product = {unit: 1}
    for right in range(4):
        for left in range(right):
            next_terms = {}
            for powers, coeff in product.items():
                for index, sign in ((right, 1), (left, -1)):
                    key = tuple(e+(k == index) for k, e in enumerate(powers))
                    next_terms[key] = next_terms.get(key, 0)+sign*coeff
            product = {k: v for k, v in next_terms.items() if v}
    def decode(poly):
        return {tuple(row['exponents']): F(int(row['coefficient']['num']), int(row['coefficient']['den']))
                for row in poly['terms']}
    require(determinant['variables'] == ['z0', 'z1', 'z2', 'z3'] and
            decode(determinant['numerator']) == product and decode(determinant['denominator']) == {unit: 1},
            'four remaining forces have the Vandermonde determinant')


def bound(anchors, mode, radius):
    op = anchors.saddle
    exact, zero = anchors.exact, anchors.zero
    root = anchors.points[mode]
    t_star = op.norm2(root[0])
    require(root[0][0].lo > 0 and root[0][1].lo == root[0][1].hi == 0, 'real inner axis')
    require(anchors.endpoint(t_star) < F(53, 400), 'anchor half-difference')
    # The real Jacobian of delta H + 2 delta B s_star at the symmetric anchor.
    columns = []
    H = (zero, zero)
    for point in root[2:]:
        r = op.norm2(point)
        require(anchors.endpoint(r) < F(1271, 1000), 'anchor outer squared radius')
        H = op.add(H, op.scale_complex(op.mul(point, point), 17-10*r))
        z3 = op.mul(op.mul(point, point), point)
        linear = op.scale_complex(point, 34-30*r)
        for unit in [(exact(1), zero), (zero, exact(1))]:
            u, v = unit
            entry = op.sub(op.mul(linear, unit), op.scale_complex(op.mul(z3, (u, -v)), 10))
            radial = point[0]*u+point[1]*v
            columns.append((entry[0]+20*t_star*radial, entry[1]))
    target_error = op.add(H, (2*(17-10*t_star)*t_star, zero))
    require(all(v.lo <= 0 <= v.hi for v in target_error), 'outer-only Ward identity at the anchor')
    xx = sum((u.square() for u, v in columns), zero)
    xy = sum((u*v for u, v in columns), zero)
    yy = sum((v.square() for u, v in columns), zero)
    spectral = anchors.sqrt(anchors.Ball(0, max(xx.hi, yy.hi)+xy.magnitude()))
    jacobian_bound = -((-spectral.hi*1000)//anchors.scale)
    # Uniform absolute quadratic remainder, factored by ||D||^2.
    outer_squared = F(1271, 1000)
    outer_modulus = F(1128, 1000)
    require(outer_modulus**2 > outer_squared, 'outer anchor modulus')
    t_max = F(53, 400)
    coeff_A = F(17)   # A=17-10|m|^2-30t is positive and <= 17.
    require(17-10*F(87, 1000)**2-30*t_max > 0, 'coefficient A')
    b_min = 17-30*F(87, 1000)**2-10*t_max
    require(b_min > F(77, 5), 'uniform reconstruction denominator')
    quadratic = (30*outer_squared-17 + 30*outer_squared + 40*outer_modulus*radius
                 + 10*radius**2 + 2*35*t_max + 4*coeff_A)
    eta = (F(jacobian_bound, 1000)*radius+quadratic*radius**2)/(2*b_min)
    require(anchors.endpoint(t_star, False) > eta, 'square root stays in the right half-plane')
    inner_error = exact(eta)/(root[0][0] + anchors.sqrt(t_star-exact(eta)))
    if mode == 'candidate':
        squared = inner_error.square()+exact(radius**2/2)
        limit = F(7, 250)
    else:
        squared = 2*inner_error.square()+exact(5*radius**2)
        limit = F(21, 500)
    resulting_norm = anchors.sqrt(squared)
    require(anchors.endpoint(resulting_norm) < limit, mode+' basin/exclusion test')
    return {'mode': mode, 'outer_euclidean_radius': str(radius),
            'ward_linear_norm_upper': str(F(jacobian_bound, 1000)),
            'resulting_norm_upper_millionths': -((-resulting_norm.hi*10**6)//anchors.scale),
            'acceptance_radius': str(limit)}


def main():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for directory, items in ((HERE, manifest['files']), (ROOT, manifest['dependencies'])):
        for name, digest in items.items():
            require(sha256((directory/name).read_bytes()).hexdigest() == digest, 'hash '+name)
    verify_vandermonde()
    anchors = inputs()
    results = [bound(anchors, 'candidate', F(1, 70)), bound(anchors, 'saddle', F(1, 80))]
    print(json.dumps({'passed': True, 'mass': 9, 'outer_points': 8,
                      'outer_force_equations_complex': 6, 'outer_ward_equations_real': 3,
                      'outer_neighborhoods': results}, indent=2))


if __name__ == '__main__':
    main()
