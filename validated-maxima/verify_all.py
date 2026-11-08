"""Verify existence and full second variation using integers only; -O safe.

The analytic gauge, normalization and constrained-inertia argument is in PROOF.md.
No floating residual, numerical eigenvalue or stored success flag is accepted.
"""
from fractions import Fraction
from pathlib import Path
import argparse, hashlib, json
from dyadic import Ball, PRECISION, SCALE, matmul, transpose

HERE = Path(__file__).resolve().parent


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def equations(variables, n):
    # y_0=0 is the rotation gauge; retain every other force component.
    need(len(variables) == 2*n-1, 'rotation-gauge dimension')
    coordinates = variables[:1]+[Ball.integer(0)]+variables[1:]
    force = [-((n-1)*q/2) for q in coordinates]
    jacobian = [[Ball.integer(0) for _ in coordinates] for _ in coordinates]
    for i in range(2*n):
        jacobian[i][i] = Ball.integer(-(n-1))/2
    product = Ball.integer(1)
    for i in range(n):
        for j in range(i):
            dx, dy = coordinates[2*i]-coordinates[2*j], coordinates[2*i+1]-coordinates[2*j+1]
            d2 = dx.square()+dy.square()
            need(d2.lo > 0, 'entire box is collision-free')
            product = product*d2
            inverse = 1/d2
            for k, diff in enumerate((dx, dy)):
                value = diff*inverse
                force[2*i+k] = force[2*i+k]+value
                force[2*j+k] = force[2*j+k]-value
            inverse2 = inverse.square()
            block = [[inverse-2*dx.square()*inverse2, -2*dx*dy*inverse2],
                     [-2*dx*dy*inverse2, inverse-2*dy.square()*inverse2]]
            for a in range(2):
                for b in range(2):
                    v = block[a][b]
                    for r, c, sign in ((2*i+a,2*i+b,1),(2*j+a,2*j+b,1),
                                       (2*i+a,2*j+b,-1),(2*j+a,2*i+b,-1)):
                        jacobian[r][c] = jacobian[r][c]+sign*v
    keep = [i for i in range(2*n) if i != 1]
    return [force[i] for i in keep], [[jacobian[i][j] for j in keep] for i in keep], product


def point_matrix(numbers, exponent, d):
    need(len(numbers) == d and all(len(row) == d for row in numbers), 'square certificate matrix')
    return [[Ball.dyadic(q, exponent) for q in row] for row in numbers]


def rational_floor(value, bits=40):
    return str(Fraction(value//(1 << (PRECISION-bits)), 1 << bits))


def power_upper(value):
    return '0' if value == 0 else '2^'+str((value-1).bit_length()-PRECISION)


def independent_checks():
    intervals = [(Fraction(a,8),Fraction(b,8)) for a in range(-4,5) for b in range(a,5)]
    checked = 0
    for a,b in intervals:
        x = Ball.dyadic(a.numerator, a.denominator.bit_length()-1)
        y = Ball.dyadic(b.numerator, b.denominator.bit_length()-1)
        u = Ball(x.lo,y.hi)
        sq = u.square()
        low = 0 if a <= 0 <= b else min(a*a,b*b)
        need(Fraction(sq.lo,SCALE) <= low and Fraction(sq.hi,SCALE) >= max(a*a,b*b), 'independent square enclosure')
        for c,d in intervals:
            l = Ball.dyadic(c.numerator,c.denominator.bit_length()-1)
            h = Ball.dyadic(d.numerator,d.denominator.bit_length()-1)
            v = Ball(l.lo,h.hi)
            for value, corners in ((u+v,[a+c,b+d]), (u*v,[a*c,a*d,b*c,b*d])):
                need(Fraction(value.lo,SCALE) <= min(corners)
                     and Fraction(value.hi,SCALE) >= max(corners), 'independent rational endpoint enclosure')
                checked += 1
            if not c <= 0 <= d:
                value = u/v; corners = [a/c,a/d,b/c,b/d]
                need(Fraction(value.lo,SCALE) <= min(corners)
                     and Fraction(value.hi,SCALE) >= max(corners), 'independent signed division enclosure')
                checked += 1
    # Independently construct the complex inverse-square Laplacian with Fraction.
    # Its conjugate-real action must enclose every real Jacobian entry.
    complex_entries = 0
    for n in (3,4,5,6):
        coords = [(Fraction(1),Fraction(0))]+[(Fraction(i+1,3),Fraction((i*i)%7+1,5)) for i in range(1,n)]
        flat = [v for xy in coords for v in xy]
        variables = [Ball.integer(q.numerator)/q.denominator for i,q in enumerate(flat) if i != 1]
        f,jac,prod = equations(variables,n)
        A = [[(Fraction(0),Fraction(0)) for _ in range(n)] for _ in range(n)]
        reference_force = [[-(n-1)*a/2,-(n-1)*b/2] for a,b in coords]
        delta = Fraction(1)
        for i in range(n):
            for j in range(i):
                dx,dy = coords[i][0]-coords[j][0],coords[i][1]-coords[j][1]
                s = dx*dx+dy*dy;delta *= s
                real,imag = (dx*dx-dy*dy)/(s*s), -2*dx*dy/(s*s)
                for r,c,sign in ((i,i,1),(j,j,1),(i,j,-1),(j,i,-1)):
                    A[r][c] = (A[r][c][0]+sign*real,A[r][c][1]+sign*imag)
                for k,diff in enumerate((dx,dy)):
                    reference_force[i][k] += diff/s
                    reference_force[j][k] -= diff/s
        keep = [i for i in range(2*n) if i != 1]
        for ri,i in enumerate(keep):
            force = reference_force[i//2][i%2]
            need(Fraction(f[ri].lo,SCALE) <= force <= Fraction(f[ri].hi,SCALE), 'exact independent force')
            for cj,j in enumerate(keep):
                a,b = A[i//2][j//2]
                value = ((-a,b),(b,a))[i%2][j%2]
                if i == j:
                    value -= Fraction(n-1,2)
                need(Fraction(jac[ri][cj].lo,SCALE) <= value <= Fraction(jac[ri][cj].hi,SCALE),
                     'complex Laplacian vs real Jacobian')
                complex_entries += 1
        need(Fraction(prod.lo,SCALE) <= delta <= Fraction(prod.hi,SCALE), 'independent distance product')
    try:
        Ball.integer(1)/Ball(-1,1)
    except ZeroDivisionError:
        pass
    else:
        raise ArithmeticError('division through zero accepted')
    return {'rational_endpoint_operations':checked, 'independent_complex_laplacian_entries':complex_entries,
            'division_through_zero_rejected':True}


def verify_row(row, common):
    n = row['n']; d = 2*n-1
    need(type(n) is int and 2 <= n <= 100, 'supported positive n')
    radius_exp = common['radius_exponent']
    need(0 < radius_exp < PRECISION, 'positive representable cube radius')
    radius = 1 << (PRECISION-radius_exp)
    mid = [Ball.dyadic(q, common['center_exponent']) for q in row['center_numerators']]
    need(len(mid) == d and all(q.lo == q.hi for q in mid), 'exact dyadic midpoint')
    box = [Ball(q.lo-radius, q.hi+radius) for q in mid]
    need(box[0].lo > 0, 'rotation gauge x_0>0 throughout the box')
    fmid, _, _ = equations(mid, n)
    _, jac, delta = equations(box, n)
    C = point_matrix(row['inverse_numerators'], common['matrix_exponent'], d)
    residual = matmul(C, [[q] for q in fmid])
    eta = max(q[0].magnitude() for q in residual)
    cj = matmul(C, jac)
    contraction = max(sum((Ball.integer(int(i == j))-cj[i][j]).magnitude() for j in range(d)) for i in range(d))
    need(4*eta < radius, 'preconditioned exact midpoint residual < radius/4')
    need(4*contraction < SCALE, 'uniform contraction norm < 1/4')
    result = {'n':n, 'passed':True, 'gauge_dimension':d,
              'unique_stationary_orbit_in_certified_slice_box':True,
              'box_radius':'2^-'+str(radius_exp),
              'preconditioned_residual_upper':power_upper(eta),
              'contraction_norm_upper':power_upper(contraction),
              'classification':row['classification'],
              'delta_interval':[str(Fraction(delta.lo,SCALE)),str(Fraction(delta.hi,SCALE))]}
    if row['classification'] == 'strict_local_maximum_mod_rotation':
        W = point_matrix(row['whitener_numerators'], common['matrix_exponent'], d)
        congruence = matmul(transpose(W), matmul([[-q for q in line] for line in jac], W))
        gap = min(congruence[i][i].lo-sum(congruence[i][j].magnitude() for j in range(d) if i != j) for i in range(d))
        need(gap > 0, 'full reduced negative Hessian: strict Gershgorin congruence')
        result['whitened_negative_hessian_lower'] = rational_floor(gap)
    elif row['classification'] == 'unstable_stationary_orbit':
        vector = [[Ball.dyadic(q, common['matrix_exponent'])] for q in row['positive_vector_numerators']]
        need(len(vector) == d, 'instability direction dimension')
        value = matmul(transpose(vector), matmul(jac, vector))[0][0]
        need(value.lo > 0, 'positive unrestricted second-variation direction')
        result['positive_quadratic_form_lower'] = rational_floor(value.lo, 80)
    else:
        raise ArithmeticError('unknown inertia claim')
    return result


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    for name, expected in manifest['files'].items():
        need('/' not in name and '\\' not in name and name not in ('.','..'), 'local fixed file path')
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, 'fixed input '+name)
    certificates = json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    need(certificates['arithmetic_precision'] == PRECISION, 'fixed integer interval precision')
    checks = independent_checks()
    rows = [verify_row(row, certificates) for row in certificates['rows']]
    report = {'passed':True, 'arithmetic':'outward integer dyadic intervals at 384 bits',
              'checks_survive_optimization':True, 'stationary_orbits_certified':len(rows),
              'independent_arithmetic_checks':checks,
              'global_optimality_claimed':False, 'rows':rows,
              'source_sha256':manifest['files']}
    if args.record_report:
        (HERE/'verification_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    run()
