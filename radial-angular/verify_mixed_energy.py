"""Exact acceptance for radial-weighted angular projection certificates."""
from mixed_energy import *
from math import factorial
import json


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def multiply(a, b):
    out = {}
    for (i, j), c in a.items():
        for (p, q), d in b.items():
            key = (i+p, j+q)
            z = cmul(c, d)
            old = out.get(key, (F(0), F(0)))
            out[key] = (old[0]+z[0], old[1]+z[1])
    return out


def independent_disk_moment(x, k, j, r2):
    # Cartesian expansion of (x+u)^(k+j) conjugate(x+u)^j,
    # integrated using independent exact Cartesian uniform-disk moments.
    re, im = map(F, x)
    z = {(0, 0): (re, im), (1, 0): (F(1), F(0)), (0, 1): (F(0), F(1))}
    zb = {(0, 0): (re, -im), (1, 0): (F(1), F(0)), (0, 1): (F(0), F(-1))}
    poly = {(0, 0): (F(1), F(0))}
    for _ in range(k+j):
        poly = multiply(poly, z)
    for _ in range(j):
        poly = multiply(poly, zb)
    total = [F(0), F(0)]
    for (a, b), coeff in poly.items():
        if a % 2 or b % 2:
            continue
        u, v = a//2, b//2
        weight = (r2**(u+v)*F(factorial(2*u)*factorial(2*v),
                    4**(u+v)*factorial(u)*factorial(v)*factorial(u+v+1)))
        total[0] += weight*coeff[0]
        total[1] += weight*coeff[1]
    return tuple(total)


def main():
    disk_checks = 0
    for x in ((F(2, 3), F(1, 4)), (F(0), F(1)), (F(-2), F(0))):
        for k in range(1, 6):
            power = (F(1), F(0))
            for _ in range(k):
                power = cmul(power, x)
            for j in range(4):
                r2 = F(2, 7)
                independent = independent_disk_moment(x, k, j, r2)
                w = disk_weight(k, j, x[0]**2+x[1]**2, r2)
                require(independent == (power[0]*w, power[1]*w), 'independent Cartesian disk mean')
                disk_checks += 1
    gram_checks = 0
    for k in range(1, 9):
        for degree in range(5):
            s = F(9, 4)
            h = gram_matrix(k, degree, s)
            for j in range(degree+1):
                for l in range(degree+1):
                    a, b = k+2*j, k+2*l
                    independent = (F(a*b+k*k, 2*(a+b))+F(k, 2))*s**(k+j+l)
                    require(independent == h[j][l], 'inside plus outside Gram integral')
            rhs = [F(j+1, j+2) for j in range(degree+1)]
            solved = solve(h, rhs)
            require(all(sum((h[j][l]*solved[l] for l in range(degree+1)), F(0)) == rhs[j]
                        for j in range(degree+1)), 'positive-pivot exact solve')
            gram_checks += 1
        # Closed Schur formula is independently compared to exact matrix solve.
        p, q, u, v = F(2, 3), F(-4, 5), F(7, 9), F(3, 8)
        value = projection_value(k, s, [p, u], [q, v])
        schur = (p*p+q*q)/(k*s**k)+(k+2)*((u-s*p)**2+(v-s*q)**2)/s**(k+2)
        require(value == schur, 'explicit first mixed level')
    example = [(3, 0), (-3, 0), (4, 0), (-4, 0), (0, 5), (0, -5)]
    cert = shape_certificate(example, support_squared=3)
    require(cert['holomorphic_delta_lower'] == 0, 'vanishing first three power sums')
    require(cert['delta_lower'] == F(64, 625), 'hidden radial-angular mode')
    # Dual certificate attains this projection for the exact second-mode vector.
    h = gram_matrix(2, 1, 3)
    actual_re, actual_im = [F(0), F(-36, 25)], [F(0), F(0)]
    ar, ai = solve(h, actual_re), solve(h, actual_im)
    coeffs = list(zip(ar, ai))
    rectangles = [(a, a, b, b) for a, b in zip(actual_re, actual_im)]
    dual = dual_interval_bound(2, 3, coeffs, rectangles)
    require(dual == F(64, 625), 'exact dual equality')
    widened = [(a-F(1, 1000), a+F(1, 1000), b-F(1, 1000), b+F(1, 1000))
               for a, b in zip(actual_re, actual_im)]
    require(0 < dual_interval_bound(2, 3, coeffs, widened) < dual, 'outward whole-class intervals')
    shapes = [example, [(0, 0), (1, 0), (0, 2)],
              [(1, 0), (0, 1), (-1, 0), (0, -1)],
              [(i, i*i % 5) for i in range(8)]]
    nested = 0
    for points in shapes:
        prior = F(-1)
        for degree in range(4):
            data = shape_certificate(points, max_mode=4, radial_degree=degree)
            require(data['angular_delta_lower'] >= prior, 'nested radial-degree projections')
            prior = data['angular_delta_lower']
            n, pairs = len(points), len(points)*(len(points)-1)//2
            sx = F(sum(a for a, b in points), n)
            sy = F(sum(b for a, b in points), n)
            q = sum((F(a)-sx)**2+(F(b)-sy)**2 for a, b in points)
            raw = F(1)
            for i, (a, b) in enumerate(points):
                for c, d in points[i+1:]:
                    raw *= (a-c)**2+(b-d)**2
            delta = F(n)**pairs*raw/q**pairs
            actual_gap_lower = (F(n, 2)*log_bounds(n)[0]
                                 +pairs*log_bounds(F(2*n, n-1))[0]
                                 -F(pairs, 2)-log_bounds(delta)[1])
            require(data['delta_lower'] < actual_gap_lower, 'certificate below actual shape gap')
            nested += 1
    for hval in (F(1, 2), F(99, 100), F(1), F(3, 2)):
        got = shape_certificate(example, h=hval)
        require(got['smear_penalty_upper'] >= 0, 'correct sign of adjustable-smear subtraction')
    bad = 0
    for call in (lambda: gram_matrix(0, 1, 2), lambda: gram_matrix(2, -1, 2),
                 lambda: gram_matrix(2, 1, 0), lambda: shape_certificate(example, h=0),
                 lambda: shape_certificate(example, support_squared=2),
                 lambda: shape_certificate([(0, 0), (0, 0)]),
                 lambda: dual_interval_bound(2, 3, [(1, 0)], [(1, 0, 0, 0)])):
        try:
            call()
        except ValueError:
            bad += 1
        else:
            raise ArithmeticError('invalid input accepted')
    report = {'passed': True, 'independent_cartesian_disk_means': disk_checks,
              'gram_systems_checked': gram_checks, 'nested_actual_shape_checks': nested,
              'invalid_inputs_rejected': bad, 'hidden_mode_example': {
                  'points': example, 'vanished_power_sums': [1, 2, 3],
                  'old_same_mode_sum': '0', 'new_mixed_sum': '64/625'},
              'dual_exact_and_interval_checks': 2,
              'scope': 'Finite algebraic implementation checks; universal analytic theorem is MIXED_ENERGY.md. No optimality or full-domain cover.'}
    Path(__file__).with_name('mixed_energy_verification.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print('PASS mixed angular Gram, disk means, nested projections and dual intervals')


if __name__ == '__main__':
    main()
