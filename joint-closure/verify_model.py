"""Exact arithmetic checks of the denominator-free physical/stability matrix."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'power-moment'))
from complex_interval import C


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def trim(p):
    p = list(p)
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def padd(p, q):
    out = list(p)+[C()] * max(0, len(q)-len(p))
    for i, v in enumerate(q):
        out[i] = out[i]+v
    return trim(out)


def pscale(p, a):
    return trim([v*a for v in p])


def pmul(p, q):
    out = [C()]*(len(p)+len(q)-1)
    for i, v in enumerate(p):
        for j, w in enumerate(q):
            out[i+j] = out[i+j]+v*w
    return trim(out)


def prem(p, modulus):
    p = trim(p)
    while len(p) >= len(modulus):
        k = len(p)-len(modulus)
        factor = p[-1]/modulus[-1]
        for j, value in enumerate(modulus):
            p[j+k] = p[j+k]-factor*value
        p = trim(p)
    return p


def psubstitute(p, x, modulus):
    out = [C()]
    for value in reversed(p):
        out = prem(padd(pmul(out, x), [value]), modulus)
    return out


def determinant(a):
    a = [list(row) for row in a]
    out = C(1)
    for j in range(len(a)):
        k = next((i for i in range(j, len(a)) if a[i][j]), None)
        if k is None:
            return C()
        if k != j:
            a[k], a[j] = a[j], a[k]
            out = -out
        pivot = a[j][j]
        out = out*pivot
        for i in range(j+1, len(a)):
            factor = a[i][j]/pivot
            for k in range(j+1, len(a)):
                a[i][k] = a[i][k]-factor*a[j][k]
    return out


def pure_sums(p):
    n = len(p)-1
    values = [C(n)]
    for k in range(1, n+1):
        value = -k*p[n-k]
        for j in range(1, k):
            value = value-p[n-j]*values[k-j]
        values.append(value)
    return values


def quotient_model(p, conjugate):
    """The conjugate map is supplied by an explicit actual root geometry."""
    sums = pure_sums(p)

    def trace(poly):
        poly = prem(poly, p)
        return sum((v*sums[k] for k, v in enumerate(poly)), C())

    def inner(f, g):
        barf = psubstitute([v.conjugate() for v in f], conjugate, p)
        return trace(pmul(barf, g))

    powers = [[C()]*k+[C(1)] for k in range(8)]
    action = [[C()], pscale(conjugate, 8)]
    for k in range(2, 8):
        value = pscale(pmul(powers[k-1], conjugate), 8*k)
        for r in range(k-1):
            value = padd(value, pscale(powers[k-2-r], -2*(k-1-r)*sums[r]))
        value = padd(value, pscale(powers[k-2], k*(k-1)))
        action.append(prem(value, p))
    H = [[inner(powers[j], powers[k]) for k in range(8)] for j in range(8)]
    K = [[64*H[j][k]-inner(action[j], action[k]) for k in range(8)] for j in range(8)]
    need(all(K[j][k] == K[k][j].conjugate() for j in range(8) for k in range(8)), 'Hermitian model')
    need(psubstitute(conjugate, conjugate, p) == [C(), C(1)], 'actual conjugation is involutive')
    need(all(not K[1][k] for k in range(8)), 'saturated row zero for these sparse/real models')
    return H, K


def model_checks():
    heptagon = [C(), C(-1)]+[C()]*6+[C(1)]
    H, K = quotient_model(heptagon, [C()]*6+[C(1)])
    need(H[1][1] == C(7), 'candidate mass seven')
    need(determinant(H) == C(7**7), 'candidate exact discriminant')
    minors = []
    for size in range(1, 9):
        for ix in combinations(range(8), size):
            value = determinant([[K[i][j] for j in ix] for i in ix])
            need(not value.im and value.re >= 0, 'candidate all principal minors')
            minors.append(value.re)
    octagon = [C(-F(7, 8)**4)]+[C()]*7+[C(1)]
    H8, K8 = quotient_model(octagon, [C()]*7+[C(F(512, 343))])
    need(H8[1][1] == C(7) and K8[4][4] == C(-F(735, 8)),
         'integer-scaled octagon instability benchmark')
    r = F(1001, 1000)
    scaled = [C(), C(-r**7)]+[C()]*6+[C(1)]
    Hr, Kr = quotient_model(scaled, [C()]*6+[C(1/r**5)])
    need(Hr[1][1] == C(7*r*r) and Hr[1][1] != C(7), 'wrong mass is explicit')
    for size in range(1, 9):
        for ix in combinations(range(8), size):
            value = determinant([[Kr[i][j] for j in ix] for i in ix])
            need(not value.im and value.re >= 0, 'mass-free matrix alone can accept wrong scaling')
    return {'candidate_principal_minors': len(minors), 'octagon_K44': str(K8[4][4].re),
            'mass_constraint_is_essential': True, 'wrong_mass_accepted_without_normalization': str(Hr[1][1].re)}


def nonstationary_check():
    base = [C(F(2, 3)), C(-F(2, 3)), C(0, F(5, 6)), C(0, -F(5, 6)),
            C(1), C(-1), C(0, F(7, 6)), C(0, -F(7, 6))]
    integer_target = [(1123, 0), (112, -1073), (427, -428), (112, 1073),
                      (-754, 727), (-694, 0), (427, 428), (-754, -727)]
    target = [C(F(a, 1000), F(b, 1000)) for a, b in integer_target]
    mean = sum(target, C())/8
    target = [z-mean for z in target]
    direction = [a-b for a, b in zip(base, target)]
    numerator = sum((a.conjugate()*b for a, b in zip(direction, base)), C()).re
    denominator = sum(v.norm() for v in direction)
    q = [z-2*numerator/denominator*v for z, v in zip(base, direction)]
    need(sum(q, C()) == C() and sum(z.norm() for z in q) == 7, 'actual rational centered mass seven')
    need(len(set(q)) == 8, 'eight distinct actual Gaussian-rational points')
    sums = [sum((z**k for z in q), C()) for k in range(9)]
    need(sums[2].norm()+F(54, 175)*sums[3].norm() < 1 and sums[2] and sums[3],
         'actual nonstationary witness is inside the nonzero-moment ellipse')
    action = [[C() for _ in range(8)] for _ in range(8)]
    for i in range(8):
        action[i][1] = 8*q[i].conjugate()
        for k in range(2, 8):
            action[i][k] = 8*k*q[i]**(k-1)*q[i].conjugate()+k*(k-1)*q[i]**(k-2)
            for r in range(k-1):
                action[i][k] = action[i][k]-2*(k-1-r)*sums[r]*q[i]**(k-2-r)
    K = [[64*sum((q[i]**j).conjugate()*q[i]**k for i in range(8))
          -sum(action[i][j].conjugate()*action[i][k] for i in range(8)) for k in range(8)] for j in range(8)]
    need(K[1][1] == C(), 'zero diagonal is an exact identity on actual coordinates')
    residual = [sum((C(1)/(q[i]-q[j]) for j in range(8) if j != i), C())
                -4*q[i].conjugate() for i in range(8)]
    nonzero = []
    for k in range(1, 7):
        Q = sum((z.norm()*z**k for z in q), C())
        G = 8*Q-(15-k)*sums[k]-sum((sums[r]*sums[k-r] for r in range(1, k)), C())
        need(K[1][k+1] == -8*k*G, 'zero-row exact Ward identity '+str(k))
        projected = sum((q[i]**(k+1)*residual[i] for i in range(8)), C())
        need(projected == -G/2, 'pointwise force projection identity '+str(k))
        if G:
            minor = -K[1][k+1].norm()
            need(minor < 0, 'nonstationarity rejected by a full matrix principal minor')
            nonzero.append(k)
    need(nonzero, 'a nonstationary actual input is not silently accepted')
    return {'actual_point_center_and_mass_exact': True, 'inside_old_ellipse': True,
            'both_low_moments_nonzero': True, 'zero_row_identities_checked': 6,
            'direct_force_projection_identities_checked': 6,
            'negative_minor_Ward_indices': nonzero}


def interval_crosschecks():
    import interval_engine as E
    E.setup(60)
    q = [C(F(2, 3)), C(-F(2, 3)), C(0, F(5, 6)), C(0, -F(5, 6)),
         C(1), C(-1), C(0, F(7, 6)), C(0, -F(7, 6))]
    free = [q[0].re]+[v for z in q[1:7] for v in [z.re, z.im]]
    box = [E.rat(v.numerator, v.denominator) for v in free]
    qi = E.positions(box)
    inv, _ = E.invert_edges(qi)

    def contains(interval, value):
        return interval[0]*value.denominator <= value.numerator*E.S <= interval[1]*value.denominator

    for (i, j), interval in zip(E.PAIRS, inv):
        exact = C(1)/(q[i]-q[j])
        need(contains(interval[0], exact.re) and contains(interval[1], exact.im),
             'Gaussian-rational inverse lies in the integer enclosure')
    residual = [sum((C(1)/(q[i]-q[j]) for j in range(8) if j != i), C())
                -4*q[i].conjugate() for i in range(8)]
    Jpoly, Jforce = E.polynomial_jacobian(qi), E.force_jacobian(qi)
    checked = 0
    for a in range(13):
        h = [C(*E.L[i][a]) for i in range(8)]
        Bh = [sum(((h[i]-h[j])/(q[i]-q[j])**2 for j in range(8) if j != i), C())
              for i in range(8)]
        dr = [-Bh[i]-4*h[i].conjugate() for i in range(8)]
        derivatives = [2*sum((q[i].conjugate()*h[i] for i in range(8)), C()).re]
        for k in range(1, 7):
            # Differentiate the DIRECT RATIONAL FORCE projection, not the polynomial implementation.
            dG = -2*sum(((k+1)*q[i]**k*h[i]*residual[i]+q[i]**(k+1)*dr[i]
                         for i in range(8)), C())
            derivatives.extend([dG.re, dG.im])
        for row, value in enumerate(derivatives):
            need(contains(Jpoly[row][a], value), 'polynomial Jacobian contains independent force derivative')
            checked += 1
        dg = [4*h[i]+Bh[i].conjugate() for i in range(8)]
        for b in range(13):
            value = sum((C(E.L[i][b][0], -E.L[i][b][1])*dg[i]
                         for i in range(8)), C()).re
            need(contains(Jforce[b][a], value), 'real Hessian contains independent directional derivative')
            checked += 1
    displacement = [[i-3, 2*i-5] for i in range(8)]
    h = [C(*z) for z in displacement]
    exact = 4*sum(z.norm() for z in h)+sum((((h[i]-h[j])/(q[i]-q[j]))**2
                for i, j in E.PAIRS), C()).re
    bound = E.rayleigh(qi, inv, displacement)
    need(contains(bound, exact), 'full-Hessian integer direction encloses an independent exact value')
    return {'exact_Gaussian_inverse_checks': 28, 'independent_Jacobian_and_Hessian_entry_checks': checked,
            'full_Hessian_direction_crosschecked': True,
            'floating_tolerances_used_for_acceptance': False}


def coverage_guard_checks():
    import tempfile
    import interval_engine as E
    from check_cover import check
    E.setup(40)
    settings = E.model('global_maximum')
    toy = {'precision_bits': 40, 'settings': settings, 'nodes': [{'kind': 'unknown'}],
           'pending': [[0, settings['initial_box']]], 'complete_cover': False,
           'full_eight_point_exclusion_proved': False}
    rejected = 0
    with tempfile.TemporaryDirectory(prefix='joint-guard-', dir=ROOT) as temp:
        path = Path(temp)/'cover.json'
        path.write_text(json.dumps(toy))
        result = check(path)
        need(result['passed'] and not result['complete_cover'], 'unknown root is accepted only as incomplete')
        for kind in range(3):
            altered = json.loads(json.dumps(toy))
            if kind == 0:
                altered['full_eight_point_exclusion_proved'] = True
            elif kind == 1:
                altered['pending'] = []
            else:
                altered['pending'] = []
                altered['nodes'] = [{'kind': 'excluded', 'reason': 'force_sign'}]
                altered['complete_cover'] = altered['full_eight_point_exclusion_proved'] = True
            path.write_text(json.dumps(altered))
            try:
                check(path)
            except ArithmeticError:
                rejected += 1
            else:
                raise ArithmeticError('false completeness or false exclusion was accepted')
    need(rejected == 3, 'all dishonest coverage variants rejected')
    return {'unknown_leaf_is_not_treated_as_empty': True, 'rejected_false_claims': rejected}


def main():
    result = {'passed': True, 'model_checks': model_checks(),
              'actual_nonstationary_witness': nonstationary_check(),
              'integer_interval_crosschecks': interval_crosschecks(),
              'coverage_guard_checks': coverage_guard_checks(),
              'polynomial_matrix_degree_upper': 14, 'full_eight_point_exclusion_proved': False}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
