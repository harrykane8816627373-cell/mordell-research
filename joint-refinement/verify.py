"""Exact checks of rotation saturation, logarithm arithmetic and continuation guards."""
from pathlib import Path
from fractions import Fraction as F
import gzip
import json
import sys
import tempfile
import refinement as R
from check import check

E = R.E
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'power-moment'))
from complex_interval import C


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def force(q):
    return [sum((1/(z-w) for j, w in enumerate(q) if i != j), C()) for i, z in enumerate(q)]


def H_action(q, h):
    return [4*h[i]+sum(((h[i]-h[j])/(q[i]-q[j])**2
        for j in range(len(q)) if i != j), C()).conjugate() for i in range(len(q))]


def inner(a, b):
    return sum((x.conjugate()*y for x, y in zip(a, b)), C()).re


def rotation_checks():
    q0 = [C(1), C(-1), C(0, 1), C(0, -1), C(1, F(1, 2)),
          C(-1, -F(1, 2)), C(F(1, 2)), C(-F(1, 2))]
    q1 = list(q0)
    q1[0], q1[1] = F(3, 5)*q0[0]-F(4, 5)*q0[1], F(4, 5)*q0[0]+F(3, 5)*q0[1]
    rows = []
    for q in [q0, q1]:
        mass = sum(z.norm() for z in q)
        need(mass == 7 and len(set(q)) == 8, 'actual distinct mass-seven rational points')
        f = force(q)
        g = [4*z-v.conjugate() for z, v in zip(q, f)]
        u, v = [C(0, 1)*z for z in q], [C(0, 1)*z for z in g]
        Hu = H_action(q, u)
        need(Hu == v and inner(u, Hu) == 0, 'rotation row has zero quadratic form and equals rotated gradient')
        residual2 = sum(z.norm() for z in g)
        need(residual2 > 0 and sum(z.norm() for z in f) == 112+residual2,
             'inverse-square norm lower bound with exact force-residual square')
        row_budget = max(sum((F(1, (q[i]-q[j]).norm()) for j in range(8) if i != j), F())
                         for i in range(8))
        bound = 4+2*row_budget
        epsilon = 1/bound
        h = [a-epsilon*b for a, b in zip(u, v)]
        quadratic = inner(h, H_action(q, h))
        Hv = inner(v, H_action(q, v))
        need(quadratic == -2*epsilon*residual2+epsilon**2*Hv < 0,
             'explicit rational negative-curvature direction at a nonstationary normalized point')
        rows.append({'centered': sum(q, C()) == C(), 'residual_squared': str(residual2),
            'negative_quadratic': str(quadratic)})
    far = [100*z for z in q0]
    row_bound = 2*max(sum((F(1, (far[i]-far[j]).norm()) for j in range(8) if i != j), F())
                      for i in range(8))
    need(row_bound < 4 and sum(z.norm() for z in far) == 70000,
         'wrong-mass actual configuration has strictly positive full Hessian')
    return {'two_actual_normalized_witnesses': rows,
        'centering_not_required_for_rotation_theorem': True,
        'wrong_mass_strict_positive_Hessian_witness_mass': 70000,
        'wrong_mass_inverse_square_operator_bound': str(row_bound)}


def reference_log(a, b):
    ratio, exponent = F(a, b), 0
    while ratio < 1:
        ratio *= 2
        exponent -= 1
    while ratio >= 2:
        ratio /= 2
        exponent += 1

    def reduced(z):
        x = (z-1)/(z+1)
        m = 64
        total = 2*sum((x**(2*k+1)/F(2*k+1) for k in range(m)), F())
        tail = 2*x**(2*m+1)/((2*m+1)*(1-x*x))
        return total, total+tail
    lo, hi = reduced(ratio)
    l2, u2 = reduced(F(2))
    return (lo+exponent*l2, hi+exponent*u2) if exponent >= 0 else (lo+exponent*u2, hi+exponent*l2)


def logarithm_checks():
    E.setup(40)
    pairs = [(1, 1), (2, 1), (1, 2), (7, 10), (11, 7), (100, 1),
             (1, 10**100), (10**100, 1), (2**127-1, 2**100), (7**7, 2**30)]
    for a, b in pairs:
        lo, hi = R.log_ratio(a, b)
        reflo, refhi = reference_log(a, b)
        need(F(lo, E.S) <= reflo <= refhi <= F(hi, E.S),
             'integer logarithm encloses independent high-precision Fraction series')
    return {'independent_fraction_logarithm_cases': len(pairs), 'floating_point_acceptance': False}


def fixed(v):
    return E.rat(v.numerator, v.denominator)


def exact_log_product(q):
    product = F(1)
    for i, j in E.PAIRS:
        product *= (q[i]-q[j]).norm()
    return reference_log(product.numerator, product.denominator*7**7)


def Taylor_checks():
    E.setup(40)
    q = [C(1), C(-1), C(0, 1), C(0, -1), C(1, F(1, 2)),
         C(-1, -F(1, 2)), C(F(1, 2)), C(-F(1, 2))]
    rotated = list(q)
    n = 10001
    phase = C(F(n*n-1, n*n+1), F(2*n, n*n+1))
    rotated[4], rotated[5] = q[4]*phase, q[5]*phase
    need(sum(z.norm() for z in rotated) == 7 and sum(rotated, C()) == C(),
         'second actual point has exact center and mass')
    data = [[q[0].re]+[v for z in q[1:7] for v in [z.re,z.im]],
            [rotated[0].re]+[v for z in rotated[1:7] for v in [z.re,z.im]]]
    box = [(min(fixed(a)[0],fixed(b)[0]), max(fixed(a)[1],fixed(b)[1])) for a,b in zip(*data)]
    center = [int(a*E.S) for a in data[0]]
    # A second, deliberately outside center tests the reference-center form.
    outside = [a+(E.S//1024 if i%3 == 0 else -E.S//2048) for i,a in enumerate(center)]
    bounds = []
    for c in [center, outside]:
        result = R.target_log_upper(box, E.positions(box), E.model()['radius'], c)
        need(result is not None and result['log_Delta_over_target_upper'] < 0,
             'mass-adjusted Taylor bound excludes this entire normalized box')
        for actual in [q, rotated]:
            lower, _ = exact_log_product(actual)
            need(lower <= F(result['log_Delta_over_target_upper'], E.S),
                 'direct exact distance product is bounded by the Taylor certificate')
        bounds.append(str(F(result['log_Delta_over_target_upper'], E.S)))
    wide = [(a-10*E.S,a+10*E.S) for a in center]
    need(R.target_log_upper(wide, E.positions(wide), E.model()['radius'], center) is None,
         'a Taylor disk reaching a pole is never treated as an exclusion')
    return {'actual_mass_seven_points_checked': 2, 'arbitrary_reference_centers_checked': 2,
        'whole_box_log_upper_bounds': bounds, 'wide_pole_box_rejected_as_inconclusive': True}


def guard_checks():
    base = json.loads(gzip.decompress((ROOT/'joint-closure/certificates/partial_cover.json.gz').read_bytes()))
    nodes = [{'kind':'unknown'} for _ in base['pending']]
    toy = {'base_sha256':R.BASE_SHA256,'precision_bits':40,'settings':base['settings'],
        'roots':[[i,j] for j,(i,b) in enumerate(base['pending'])], 'nodes':nodes,
        'pending':[[j,b] for j,(i,b) in enumerate(base['pending'])],
        'complete_cover':False,'full_eight_point_exclusion_proved':False}
    rejected=0
    with tempfile.TemporaryDirectory(prefix='refinement-guard-',dir=ROOT) as temp:
        path=Path(temp)/'toy.gz'

        def save(data):
            path.write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        save(toy)
        need(check(path,False)['unresolved_regions']==21,'initial frontier remains unresolved')
        for kind in range(4):
            altered=json.loads(json.dumps(toy))
            if kind==0:
                altered['roots'].pop()
            elif kind==1:
                altered['full_eight_point_exclusion_proved']=True
            elif kind==2:
                altered['pending'].pop()
            else:
                altered['pending'].pop(0)
                altered['nodes'][0]={'kind':'excluded','reason':'mass_log_Taylor'}
            save(altered)
            try:
                check(path,False)
            except ArithmeticError:
                rejected+=1
            else:
                raise ArithmeticError('invalid frontier accounting or unjustified Taylor exclusion accepted')
    need(rejected==4,'all four false-proof variants rejected')
    return {'initial_21_regions_preserved':True,'false_proof_variants_rejected':rejected}


if __name__=='__main__':
    print(json.dumps({'passed':True,'rotation_saturation':rotation_checks(),
        'logarithm_arithmetic':logarithm_checks(),'Taylor_arithmetic':Taylor_checks(),
        'coverage_guards':guard_checks(),'full_eight_point_exclusion_proved':False},indent=2))
