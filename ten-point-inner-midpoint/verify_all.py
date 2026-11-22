"""Replay a first-harmonic inner-midpoint bound on the latest necessary domain."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from math import comb, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BANDS = ((F('0.084'), F('0.159')), (F('0.89'), F('1.271')))
COEFFICIENTS = tuple(map(F, ('-.00400788', '-.08575218', '.66195744',
                            '-.94181080', '-.17809273', '.18855556')))
ERRORS = (F('.00604'), F('.013535'))
SLOPE = F('-.03056048')
CAP = F(87, 1000)
W_BAND = (F(41, 400), F(11, 80))
CURVATURE = F(219, 10)
ANGULAR_CURVATURE = F(161, 50)
GRAM_BUDGET = F(10981, 2500000)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require(test, detail):
    if not test:
        raise ArithmeticError(detail)


def read(path):
    return json.loads(path.read_bytes())


def scalar(entry):
    return F(int(entry['num']), int(entry['den']))


def inputs(repo):
    manifest = read(HERE/'source_manifest.json')
    for group, base in (('files', HERE), ('dependencies', repo)):
        for name, digest in manifest[group].items():
            require(sha256((base/name).read_bytes()).hexdigest() == digest, 'fixed input '+name)
    radii = read(repo/'ten-point-angular-radii/verification_report.json')
    inner = read(repo/'ten-point-conditioned-inner/verification_report.json')
    for report in radii, inner:
        require(report['passed'] and report['applies_to_every_global_maximizer']
                and report['normalization_mass']==9 and report['unresolved']==0, 'accepted necessary range')
    require((radii['exact_inner_count'], radii['exact_outer_count']) == (2, 8)
            and tuple(tuple(map(F, radii[k])) for k in ('inner_squared_radius_open_interval',
                                                       'outer_squared_radius_open_interval')) == BANDS,
            'current radial bands and population')
    require(tuple(F(x)/4 for x in inner['inner_squared_distance_open_interval']) == W_BAND,
            'current inner half-difference band')
    require(F(inner['inner_angle_cosine_strict_upper']) == -F(77,100), 'accepted negative inner-angle cosine')
    curved = load(repo/'ten-point-curved-kernel/verify_all.py', 'new_midpoint_kernel')
    upper = curved.upper
    kernel = read(curved.HERE/'kernel.json')
    report = read(repo/'ten-point-angular-upper/verification_report.json')
    require(report['passed'] and report['normalization_mass']==9 and report['unresolved']==0,
            'accepted angular-class deficits')
    classes = [(x['class'], x['weight']) for x in report['classes']]
    require(classes == [([0,0,0],1),([0,1,0],16),([1,1,1],8),([1,1,2],8),([1,1,3],8),([1,1,4],4)],
            'every unordered deficiency counted once')
    allowance = -sum(x['weight']*F(x['deficit_floor']) for x in report['classes'])
    require(2*(upper.hx(kernel['B']) - upper.log.log_lower(upper.actual_root_lower()) + allowance) < GRAM_BUDGET,
            'strict total Gram budget')
    return upper, kernel


def gram_norm(upper, kernel):
    E = [[upper.hx(x) for x in row] for row in kernel['transforms'][0]]
    require(E == [[F(i==j) for j in range(6)] for i in range(6)], 'original monomial first-harmonic basis')
    Q = [[upper.hx(x) for x in row[2:]] for row in kernel['Q_blocks'][1][2:]]
    upper.gram.positive_ldl(Q)
    native = read(HERE/'jacobian_solve.json')
    proof = native['output']
    c = COEFFICIENTS[2:]
    require(native['operation_id']=='matrix.rational_linear_system.solve' and proof['outcome']=='UNIQUE'
            and [[scalar(x) for x in row] for row in proof['matrix']['entries']]==Q
            and tuple(map(scalar, proof['rhs']))==c, 'exact Jacobian Gram-solve source')
    solution = tuple(map(scalar, proof['solution']))
    require(len(solution)==len(c) and tuple(sum(a*b for a,b in zip(row, solution)) for row in Q)==c,
            'exact rational solution substitution')
    return sum(a*b for a,b in zip(c, solution))


def radius_errors(repo):
    radial = load(repo/'ten-point-ward-inner-count/verify_radial.py', 'new_midpoint_radial')
    g = [sum(COEFFICIENTS[j]*comb(j,k)*(-1)**(j-k) for j in range(k,6)) for k in range(6)]
    covers = []
    for band_index, (band,error) in enumerate(zip(BANDS,ERRORS)):
        h = g.copy()
        h[0] -= band_index==0
        intercept = error-F(9,10)*SLOPE
        require(all(intercept+SLOPE*r > 0 for r in band), 'positive affine error on full band')
        lhs = radial.mul([intercept, SLOPE], [intercept, SLOPE])
        difference = radial.mul([F(0), F(-1)], radial.mul(h,h))
        for index,x in enumerate(lhs):
            difference[index] += x
        covers.append(radial.strictly_above(difference,band,F(0)))
    return covers


def inner_curvature():
    # For F(w) with fixed outer points, axis and midpoint, F'' <= -1/W^2 + 2 sum g(u)^-2.
    a,b = BANDS[1]
    W = W_BAND[1]
    require(F(9,10)**2<a and CAP<F(9,100) and b+CAP**2<F(13,10), 'positive g derivative and upper envelope')
    require(CAP/(2*F(9,10)**3)<F(1,10) and 3*F(9,10)**2>F(13,100), 'g^-2 is strictly decreasing and convex')

    def inverse_cap(u):
        scale = 10**7
        r = F(isqrt(u.numerator*scale**2//u.denominator),scale)
        g = (r-CAP)**2-W
        require(r>CAP and g>0, 'positive denominator for every inner-outer pair product')
        return 1/g**2

    fa,fb = map(inverse_cap,(a,b))
    outer_mass = 9-2*(W+CAP**2)
    require(fa>fb and 8*a<outer_mass<8*b, 'decreasing chord applied to total outside mass')
    sum_upper = 8*fa+(fb-fa)*(outer_mass-8*a)/(b-a)
    require(1/W**2-2*sum_upper > CURVATURE, 'strict inner squared half-difference curvature')


def angular_curvature(repo):
    gaps = read(repo/'ten-point-outer-gap/verification_report.json')
    require(gaps['passed'] and gaps['applies_to_every_global_maximizer'] and gaps['normalization_mass']==9
            and gaps['unresolved']==0, 'accepted successive-gap convex angular domain')
    lo,hi = map(F,gaps['outer_successive_angle_cosine_open_interval'])
    require((lo,hi)==(F(59,100),F(39,50)), 'fixed successive cosine band')

    def K(q,c):
        return (1-q*c)/(q-c)**2

    (ulo,uhi),(vlo,vhi) = BANDS
    q = F(127,125)
    require((vlo+vhi)**2<4*q*q*vlo*vhi and 2-q*q-q*hi>0 and hi*q+hi**2<2,
            'outer radial ratio and monotonicity in q and c')
    base,w1,w2 = 1/(q+1),F(11,5),F(3,4)
    require(lo>F(1,2) and K(q,lo)>w1 and K(q,2*lo*lo-1)>w2,
            'outer cycle and two-step curvature bounds')
    require(F(283,200)**2>2, 'first cycle eigenvalue exceeds 117/200')
    outer = 8*base+min((w1-base)*F(117,200)+2*(w2-base), 2*(w1-base))

    s = F(423,1000)
    require(F(1,4)**2<ulo/vhi and uhi/vlo<s*s<1, 'inner-outer radial ratio')
    q0,q1 = (1+s*s)/(2*s),F(17,8)
    inner_cos = F(77,100)
    far = F(17,50)
    require((1-inner_cos)/2<far**2 and q1*far<1, 'one inner direction contributes nonpositive curvature')
    t = F(49,50)
    require(t*q0+t*t>2 and -K(q0,t)<F(9,4), 'combined curvature bound away from the inner ray')
    c = F(31,50)
    require(inner_cos*t>c and (inner_cos*t-c)**2>(1-inner_cos**2)*(1-t*t),
            'opposite inner edge cosine below -31/50 near the inner ray')
    require(2-q1*q1+q1<0 and K(q1,-c)>F(3,10), 'negative curvature from the opposite inner edge')
    require(2*s/(1-s)**2<F(51,20) and outer-F(9,4)>ANGULAR_CURVATURE,
            'strengthened mean-zero angular curvature after subtracting both inner edges')


def verify(repo):
    upper,kernel = inputs(repo)
    dual = gram_norm(upper,kernel)
    covers = radius_errors(repo)
    affine = 2*ERRORS[0]+8*ERRORS[1]
    spare = 2*CAP-affine
    require(0<dual and spare>0 and spare**2>GRAM_BUDGET*dual, 'bound half the inner point sum')
    inner_curvature()
    angular_curvature(repo)
    return {'passed':True,'normalization_mass':9,'applies_to_every_global_maximizer':True,
            'inner_midpoint_norm_strict_upper':str(CAP),'inner_midpoint_squared_strict_upper':str(CAP**2),
            'outer_centroid_norm_strict_upper':str(CAP/4),
            'inner_squared_half_difference_negative_hessian_lower':str(CURVATURE),
            'outer_mean_zero_negative_angle_hessian_lower_on_successive_gap_polytope':str(ANGULAR_CURVATURE),
            'first_harmonic_gram_budget_strict_upper':str(GRAM_BUDGET),
            'indicator_error_covers':covers, 'unresolved':0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=ROOT)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    report = verify(args.repo.resolve())
    path = HERE/'verification_report.json'
    out = json.dumps(report,indent=2)+'\n'
    if args.record_report:
        path.write_text(out,encoding='ascii',newline='\n')
    else:
        require(report==read(path), 'saved midpoint report')
    print(out,end='')


if __name__=='__main__':
    main()
