"""Finite exact interval cover of D2-symmetric ten-point maximizers."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('d2_basin', ROOT/'ten-point-weighted-pair-basin'/'verify_all.py')
basin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(basin)
Ball, SCALE, ZERO, exact, sqrt = basin.Ball, basin.SCALE, basin.ZERO, basin.exact, basin.sqrt

INNER = (F(41, 400), F(11, 80))
OUTER = (F(169, 200), F(164, 125))
WITNESS = 469765409


def enclosure(lo, hi):
    return Ball(exact(lo).lo, exact(hi).hi)


def intersect(v, lim):
    low, high = max(v.lo, lim.lo), min(v.hi, lim.hi)
    return Ball(low, high) if low <= high else None


def misses_zero(v):
    return v.lo > 0 or v.hi < 0


def real_powers(x):
    sq = x.square()
    return [exact(1), x, sq, sq*x, sq.square()]


def paired_powers(radius, real):
    return [exact(2), 2*real, 4*real.square()-2*radius.square(),
            8*real*real.square()-6*real*radius.square(),
            16*real.square().square()-16*real.square()*radius.square()+2*radius.square().square()]


def box_objective(a, factors):
    # Each factor is represented by its radius, real part, and degree (1 or 2).
    disc = exact(1)
    radius_product = a
    polynomials = [(a, a, 1)] + factors
    for i, (r, x, deg) in enumerate(polynomials):
        if i > 0:
            radius_product *= r if deg == 1 else r.square()
        if deg == 2:
            disc *= 4*(r.square()-x.square())
        for s, y, deg_y in polynomials[:i]:
            if deg == deg_y == 1:
                res = x-y
            elif deg == 2 and deg_y == 1:
                res = y.square()-2*x*y+r.square()
            elif deg == 1 and deg_y == 2:
                res = x.square()-2*y*x+s.square()
            else:
                res = (s.square()-r.square()).square()+4*(x-y)*(x*s.square()-y*r.square())
            disc *= res.square()
    return 4**5*radius_product*disc.square()


def dimensions(branch, box):
    a, b, h = box
    lim = enclosure(*OUTER)
    if branch == 'rectangles':
        c = intersect(exact(F(9, 4))-a/2-b, lim)
        if c is None or b.lo > c.hi:
            return None
        b = Ball(b.lo, min(b.hi, c.hi))
        c = Ball(max(c.lo, b.lo), c.hi)
        x = intersect(h, Ball(-b.hi+exact(F(1, 5)).lo, b.hi-exact(F(1, 5)).lo))
        if x is None:
            return None
        y = (-(17-10*a)*a/2-(17-10*b)*x)/(17-10*c)
        y = intersect(y, Ball(-c.hi+exact(F(1, 5)).lo, c.hi-exact(F(1, 5)).lo))
        if y is None:
            return None
        if (c-b).hi < exact(F(1, 100)).lo:
            return None
        return a, [(b, x, 2), (c, y, 2)]
    t = h
    c = intersect(exact(F(9, 2))-a-b-2*t, lim)
    if c is None or (sqrt(b)-sqrt(a)).square().hi < exact(F(2, 5)).lo:
        return None
    x = (-(17-10*a)*a-(17-10*b)*b+(17-10*c)*c)/(2*(17-10*t))
    x = intersect(x, Ball(-t.hi+exact(F(1, 5)).lo, t.hi-exact(F(1, 5)).lo))
    if x is None:
        return None
    return a, [(b, b, 1), (c, -c, 1), (t, x, 2)]


def check_force(a, factors):
    groups = [(a, real_powers(a))]
    for r, x, degree in factors:
        groups.append((r, real_powers(x) if degree == 1 else paired_powers(r, x)))
    sums = [sum((p[k] for _, p in groups), ZERO) for k in range(5)]
    for k in (2, 3, 4):
        weighed = sum((r*p[k] for r, p in groups), ZERO)
        conv = sum((sums[j]*sums[k-j] for j in range(1, k)), ZERO)
        if misses_zero((19-2*k)*sums[k]+2*conv-10*weighed):
            return False
    return True


def reciprocal_pair_force(w, radius, real):
    offset = w - real
    denominator = offset.square() + radius.square() - real.square()
    # The rectangle has r+x and r-x both at least 1/5 by separation.
    lower = max(denominator.lo, exact(F(1, 25)).lo)
    if lower > denominator.hi:
        return None
    return 2*offset/Ball(lower, denominator.hi)


def check_real_root_force(a, factors):
    singles = [(a, a)] + [(r, x) for r, x, degree in factors if degree == 1]
    for i, (r, w) in enumerate(singles):
        force = sum((exact(1)/(w-v) for j, (_, v) in enumerate(singles) if i != j), ZERO)
        for radius, real, degree in factors:
            if degree == 2:
                term = reciprocal_pair_force(w, radius, real)
                if term is None:
                    return False
                force += term
        if misses_zero(exact(F(1, 2))-5*r+2*w*force):
            return False
    return True


def candidate_neighborhood(a, factors, reps):
    coords = [(sqrt(a), ZERO)]
    for r, x, _ in factors:
        # Each factor is nondegenerate on the relevant feasible subset.
        plus, minus = (r+x)/2, (r-x)/2
        if plus.lo <= 0 or minus.lo <= 0:
            return False
        u, v = sqrt(plus), sqrt(minus)
        coords.extend([(u, v), (u, -v)])
    error = sum(((x-s).square()+(y-t).square() for (x,y), (s,t) in zip(coords,reps)), ZERO)
    return F(error.hi, SCALE) <= basin.AMPLITUDE**2


def cover(branch, reps, period=None, max_depth=100):
    lim = enclosure(*OUTER)
    extra = enclosure(-OUTER[1]+F(1,5), OUTER[1]-F(1,5)) if branch == 'rectangles' else lim
    stack = [((enclosure(*INNER), lim, extra), 0)]
    counts = Counter()
    deepest = 0
    while stack:
        box, depth = stack.pop()
        deepest = max(deepest, depth)
        counts['nodes'] += 1
        if period and counts['nodes'] % period == 0:
            print(branch, counts['nodes'], depth, flush=True)
        physical = dimensions(branch, box)
        if physical is None:
            counts['geometry'] += 1
            continue
        a, factors = physical
        if not check_real_root_force(a, factors) or not check_force(a, factors):
            counts['stationarity'] += 1
            continue
        if box_objective(a, factors).hi < exact(WITNESS).lo:
            counts['value'] += 1
            continue
        if branch == 'rectangles' and candidate_neighborhood(a, factors, reps):
            counts['in_sharp_basin'] += 1
            continue
        if depth == max_depth:
            print([tuple(F(t, SCALE) for t in (v.lo,v.hi)) for v in box], flush=True)
            raise ArithmeticError('unresolved box')
        # Relative sensitivity of the small-radius and angle-coordinate intervals.
        index = max(range(3), key=lambda i: (box[i].hi-box[i].lo)*(4 if i == 0 else 1))
        mid = (box[index].lo+box[index].hi)//2
        basin.need(box[index].lo < mid < box[index].hi, 'checked split width')
        for edge in (Ball(box[index].lo, mid), Ball(mid, box[index].hi)):
            children = box[:index] + (edge,) + box[index+1:]
            stack.append((children, depth+1))
        counts['split'] += 1
    basin.need(counts['nodes'] == 2*counts['split'] + 1, 'complete binary cover')
    return dict(sorted(counts.items()), maximum_depth=deepest)


def sources():
    files = json.loads((HERE/'source_manifest.json').read_bytes())
    for name, digest in files.items():
        basin.need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, f'fixed input {name}')

    def read(folder):
        report = json.loads((ROOT/folder/'verification_report.json').read_bytes())
        basin.need(report['passed'], f'upstream certificate {folder}')
        return report

    radii = read('ten-point-refined-radii')
    basin.need(radii['applies_to_every_global_maximizer'] and radii['normalization_mass'] == 9
               and (radii['exact_inner_count'], radii['exact_outer_count']) == (2, 8)
               and tuple(map(F, radii['outer_squared_radius_open_interval'])) == OUTER, 'global outer band')
    inner = read('ten-point-conditioned-inner')
    basin.need(inner['applies_to_every_global_maximizer'] and inner['normalization_mass'] == 9
               and tuple(F(x)/4 for x in inner['inner_squared_distance_open_interval']) == INNER, 'global inner edge')
    dispersion = read('ten-point-radial-dispersion')
    basin.need(dispersion['applies_to_every_global_maximizer'] and dispersion['normalization_mass'] == 9
               and F(dispersion['outer_squared_radius_span_strict_lower']) == F(1, 100), 'global outer dispersion')
    domain = read('maximizer-domain')
    basin.need(domain['symmetry_assumed'] is False and F(domain['pair_squared_distance_strict_lower']) == F(2, 5),
               'global pair separation')
    basin.need(F(9, 10)**2 < OUTER[0] and OUTER[1] < F(23, 20)**2 and F(1, 4)**2 < F(2, 5),
               'at most one outer point on each ray')


def check_witness():
    _, _, _, _, contacts = basin.contact.candidate()
    product = exact(1)
    for contact in contacts:
        product *= basin.contact.power(contact['d'], contact['weight'])
    basin.need(product.lo > exact(WITNESS).hi, 'certified feasible value strictly above cover threshold')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--progress', type=int)
    parser.add_argument('--branch', choices=('all','rectangles','axes'), default='all')
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    sources()
    reps = basin.algebraic_data()
    basin.matrix_bounds(reps)
    basin.remainder_bounds()
    check_witness()
    results = {}
    for branch in ('axes','rectangles') if args.branch=='all' else (args.branch,):
        results[branch] = cover(branch, reps, period=args.progress)
    report = {
        'passed': True,
        'normalization_mass': 9,
        'assumption': 'D2 symmetry of the global maximizer',
        'conditional_candidate_orbit_classified': args.branch == 'all',
        'unresolved': 0,
        'branches': results,
    }
    if args.record_report:
        basin.need(args.branch == 'all', 'both branches required for frozen report')
        (HERE/'verification_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    elif args.branch == 'all':
        basin.need(report == json.loads((HERE/'verification_report.json').read_bytes()), 'frozen conditional cover report')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
