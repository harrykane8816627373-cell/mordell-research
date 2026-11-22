"""Check all phases of half-turn neighborhoods defined by relative angles."""

from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TOLERANCES = {'candidate': F(3, 5000), 'saddle': F(3, 250), 'square': F(3, 125)}
ANCHOR_NEIGHBORHOODS = {'candidate': F(11, 500), 'saddle': F(3, 125)}


def read(path):
    return json.loads(path.read_bytes())


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


class Sweep:
    def __init__(self):
        spec = importlib.util.spec_from_file_location('half_turn_phase', ROOT/'ten-point-phase-sweep/verify_all.py')
        self.phase = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.phase)
        self.box = self.phase.AngularBox()
        self.op = self.box.saddle
        self.directions = {}
        for name, points in self.box.points.items():
            self.directions[name] = tuple(self.op.scale_complex(self.op.mul(a, a), 1/self.op.norm2(a))
                                          for a in (points[i] for i in (0, 2, 3, 4, 5)))
        self.directions['square'] = tuple((self.box.exact(x), self.box.exact(y))
                                          for x, y in [(1,0), (1,0), (0,1), (-1,0), (0,-1)])

    @staticmethod
    def near_orbit(name, lo, hi):
        return (name in ANCHOR_NEIGHBORHOODS and
                2*max(abs(lo), abs(hi))+TOLERANCES[name] <= ANCHOR_NEIGHBORHOODS[name])

    def check(self, name, lo, hi, radii):
        tolerance = TOLERANCES[name]
        if self.near_orbit(name, lo, hi):
            return 'candidate_basin' if name == 'candidate' else 'saddle_orbit_excluded'
        require(radii is not None and len(radii) == 5 and sum(radii) == F(9, 2), 'reference mass')
        require(F(43, 400) <= radii[0] <= F(53, 400) and
                all(F(89, 100) <= r <= F(1271, 1000) for r in radii[1:]), 'reference radial bands')
        require(self.prod(radii[1:]) >= F(43, 40)**4, 'reference log radial domain')
        box = self.box
        op = self.op
        half = box.rotation(box.exact((hi+lo)/2))
        anchor = [self.directions[name][0],
                  *(op.mul(u, half) for u in self.directions[name][1:])]
        phases = [F(0), F(0), tolerance, tolerance, tolerance]
        grad = [box.exact(9-10*r) for r in radii]
        angle_slopes = [box.zero]*5
        quadratic = box.zero
        product = box.exact(4**5*self.prod(radii))
        for i, j in combinations(range(5), 2):
            real, imag = anchor[j]
            cosine, sine = op.mul(anchor[i], (real, -imag))
            ri, rj = radii[i], radii[j]
            scales = box.exact(2*ri*rj)
            common_width = hi-lo if i == 0 else F(0)
            edge_error = common_width+phases[i]+phases[j]
            edge_rot = (box.interval(1-edge_error**2/2, 1), box.interval(-edge_error, edge_error))
            cos_full, _ = op.mul((cosine, sine), edge_rot)
            base = box.exact(ri**2+rj**2)
            d_center = base-scales*cosine
            d_full = base-scales*cos_full
            if d_full.lo <= 0 or d_center.lo <= 0:
                return None
            term = box.exact(2*(ri**2-rj**2))/d_full
            grad[i] += term
            grad[j] -= term
            slope = 2*scales*sine/d_center
            angle_slopes[i] += slope
            angle_slopes[j] -= slope
            if (scales-base*cos_full).lo < 0:
                quadratic += scales*box.exact(edge_error**2)/d_center
            product *= d_center.square()
        phase_slope = sum(angle_slopes[1:], box.zero)
        correction = self.mag(phase_slope)*(hi-lo)+sum(self.mag(v)*e for v,e in zip(angle_slopes, phases))
        correction += box.endpoint(quadratic+2*sum((g.square() for g in grad), box.zero))
        if 0 <= correction < 2 and box.endpoint(product)*self.phase.exp_upper(correction) < box.lower:
            return 'energy_excluded'
        return None

    @staticmethod
    def prod(items):
        result = F(1)
        for item in items:
            result *= item
        return result

    def mag(self, interval):
        return F(interval.magnitude(), self.box.scale)


def main():
    manifest = read(HERE/'source_manifest.json')
    for directory, items in [(HERE, manifest['files']), (ROOT, manifest['dependencies'])]:
        for name, digest in items.items():
            require(sha256((directory/name).read_bytes()).hexdigest() == digest, 'hash '+name)
    sweep = Sweep()
    half_turn_proof = read(ROOT/'ten-point-half-turn-radial-concavity/verification_report.json')
    angle_proof = read(ROOT/'ten-point-half-turn-angle-neighborhoods/verification_report.json')
    require(half_turn_proof['passed'] and half_turn_proof['cover']['unresolved']==0 and
            F(half_turn_proof['negative_log_radial_hessian_strict_lower'])==F(1,4),
            'global half-turn log radial concavity')
    require(angle_proof['passed'] and all(
        F(angle_proof['angle_neighborhoods'][name]['squared_image_angle_tolerance'])==tol
        for name,tol in ANCHOR_NEIGHBORHOODS.items()), 'two angle neighborhoods')
    report = {}
    for name in TOLERANCES:
        cursor = F(-1)
        counts = {}
        for left, right, expected, reference in read(HERE/(name+'.json')):
            lo, hi = F(left), F(right)
            require(lo == cursor < hi <= 1, 'complete phase interval')
            result = sweep.check(name, lo, hi, list(map(F, reference)) if reference else None)
            require(result is not None and result == expected, 'phase cell accepted')
            counts[result] = counts.get(result, 0)+1
            cursor = hi
        require(cursor == 1, 'full cyclic phase')
        report[name] = {'relative_squared_angle_tolerance':str(TOLERANCES[name]), 'cells':sum(counts.values()),
                        'acceptance_counts':counts}
    print(json.dumps({'passed': True, 'normalization_mass':9, 'assumption':'half-turn symmetry',
                      'phase_covers':report, 'unresolved_inside_the_tubes':0}, indent=2))


if __name__ == '__main__':
    main()
