"""Bound the outer Cauchy field on inner disks by a circular dynamic program."""

from collections import deque
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
N = 512
RADIAL_CELLS = 24
SUPPORTS = 32
OUTPUT_SCALE = 10**6


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dyadic = load('inner_field_dyadic', 'validated-maxima/dyadic.py')
symmetry = load('inner_field_symmetry', 'candidate-global-screen/verify_symmetry.py')
Ball, SCALE = dyadic.Ball, dyadic.SCALE


def rational(num, den):
    return Ball.integer(num)/den


def sqrt(x):
    return symmetry.sqrt_ball(x, Ball, SCALE)


def circle():
    """Unit circle at uniform angles, starting with the quarter-turn root."""
    c, s = Ball.integer(0), Ball.integer(1)
    order = 4
    while order < N:
        c, s = sqrt((1+c)/2), sqrt((1-c)/2)
        order *= 2
    pts = [(Ball.integer(1), Ball.integer(0))]
    for _ in range(N-1):
        x, y = pts[-1]
        pts.append((c*x-s*y, s*x+c*y))
    return pts


def one_site_bounds(pts, radius):
    low, high = rational(943, 1000), rational(1128, 1000)
    require(low.square().hi < rational(89, 100).lo and
            high.square().lo > rational(1271, 1000).hi, 'outer radial enclosure')
    width = (high-low)/RADIAL_CELLS
    # Half the radial cell width plus the chord to the central direction.
    delta = width/2 + high*rational(22, 7*N)
    rows = []
    for c, s in pts:
        cells = []
        for j in range(RADIAL_CELLS):
            r = low+(2*j+1)*width/2
            x, y = r*c-radius, r*s
            d = x.square()+y.square()
            gap = sqrt(d)-delta
            require(gap.lo > 0, 'pole outside disk about central sample')
            denom = d.square()
            real = (x.square()-y.square())/denom
            imag = (-2*x*y)/denom
            error = 2*delta/(gap*gap.square())
            cells.append((real, imag, error))
        rows.append(cells)
    return rows


def cyclic_support(weights, minimum, maximum):
    """Maximize eight weights over lifted indices with the closing gap enforced."""
    best = None
    for first in range(maximum):
        previous = {first: weights[first]}
        for step in range(1, 8):
            left = first + max(step*minimum, N-(8-step)*maximum)
            right = first + min(step*maximum, N-(8-step)*minimum)
            window = deque()
            current = {}
            for prev in range(left-maximum, left-minimum):
                if prev in previous:
                    while window and previous[window[-1]] <= previous[prev]:
                        window.pop()
                    window.append(prev)
            for pos in range(left, right+1):
                prev = pos-minimum
                if prev in previous:
                    while window and previous[window[-1]] <= previous[prev]:
                        window.pop()
                    window.append(prev)
                while window and window[0] < pos-maximum:
                    window.popleft()
                if window:
                    current[pos] = weights[pos % N]+previous[window[0]]
            previous = current
        require(previous, 'nonempty set of closed eight-point walks')
        value = max(previous.values())
        best = value if best is None else max(value, best)
    return best


def certify(pts, numerator, denominator, bound):
    rows = one_site_bounds(pts, rational(numerator, denominator))
    minimum, maximum = 55, 77
    require(pts[minimum][0].lo > rational(39, 50).hi and
            pts[maximum][0].hi < rational(59, 100).lo, 'rounded gap limits')
    results = []
    for index in range(SUPPORTS//2+1):
        c, s = pts[index*N//SUPPORTS]
        weights = []
        for cells in rows:
            raw = max((c*x+s*y+err).hi for x, y, err in cells)
            weights.append(-((-raw*OUTPUT_SCALE)//SCALE))
        results.append(cyclic_support(weights, minimum, maximum))
    cosine = pts[N//(2*SUPPORTS)][0]
    # The other half of the support grid follows by reflecting the shell.
    upper = Fraction(max(results)*SCALE, OUTPUT_SCALE*cosine.lo)
    require(upper < bound, 'complex modulus threshold')
    return {'disk_radius': str(Fraction(numerator, denominator)),
            'projected_upper_numerators': results, 'denominator': OUTPUT_SCALE,
            'field_modulus_strict_upper': str(bound)}


def main():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for directory, items in ((HERE, manifest['files']), (ROOT, manifest['dependencies'])):
        for name, digest in items.items():
            require(sha256((directory/name).read_bytes()).hexdigest() == digest, 'hash '+name)
    pts = circle()
    small = certify(pts, 2, 5, Fraction(29, 10))
    large = certify(pts, 23, 50, Fraction(18, 5))
    require(Fraction(159, 1000) < Fraction(2, 5)**2 and
            (Fraction(23, 50)-Fraction(87, 1000))**2 > Fraction(53, 400),
            'inner disks contain the admissible points')
    print(json.dumps({'passed': True, 'angular_cells': N, 'radial_cells': RADIAL_CELLS,
                      'support_directions': SUPPORTS, 'disks': [small, large]}, indent=2))


if __name__ == '__main__':
    main()
