"""Integer acceptance of mass-adjusted logarithmic Taylor exclusion."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'joint-closure'))
import interval_engine as E

BASE_SHA256 = '9b40704b8be9b15a219c89b64989bd71caddf74564386570900164494f857b58'
REFERENCES = json.loads((Path(__file__).resolve().parent/'reference_centers.json').read_text())


def _log_reduced(a, b):
    """Enclose log(a/b) for 1 <= a/b <= 2 by a positive atanh series."""
    if not (0 < b <= a <= 2*b):
        raise ValueError('reduced positive ratio')
    if a == b:
        return E.ZERO
    x = E.rat(a-b, a+b)
    x2, power, result = E.mul(x, x), x, E.ZERO
    terms = 40
    for j in range(terms):
        result = E.add(result, E.divide(E.scale(power, 2), 2*j+1))
        power = E.mul(power, x2)
    # x <= 1/3 exactly. Positive tail <= 9/[4(2m+1)3^(2m+1)].
    tail = E.rat(9, 4*(2*terms+1)*3**(2*terms+1))[1]
    return result[0], result[1]+tail


def log_ratio(a, b):
    """Integer fixed-point interval for log(a/b), for arbitrary positive integers."""
    if a <= 0 or b <= 0:
        raise ValueError('positive logarithm arguments')
    exponent = a.bit_length()-b.bit_length()
    num, den = (a, b << exponent) if exponent >= 0 else (a << -exponent, b)
    if num < den:
        exponent -= 1
        num *= 2
    if not den <= num < 2*den:
        raise ArithmeticError('binary logarithm reduction')
    return E.add(_log_reduced(num, den), E.scale(_log_reduced(2, 1), exponent))


def target_log_upper(box, qbox, radius, midpoint=None):
    """Upper log(Delta/7^7) for every mass-seven actual point in these slabs.

    None means the Taylor disks were too wide; it never means exclusion.
    """
    if midpoint is None:
        midpoint, _ = E.slab_center(box, radius, qbox[7])
    if len(midpoint) != 13 or not all(type(v) is int for v in midpoint):
        raise ValueError('exact rational center in 13 centered-coordinate variables')
    center = E.positions([E.point(x) for x in midpoint])
    raw_norms = [z[0][0]**2+z[1][0]**2 for z in center]
    product = 1
    remainder_upper = 0
    for i, j in E.PAIRS:
        dx = center[i][0][0]-center[j][0][0]
        dy = center[i][1][0]-center[j][1][0]
        raw_distance = dx*dx+dy*dy
        if not raw_distance:
            return None
        product *= raw_distance
        distance_lower = E.sqrt_lo(raw_distance//E.S)
        if distance_lower <= 0:
            return None
        displacement = E.cs(E.cs(qbox[i], qbox[j]), E.cs(center[i], center[j]))
        displacement_upper = E.sqrt_hi(E.cnorm(displacement)[1])
        rho = -(-displacement_upper*E.S//distance_lower)
        if rho >= E.S:
            return None
        remainder_upper += E.mul(E.square((0, rho)),
            E.reciprocal((E.S-rho, E.S-rho)))[1]
    inverse, _ = E.invert_edges(center)
    gradient, _ = E.force_gradient(center, inverse)
    gradient = E.reduced_gradient(gradient)
    linear_upper = 0
    for g, interval, c in zip(gradient, box, midpoint):
        linear_upper += E.scale(E.mul(E.neg(g), E.sub(interval, E.point(c))), 2)[1]
    mass_correction = E.scale(E.sub(E.rat(7), E.rat(sum(raw_norms), E.S**2)), 4)
    center_log = log_ratio(product, 7**7*E.S**56)
    upper = center_log[1]+mass_correction[1]+linear_upper+remainder_upper
    return {'log_Delta_over_target_upper': upper,
        'center_log_upper': center_log[1], 'mass_correction_upper': mass_correction[1],
        'linear_residual_upper': linear_upper, 'Taylor_remainder_upper': remainder_upper}


def contract(box, settings, passes=5):
    box, reason, aux = E.contract(box, settings, passes=passes)
    if box is None:
        return box, reason, aux
    upper = target_log_upper(box, aux[0], settings['radius'])
    if upper is not None and upper['log_Delta_over_target_upper'] < 0:
        upper['center_source'] = 'midpoint'
        return None, 'mass_log_Taylor', upper
    if REFERENCES['precision_bits'] != E.BITS:
        raise ArithmeticError('reference grid matches interval grid')
    midpoint = [(a+b)//2 for a, b in box]
    nearest = sorted(enumerate(REFERENCES['centers']),
        key=lambda item: sum((a-b)**2 for a, b in zip(midpoint, item[1])))[:2]
    for index, center in nearest:
        upper = target_log_upper(box, aux[0], settings['radius'], center)
        if upper is not None and upper['log_Delta_over_target_upper'] < 0:
            upper['center_source'] = 'reference_'+str(index)
            return None, 'mass_log_Taylor', upper
    return box, None, aux
