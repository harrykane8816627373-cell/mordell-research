"""Integer outward intervals for actual eight-point coordinates; no float acceptance."""
from math import isqrt
from itertools import combinations

BITS = 40
S = 1 << BITS
ZERO = (0, 0)
ONE = (S, S)
PAIRS = [(i, j) for i in range(8) for j in range(i+1, 8)]
PAIR_INDEX = {v: i for i, v in enumerate(PAIRS)}
SUBSETS = [(m, [PAIR_INDEX[p] for p in combinations(ix, 2)])
           for m in range(3, 8) for ix in combinations(range(8), m)]
# q=(x0, x1+iy1,...,x6+iy6, -sum(q0..q6)); y0=0.
L = [[(0, 0)]*13 for _ in range(8)]
L[0][0] = (1, 0)
L[7][0] = (-1, 0)
for k in range(1, 7):
    L[k][2*k-1] = (1, 0)
    L[k][2*k] = (0, 1)
    L[7][2*k-1] = (-1, 0)
    L[7][2*k] = (0, -1)
EDGE_L = [[(L[i][k][0]-L[j][k][0], L[i][k][1]-L[j][k][1]) for k in range(13)] for i, j in PAIRS]
GRAM = [[sum(L[i][a][0]*L[i][b][0]+L[i][a][1]*L[i][b][1] for i in range(8))
         for b in range(13)] for a in range(13)]


def setup(bits):
    global BITS, S, ZERO, ONE
    BITS = bits
    S = 1 << bits
    ZERO = (0, 0)
    ONE = (S, S)


def rat(n, d=1):
    return n*S//d, -(-n*S//d)


def point(n):
    return n, n


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    values = [a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]]
    return min(values)//S, -(-max(values)//S)


def scale(a, n):
    return (a[0]*n, a[1]*n) if n >= 0 else (a[1]*n, a[0]*n)


def divide(a, n):
    if n <= 0:
        raise ValueError('positive integer divisor')
    return a[0]//n, -(-a[1]//n)


def square(a):
    low = 0 if a[0] <= 0 <= a[1] else min(a[0]*a[0], a[1]*a[1])
    high = max(a[0]*a[0], a[1]*a[1])
    return low//S, -(-high//S)


def reciprocal(a):
    if a[0] <= 0:
        raise ZeroDivisionError('positive denominator required')
    return S*S//a[1], -(-S*S//a[0])


def meet(a, b):
    out = max(a[0], b[0]), min(a[1], b[1])
    if out[0] > out[1]:
        raise ValueError('empty contraction')
    return out


def sqrt_hi(value):
    if value < 0:
        raise ValueError('negative squared bound')
    v = value*S
    k = isqrt(v)
    return k+(k*k != v)


def sqrt_lo(value):
    return isqrt(max(0, value)*S)


def ca(a, b):
    return add(a[0], b[0]), add(a[1], b[1])


def cn(a):
    return neg(a[0]), neg(a[1])


def cs(a, b):
    return ca(a, cn(b))


def cm(a, b):
    return sub(mul(a[0], b[0]), mul(a[1], b[1])), add(mul(a[0], b[1]), mul(a[1], b[0]))


def cc(a):
    return a[0], neg(a[1])


def ci(n):
    return rat(n), ZERO


def cscale(a, n):
    return scale(a[0], n), scale(a[1], n)


def cnorm(a):
    return add(square(a[0]), square(a[1]))


def cpow(a, n):
    out = ci(1)
    for _ in range(n):
        out = cm(out, a)
    return out


def positions(box):
    q = [(box[0], ZERO)]+[(box[2*k-1], box[2*k]) for k in range(1, 7)]
    q.append(cn((tuple(sum(v[0][side] for v in q) for side in [0, 1]),
                 tuple(sum(v[1][side] for v in q) for side in [0, 1]))))
    return q


def variables(q):
    return [q[0][0]]+[q[k][component] for k in range(1, 7) for component in [0, 1]]


def model(scope='global_maximum'):
    if scope not in ['global_maximum', 'user_physical_system']:
        raise ValueError('known scope')
    radius2 = rat(49, 32)[1] if scope == 'global_maximum' else rat(49, 8)[1]
    radius = sqrt_hi(radius2)
    lower = sqrt_lo(rat(7, 8)[0])
    box = [(lower, radius)]+[(-radius, radius) for _ in range(12)]
    box[2] = (0, radius)  # reflection: y1>=0
    return {'scope': scope, 'radius2': radius2, 'radius': radius,
            'propagation_version': 2,
            'separation': rat(1, 2)[0] if scope == 'global_maximum' else rat(2, 7)[0],
            'initial_box': box}


def close_center(q):
    # Linear equation sum q_i=0, propagated to every component.
    for component in [0, 1]:
        for i in range(8):
            lo = -sum(q[j][component][1] for j in range(8) if j != i)
            hi = -sum(q[j][component][0] for j in range(8) if j != i)
            new = meet(q[i][component], (lo, hi))
            q[i] = (new, q[i][1]) if component == 0 else (q[i][0], new)


def invert_edges(q, separation=None, ordered=False):
    result, distances = [], []
    for i, j in PAIRS:
        z = cs(q[i], q[j])
        if ordered:
            z = (meet(z[0], (0, z[0][1])) if i == 0 else
                 meet(z[0], (z[0][0], 0))), z[1]
        distance = cnorm(z)
        distances.append(distance)
        lo = max(distance[0], separation or 0)
        invden = reciprocal((lo, distance[1]))
        result.append((mul(z[0], invden), neg(mul(z[1], invden))))
    return result, distances


def force_gradient(q, inverse):
    force = [ci(0) for _ in range(8)]
    for (i, j), value in zip(PAIRS, inverse):
        force[i] = ca(force[i], value)
        force[j] = cs(force[j], value)
    return [cs(cscale(q[i], 4), cc(force[i])) for i in range(8)], force


def reduced_gradient(g):
    return [sub(g[0][0], g[7][0])]+[
        sub(g[k][component], g[7][component]) for k in range(1, 7) for component in [0, 1]]


def force_jacobian(q):
    inv, _ = invert_edges(q)
    matrix = [[rat(4*GRAM[a][b]) for b in range(13)] for a in range(13)]
    for values, zinv in zip(EDGE_L, inv):
        z2 = cm(zinv, zinv)
        active = [k for k, value in enumerate(values) if value != (0, 0)]
        for a in active:
            for b in active:
                ar, ai = values[a]
                br, bi = values[b]
                coeff = ar*br-ai*bi, ar*bi+ai*br
                term = sub(scale(z2[0], coeff[0]), scale(z2[1], coeff[1]))
                matrix[a][b] = add(matrix[a][b], term)
    return matrix


def reduced_hessian(box):
    return force_jacobian(positions(box))


def power_columns(q):
    powers = [[ci(1)] for _ in range(8)]
    for i in range(8):
        for k in range(1, 9):
            powers[i].append(cm(powers[i][-1], q[i]))
    sums = [ci(8), ci(0)]
    for k in range(2, 9):
        total = ci(0)
        for i in range(8):
            total = ca(total, powers[i][k])
        sums.append(total)
    return powers, sums


def polynomial_columns(q):
    powers, sums = power_columns(q)
    action = [[ci(0), cscale(cc(q[i]), 8)] for i in range(8)]
    for i in range(8):
        t = cnorm(q[i])
        for k in range(2, 8):
            value = (mul(scale(t, 8*k), powers[i][k-2][0]),
                     mul(scale(t, 8*k), powers[i][k-2][1]))
            for r in range(k-1):
                value = cs(value, cscale(cm(sums[r], powers[i][k-2-r]), 2*(k-1-r)))
            value = ca(value, cscale(powers[i][k-2], k*(k-1)))
            action[i].append(value)
    return powers, action, sums


def polynomial_rayleigh(q, coefficients):
    powers, action, _ = polynomial_columns(q)
    value = ZERO
    for i in range(8):
        v, w = ci(0), ci(0)
        for k, (a, b) in enumerate(coefficients):
            vp = (sub(scale(powers[i][k][0], a), scale(powers[i][k][1], b)),
                  add(scale(powers[i][k][1], a), scale(powers[i][k][0], b)))
            wp = (sub(scale(action[i][k][0], a), scale(action[i][k][1], b)),
                  add(scale(action[i][k][1], a), scale(action[i][k][0], b)))
            v, w = ca(v, vp), ca(w, wp)
        value = add(value, sub(scale(cnorm(v), 64), cnorm(w)))
    return value


def polynomial_equations(q):
    powers, sums = power_columns(q)
    norms = [cnorm(z) for z in q]
    out = [sub((sum(v[0] for v in norms), sum(v[1] for v in norms)), rat(7))]
    for k in range(1, 7):
        mixed = ci(0)
        for i in range(8):
            mixed = ca(mixed, (mul(norms[i], powers[i][k][0]), mul(norms[i], powers[i][k][1])))
        residual = cs(cscale(mixed, 8), cscale(sums[k], 15-k))
        for ell in range(1, k):
            residual = cs(residual, cm(sums[ell], sums[k-ell]))
        out.extend(residual)
    return out


def polynomial_jacobian(q):
    powers, sums = power_columns(q)
    norms = [cnorm(z) for z in q]
    matrix = [[ZERO for _ in range(13)] for _ in range(13)]
    for column in range(13):
        i = 0 if column == 0 else (column+1)//2
        imaginary = column != 0 and column % 2 == 0

        def rotate(value):
            return (neg(value[1]), value[0]) if imaginary else value

        mass_difference = cs(q[i], q[7])
        matrix[0][column] = scale(mass_difference[1] if imaginary else mass_difference[0], 2)
        ds = [ci(0)]
        for k in range(1, 7):
            ds.append(rotate(cscale(cs(powers[i][k-1], powers[7][k-1]), k)))
        for k in range(1, 7):
            a_i = (mul(scale(norms[i], k+1), powers[i][k-1][0]),
                   mul(scale(norms[i], k+1), powers[i][k-1][1]))
            a_7 = (mul(scale(norms[7], k+1), powers[7][k-1][0]),
                   mul(scale(norms[7], k+1), powers[7][k-1][1]))
            if imaginary:
                dQ = rotate(cs(cs(a_i, powers[i][k+1]), cs(a_7, powers[7][k+1])))
            else:
                dQ = cs(ca(a_i, powers[i][k+1]), ca(a_7, powers[7][k+1]))
            derivative = cs(cscale(dQ, 8), cscale(ds[k], 15-k))
            for ell in range(1, k):
                derivative = cs(derivative, ca(cm(ds[ell], sums[k-ell]),
                                               cm(sums[ell], ds[k-ell])))
            matrix[2*k-1][column], matrix[2*k][column] = derivative
    return matrix


def slab_center(box, radius, extra_slab=None):
    """A rational center in the convex centroid/radius-coordinate slabs."""
    qbox = positions(box)
    qbox[7] = meet(qbox[7][0], (-radius, radius)), meet(qbox[7][1], (-radius, radius))
    if extra_slab is not None:
        qbox[7] = meet(qbox[7][0], extra_slab[0]), meet(qbox[7][1], extra_slab[1])
    midpoint = [(a+b)//2 for a, b in box]
    for component in [0, 1]:
        columns = ([0]+list(range(1, 13, 2))) if component == 0 else list(range(2, 13, 2))
        target_q7 = sum(qbox[7][component])//2
        wanted = -target_q7
        difference = wanted-sum(midpoint[j] for j in columns)
        for j in columns:
            move = min(difference, box[j][1]-midpoint[j]) if difference >= 0 else max(difference, box[j][0]-midpoint[j])
            midpoint[j] += move
            difference -= move
        if difference:
            raise ValueError('empty centroid slab')
    qm = positions([point(v) for v in midpoint])
    if not all(qbox[i][c][0] <= qm[i][c][0] <= qbox[i][c][1] for i in range(8) for c in [0, 1]):
        raise ArithmeticError('center not in convex enclosure')
    return midpoint, qbox


def newton_contract(box, preconditioner, radius=None, kind='polynomial', extra_slab=None):
    if radius is None:
        raise ValueError('explicit convex centroid slab required')
    midpoint, qbox = slab_center(box, radius, extra_slab)
    qm = positions([point(v) for v in midpoint])
    if kind == 'polynomial':
        matrix, grad = polynomial_jacobian(qbox), polynomial_equations(qm)
    elif kind == 'force':
        matrix = force_jacobian(qbox)
        inverse, _ = invert_edges(qm)
        grad, _ = force_gradient(qm, inverse)
        grad = reduced_gradient(grad)
    else:
        raise ValueError('known Newton system')
    result, contraction = [], 0
    for i in range(13):
        correction = ZERO
        for j in range(13):
            correction = add(correction, mul(point(preconditioner[i][j]), grad[j]))
        value = sub(point(midpoint[i]), correction)
        row_bound = 0
        for k in range(13):
            coefficient = rat(int(i == k))
            for j in range(13):
                coefficient = sub(coefficient, mul(point(preconditioner[i][j]), matrix[j][k]))
            row_bound += max(abs(coefficient[0]), abs(coefficient[1]))
            value = add(value, mul(coefficient, sub(box[k], point(midpoint[k]))))
        contraction = max(contraction, row_bound)
        result.append(meet(box[i], value))
    return result, contraction


def rayleigh(q, inverse, direction):
    # Integer Gaussian displacement; a positive common denominator is irrelevant to sign.
    value = rat(4*sum(a*a+b*b for a, b in direction))
    for (i, j), zinv in zip(PAIRS, inverse):
        dr, di = direction[i][0]-direction[j][0], direction[i][1]-direction[j][1]
        z2 = cm(zinv, zinv)
        term = sub(scale(z2[0], dr*dr-di*di), scale(z2[1], 2*dr*di))
        value = add(value, term)
    return value


def contract(box, settings, passes=5, direction=None):
    """All intersections preserve every hypothetical feasible actual configuration."""
    box = [tuple(v) for v in box]
    sep, radius = settings['separation'], settings['radius']
    tag = 'geometry'
    try:
        for _ in range(passes):
            old = list(box)
            q = positions(box)
            for i in range(8):
                q[i] = meet(q[i][0], (-radius, radius)), meet(q[i][1], (-radius, radius))
            tag = 'centering'
            close_center(q)
            tag = 'ordering'
            order = [1, 2, 3, 4, 5, 6, 7, 0]
            for before, after in zip(order, order[1:]):
                lo = max(q[after][0][0], q[before][0][0])
                q[after] = (meet(q[after][0], (lo, q[after][0][1])), q[after][1])
            for before, after in reversed(list(zip(order, order[1:]))):
                hi = min(q[before][0][1], q[after][0][1])
                q[before] = (meet(q[before][0], (q[before][0][0], hi)), q[before][1])
            for rank, i in enumerate(range(1, 8), 1):
                upper = -(-((rank-1)*radius-q[0][0][0])//(8-rank))
                lower = (-q[0][0][1]-(7-rank)*radius)//rank
                q[i] = meet(q[i][0], (lower, upper)), q[i][1]
            tag = 'ordered_separation'
            for left_rank, i in enumerate(order):
                for j in order[left_rank+1:]:
                    vertical = square(sub(q[j][1], q[i][1]))[1]
                    if vertical < sep:
                        gap = sqrt_lo(sep-vertical)
                        q[i] = meet(q[i][0], (q[i][0][0], q[j][0][1]-gap)), q[i][1]
                        q[j] = meet(q[j][0], (q[i][0][0]+gap, q[j][0][1])), q[j][1]
            tag = 'ellipse_covariance'
            # mass=7 and |s2|<1 imply 3<sum x_i^2<4, 3<sum y_i^2<4.
            for component in [0, 1]:
                squares = [square(z[component]) for z in q]
                low, high = sum(v[0] for v in squares), sum(v[1] for v in squares)
                if low > 4*S or high < 3*S:
                    return None, tag, None
                for i in range(8):
                    upper = 4*S-low+squares[i][0]
                    bound = sqrt_hi(upper)
                    value = meet(q[i][component], (-bound, bound))
                    lower = max(0, 3*S-high+squares[i][1])
                    minimum = sqrt_lo(lower)
                    if value[0] >= 0:
                        value = meet(value, (minimum, value[1]))
                    elif value[1] <= 0:
                        value = meet(value, (value[0], -minimum))
                    q[i] = (value, q[i][1]) if component == 0 else (q[i][0], value)
            tag = 'mass_or_radius'
            norms = [cnorm(z) for z in q]
            lo_sum, hi_sum = sum(t[0] for t in norms), sum(t[1] for t in norms)
            if lo_sum > 7*S or hi_sum < 7*S:
                return None, 'mass', None
            for i in range(8):
                lower = max(0, 7*S-hi_sum+norms[i][1])
                upper = min(settings['radius2'], 7*S-lo_sum+norms[i][0], norms[0][1])
                if i == 0:
                    lower = max(lower, rat(7, 8)[0])
                if lower > upper:
                    return None, 'maximum_radius', None
                for component in [0, 1]:
                    other = square(q[i][1-component])
                    bound = sqrt_hi(upper-other[0])
                    value = meet(q[i][component], (-bound, bound))
                    if lower > other[1]:
                        minimum = sqrt_lo(lower-other[1])
                        if value[0] >= 0:
                            value = meet(value, (minimum, value[1]))
                        elif value[1] <= 0:
                            value = meet(value, (value[0], -minimum))
                    q[i] = (value, q[i][1]) if component == 0 else (q[i][0], value)
            distances = [cnorm(cs(q[i], q[j])) for i, j in PAIRS]
            if any(d[1] <= sep for d in distances):
                return None, 'separation', None
            if settings['scope'] == 'global_maximum':
                for m, edges in SUBSETS:
                    if 8*sum(distances[k][1] for k in edges) <= m*m*(m-1)*S:
                        return None, 'strict_subset', None
            product = 1
            for d in distances:
                product *= d[1]
            if product < 7**7*S**28:
                return None, 'target_product', None
            s2, s3 = ci(0), ci(0)
            for z in q:
                z2 = cm(z, z)
                s2 = ca(s2, z2)
                s3 = ca(s3, cm(z2, z))
            if 175*cnorm(s2)[0]+54*cnorm(s3)[0] >= 175*S:
                return None, 'old_ellipse', None
            deficit = ZERO
            for z in q:
                t = cnorm(z)
                deficit = add(deficit, mul(t, square(sub(t, ONE))))
            if 1200*deficit[1] <= S:
                return None, 'radial_gap', None
            powers, sums = power_columns(q)
            for k in range(1, 7):
                mixed = ci(0)
                for i in range(8):
                    mixed = ca(mixed, (mul(cnorm(q[i]), powers[i][k][0]),
                                       mul(cnorm(q[i]), powers[i][k][1])))
                residual = cs(cscale(mixed, 8), cscale(sums[k], 15-k))
                for ell in range(1, k):
                    residual = cs(residual, cm(sums[ell], sums[k-ell]))
                if any(v[0] > 0 or v[1] < 0 for v in residual):
                    return None, 'polynomial_Ward_'+str(k), None
            inverse, distances = invert_edges(q, sep, ordered=True)
            # B*B <=16(I-11*/8): every row norm square is <=14.
            diagonal = [ci(0) for _ in range(8)]
            row_lower = [0]*8
            for (i, j), zinv, d in zip(PAIRS, inverse, distances):
                z2 = cm(zinv, zinv)
                diagonal[i] = ca(diagonal[i], z2)
                diagonal[j] = ca(diagonal[j], z2)
                lower = S**3//(d[1]*d[1])
                row_lower[i] += lower
                row_lower[j] += lower
            if any(row_lower[i]+cnorm(diagonal[i])[0] > 14*S for i in range(8)):
                return None, 'Hessian_row_budget', None
            if direction is not None and rayleigh(q, inverse, direction)[1] < 0:
                return None, 'full_Hessian_Rayleigh', None
            tag = 'stationarity'
            gradient, force = force_gradient(q, inverse)
            if any(v[0] > 0 or v[1] < 0 for g in gradient for v in g):
                return None, 'force_sign', None
            for i in range(8):
                q[i] = meet(q[i][0], divide(force[i][0], 4)), meet(q[i][1], neg(divide(force[i][1], 4)))
            close_center(q)
            box = [meet(a, b) for a, b in zip(box, variables(q))]
            if box == old:
                break
        return box, None, (q, inverse)
    except (ValueError, ZeroDivisionError):
        return None, tag, None
