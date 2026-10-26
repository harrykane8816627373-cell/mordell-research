"""Exact rational certificate for the five-point Mordell inequality.\nAuthor: Fengcan Jiang.

This certifies only the explicitly encoded scalar inequalities.  The mathematical
reduction from the original problem must be audited separately.
All mathematical bounds use rational arithmetic and integer square roots.
No third-party packages or external solver is used. Floating-point timestamps
are used only for progress reporting, never for a mathematical decision.
"""
from fractions import Fraction as Q
from math import isqrt
from time import monotonic
import argparse
import json


def require(condition, message):
    """Checks remain active under python -O."""
    if not condition:
        raise ValueError(message)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if isinstance(lo, I):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo = Q(lo)
            self.hi = Q(lo if hi is None else hi)
        require(self.lo <= self.hi, "ordered interval endpoints")

    def __add__(self, other):
        other = I(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-I(other))

    def __rsub__(self, other):
        return I(other) + (-self)

    def __mul__(self, other):
        other = I(other)
        products = (self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi)
        return I(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I(other)
        require(other.lo > 0, "positive interval divisor")
        return self * I(1 / other.hi, 1 / other.lo)

    def __rtruediv__(self, other):
        return I(other) / self

    def square(self):
        if self.lo >= 0:
            return I(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return I(self.hi * self.hi, self.lo * self.lo)
        return I(0, max(self.lo * self.lo, self.hi * self.hi))

    def positive(self):
        return I(max(0, self.lo), max(0, self.hi))

    def sqrt(self, bits=40):
        require(self.lo >= 0, "nonnegative square-root interval")
        scale = 1 << bits
        def floor_root(v):
            return isqrt((v.numerator * scale * scale) // v.denominator)
        low = floor_root(self.lo)
        high = floor_root(self.hi)
        if Q(high * high, scale * scale) < self.hi:
            high += 1
        return I(Q(low, scale), Q(high, scale))


C = Q(278365, 10000000)  # proved above the exact Hermite coefficient separately
K = Q(9, 2500)
AQ = Q(1, 10) - 3 * C + Q(102, 7) * K
BQ = Q(36, 1225) * K
DQ = Q(3, 7) * K


def log_interval(x, terms=120):
    """Enclose log(x), x>=1, using the positive atanh series."""
    x = Q(x)
    require(x >= 1, "log series domain")
    y = (x - 1) / (x + 1)
    y2 = y * y
    power = y
    total = Q(0)
    for j in range(terms):
        total += power / (2*j + 1)
        power *= y2
    tail = power / ((2*terms + 1) * (1-y2))
    return I(2*total, 2*(total+tail))


def certify_constants():
    n = 24
    series = sum((Q(1, (2*j+1)*5**j) for j in range(n)), Q(0))
    tail = Q(5, 4*(2*n+1)*5**n)
    c = Q(1, 5) - Q(4, 25)*I(series, series+tail)
    require(c.lo > Q(22, 3)*K, "p3 elimination")
    require(c.hi < C, "rational upper bound for c")
    require(c.hi < Q(17, 2)*K, "monotonicity in U")
    require(Q(800, 129) < Q(25, 4), "q cutoff")
    radius = Q(41, 6)
    require(((radius-4)/2)**2 > 2, "distance upper bound")
    h_at_radius = (log_interval(5)/2 - Q(7, 4) + radius-radius**2/10
                   + c*(-Q(25, 4)+15*radius-Q(15, 2)*radius**2+radius**3))
    margin = h_at_radius-log_interval(radius)-K*(radius**2-5*radius+5)**2
    require(margin.lo > 0, "quartic logarithm majorant")

    # Improved local lemma: 0<s<=1/20, pair radius products >7/8,
    # radius variance <1/70. Two angular bootstrap steps suffice.
    require(Q(1, 8)**2 > Q(3, 200), "pair radius product bounds")
    sqrt2 = I(2).sqrt()
    sqrt5 = I(5).sqrt()
    gamma = 2*I(35).sqrt()/5
    # Exact algebraic values of sin/cos(pi/5) and sin/cos(2*pi/5).
    # No pi approximation or trigonometric series is needed.
    sin1 = ((5-sqrt5)/8).sqrt()
    cos1 = (1+sqrt5)/4
    sin2 = ((5+sqrt5)/8).sqrt()
    cos2 = (sqrt5-1)/4
    def sine_lower(sine, cosine, u):
        return sine*(1-u*u/4)-cosine*sqrt2*u/2
    def cosine_upper(u):
        return cos1+sin1*sqrt2*u/2
    require((3*gamma/35).hi < Q(21, 100), "initial angular bound")
    steps = [Q(21,100), Q(9,100), Q(3,50)]
    local_bounds=[]
    for before, after in zip(steps, steps[1:]):
        sine = sine_lower(sin1,cos1,before)
        require(sine.lo>0,"positive algebraic sine lower bound")
        improved = gamma*cosine_upper(before)/(245*sine*sine*sine)
        require(improved.hi < after, "algebraic angular bootstrap")
        local_bounds.append({"input":str(before),"upper_below":str(after)})
    edge_weight = 1/(Q(7,2)*sine_lower(sin1,cos1,Q(3,50)).square())
    diagonal_weight = 1/(Q(7,2)*sine_lower(sin2,cos2,Q(3,50)).square())
    require(edge_weight.hi < Q(19,20), "edge weight")
    require(diagonal_weight.hi < Q(1,3), "diagonal weight")
    spectral_bound = Q(19,20)*(5+sqrt5)/2 + (5-sqrt5)/6
    require(spectral_bound.hi < 4, "local spectral contradiction")
    require(20*Q(1,20)**2 == Q(1,20), "exact local/global endpoint")
    return {"status": "CERTIFIED",
            "checks": "Hermite coefficient; quartic log majorant; two algebraic local angular steps; spectral gap; exact domain endpoint",
            "local_s_max": "1/20", "local_angular_steps": 2, "local_bounds":local_bounds}


def upper_expression(branch, x, q):
    """Interval enclosure of F after eliminating U, V, |p3|, |p4|."""
    if branch == "A":
        # s=20*x^2, 0<s<=5/4.
        x2 = x.square()
        s = 20 * x2
        t = 5 + s
        f = x2 * ((-10 + 450 * C - 500 * K)
                  + x * ((-300 * C + 3000 * K) - 8500 * K * x))
        d0 = 15 + x2 * (-180 + x * (5760 - 15600 * x))
    else:
        # x=s, 5/4<=s<=5/3.
        s = x
        t = 5 + s
        u = t.square() / 5
        radial = t * t.square() / 5 + 24 * t.square() - 275 * t + 750
        f = (-Q(1, 2) + Q(15, 2) * C) * s + 5 * C * (u - 5) - K * radial
        d0 = t * (-12 * t.square() + 144 * t - 405) / 5
    d22 = 20 - 12 * s - Q(4, 5) * q
    d44 = d0 + 5 * q
    if d22.hi < 0 or d44.hi < 0:
        return None
    # Restriction to the necessary domain d22,d44>=0 is valid here.
    radical = (q * d22.positive() * d44.positive()).sqrt()
    cutoff = (t - (Q(4, 7) + Q(9, 8) / s) * q).positive()
    return f - AQ * q + BQ * q.square() + DQ * radical - K * cutoff.square()


def certify_branch_b():
    """Uniform analytic bound on 5/4<=s<=5/3; no subdivision."""
    # On this branch f(s)=a*s+b*s^2-K*s^3/5 and
    # D_*(s)=15+27*s-36*s^2/5-12*s^3/5.
    a = -Q(1, 2) + Q(35, 2)*C + 20*K
    b = C - 27*K
    require(b < 0, "negative quadratic coefficient in branch B")
    require(a + Q(5, 2)*b < 0, "f decreasing for s>=5/4")
    require(27-Q(72,5)*Q(5,4)-Q(36,5)*Q(25,16) == -Q(9,4),
            "D_* decreasing for s>=5/4")
    # Hence f<=f(5/4), D_*<=525/16, q<=25/4, D22<=5,
    # and D_*+5*q<=1025/16. Drop the nonpositive cutoff term.
    f_left = a*Q(5,4) + b*Q(25,16) - K*Q(125,64)/5
    require(f_left == -Q(9181,256000), "branch B radial endpoint")
    alpha = AQ-Q(25,4)*BQ
    require(alpha > 0, "positive penalty after bounding q squared")
    # Completing the square in sqrt(q) bounds
    # -alpha*q+DQ*sqrt((5125/16)*q) by this rational number.
    correction = DQ*DQ*Q(5125,16)/(4*alpha)
    upper = f_left+correction
    require(upper == -Q(56657645489,1713169664000),
            "branch B exact uniform bound")
    require(upper < -Q(3,100), "strict branch B certificate")
    return {"branch":"B", "status":"CERTIFIED", "visited":1,
            "certified_boxes":1, "infeasible_boxes":0, "max_depth":0,
            "method":"uniform analytic bound by completing the square",
            "largest_leaf_upper":str(upper), "seconds":0.0}


def certify(branch, max_boxes=2000000):
    if branch == "B":
        return certify_branch_b()
    if branch == "A":
        left, right = Q(1, 20), Q(1, 4)
        xscale = Q(1, 5)
    else:
        left, right = Q(5, 4), Q(5, 3)
        xscale = Q(5, 12)
    todo = [(left, right, Q(0), Q(25), 0)]
    visited = certified = pruned = 0
    deepest = 0
    largest_leaf_upper = None
    start = monotonic()
    next_report = start + 20
    while todo:
        xl, xh, ql, qh, depth = todo.pop()
        visited += 1
        deepest = max(deepest, depth)
        smin = 20 * xl * xl if branch == "A" else xl
        smax = 20 * xh * xh if branch == "A" else xh
        # Projection Cauchy--Schwarz gives q<=80*s*(5+s)/(45+16*s).
        # This upper bound increases for s>=0, so smax is safe for a box.
        projection_cap = 80*smax*(5+smax)/(45+16*smax)
        qh = min(qh, 25 - 15 * smin, projection_cap)
        if ql > qh:
            pruned += 1
            continue
        bound = upper_expression(branch, I(xl, xh), I(ql, qh))
        if bound is None:
            pruned += 1
            continue
        if bound.hi < 0:
            certified += 1
            if largest_leaf_upper is None or bound.hi > largest_leaf_upper:
                largest_leaf_upper = bound.hi
            continue
        if visited >= max_boxes or depth >= 70:
            return {"branch": branch, "status": "UNRESOLVED", "visited": visited,
                    "remaining": len(todo), "box": list(map(str, (xl,xh,ql,qh))),
                    "upper": str(bound.hi)}
        if (xh - xl) / xscale >= (qh - ql):
            xm = (xl + xh) / 2
            todo.append((xm, xh, ql, qh, depth + 1))
            todo.append((xl, xm, ql, qh, depth + 1))
        else:
            qm = (ql + qh) / 2
            todo.append((xl, xh, qm, qh, depth + 1))
            todo.append((xl, xh, ql, qm, depth + 1))
        now = monotonic()
        if now >= next_report:
            print(json.dumps({"branch": branch, "visited": visited,
                              "pending": len(todo), "certified": certified,
                              "seconds": round(now-start,1)}), flush=True)
            next_report = now + 20
    require(largest_leaf_upper < -Q(1,200000),
            "uniform branch A certificate below -1/200000")
    return {"branch": branch, "status": "CERTIFIED", "visited": visited,
            "certified_boxes": certified, "infeasible_boxes": pruned,
            "max_depth": deepest, "largest_leaf_upper": str(largest_leaf_upper),
            "seconds": round(monotonic()-start, 3)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", choices=["A", "B", "both"], default="both")
    parser.add_argument("--max-boxes", type=int, default=2000000)
    parser.add_argument("--output")
    args = parser.parse_args()
    branches = ["A", "B"] if args.branch == "both" else [args.branch]
    result = [certify_constants()] + [certify(b, args.max_boxes) for b in branches]
    print(json.dumps(result, indent=2), flush=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
    if any(r["status"] != "CERTIFIED" for r in result):
        raise SystemExit(1)
