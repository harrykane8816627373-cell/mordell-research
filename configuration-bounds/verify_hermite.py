"""Exact supporting check for the six-point Hermite obstruction.

Only the Python standard library is required. This checks the explicitly
stated interval comparisons, not unrestricted six-point global optimality.
"""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
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



def certify_six_hermite_gap():
    def log_i(v):
        v=I(v)
        require(v.lo>=1,'log interval domain')
        return I(log_interval(v.lo,terms=24).lo,
                 log_interval(v.hi,terms=24).hi)
    rt=I(5).sqrt()
    nodes=[I(Q(6,5)),3*(5-rt)/5,3*(5+rt)/5]
    z=[a for a in nodes for _ in range(2)]
    table=[[log_i(a)] for a in z]
    for row in range(1,6):
        table[row].append(1/z[row] if row%2 else
            (table[row][0]-table[row-1][0])/(z[row]-z[row-1]))
        for col in range(2,row+1):
            table[row].append((table[row][col-1]-table[row-1][col-1])/
                              (z[row]-z[row-col]))
    coeff=[table[j][j] for j in range(6)]
    def H(x):
        x=I(x)
        y=coeff[5]
        for j in range(4,-1,-1):
            y=coeff[j]+(x-z[j])*y
        return y
    def P(x):
        x=Q(x)
        return (x-Q(6,5))**2*(x*x-6*x+Q(36,5))**2
    target=15*log_i(6)-10*log_i(5)
    gap=6*H(1)+6*H(3)+3*H(4)-target
    require(gap.lo>Q(2,25),'uncorrected Hermite obstruction')
    require(gap.hi<Q(83,1000),'upper obstruction enclosure')
    psum=6*P(1)+6*P(3)+3*P(4)
    require(psum==Q(396,5),'exact hexagon correction sum')
    radius=Q(15,2)
    kr=(H(radius)-log_i(radius))/P(radius)
    require(kr.lo>0,'positive finite-radius correction')
    corrected=gap-kr*psum
    require(corrected.lo>Q(3,100),'corrected Hermite obstruction')
    return {
        'status':'CERTIFIED',
        'scope':'Hermite obstruction checks only; no global-optimality certificate',
        'Hermite_gap_enclosure':[str(gap.lo),str(gap.hi)],
        'Hermite_gap_lower':'2/25','Hermite_gap_upper':'83/1000',
        'R':'15/2','corrected_gap_lower':'3/100',
        'corrected_gap_enclosure':[str(corrected.lo),str(corrected.hi)]}

if __name__=='__main__':
    result=certify_six_hermite_gap()
    Path(__file__).with_name('hermite_result.json').write_text(
        json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('CERTIFIED: Hermite gap > 2/25; corrected gap > 3/100.')
