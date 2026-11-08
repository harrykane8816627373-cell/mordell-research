"""Integer-only outward intervals on the fixed 2**384 lattice."""
PRECISION = 384
SCALE = 1 << PRECISION


class Ball:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi
        if type(self.lo) is not int or type(self.hi) is not int or self.lo > self.hi:
            raise ValueError('ordered integer endpoints required')

    @classmethod
    def dyadic(cls, numerator, exponent):
        if type(numerator) is not int or type(exponent) is not int or exponent < 0:
            raise ValueError('integer numerator and nonnegative exponent required')
        if exponent <= PRECISION:
            return cls(numerator << (PRECISION-exponent))
        divisor = 1 << (exponent-PRECISION)
        return cls(numerator//divisor, -((-numerator)//divisor))

    @classmethod
    def integer(cls, n):
        return cls(n*SCALE)

    @staticmethod
    def convert(value):
        if isinstance(value, Ball):
            return value
        if type(value) is int:
            return Ball.integer(value)
        raise TypeError('only integers and dyadic balls are accepted')

    def __add__(self, other):
        b = self.convert(other)
        return Ball(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Ball(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.convert(other))

    def __rsub__(self, other):
        return self.convert(other) + (-self)

    def __mul__(self, other):
        b = self.convert(other)
        products = [x*y for x in (self.lo, self.hi) for y in (b.lo, b.hi)]
        return Ball(min(products)//SCALE, -((-max(products))//SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.convert(other)
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError('denominator interval contains zero')
        fractions = [(a*SCALE, d) for a in (self.lo, self.hi) for d in (b.lo, b.hi)]
        return Ball(min(a//d for a, d in fractions), max(-((-a)//d) for a, d in fractions))

    def __rtruediv__(self, other):
        return self.convert(other)/self

    def square(self):
        small = 0 if self.lo <= 0 <= self.hi else min(self.lo*self.lo, self.hi*self.hi)
        large = max(self.lo*self.lo, self.hi*self.hi)
        return Ball(small//SCALE, -((-large)//SCALE))

    def magnitude(self):
        return max(abs(self.lo), abs(self.hi))


def matmul(a, b):
    rows, shared, columns = len(a), len(b), len(b[0])
    if any(len(row) != shared for row in a) or any(len(row) != columns for row in b):
        raise ValueError('rectangular conforming matrices required')
    out = [[Ball.integer(0) for _ in range(columns)] for _ in range(rows)]
    for i in range(rows):
        for k in range(shared):
            v = a[i][k]
            if v.lo == v.hi == 0:
                continue
            for j in range(columns):
                out[i][j] = out[i][j]+v*b[k][j]
    return out


def transpose(a):
    return [list(row) for row in zip(*a)]
