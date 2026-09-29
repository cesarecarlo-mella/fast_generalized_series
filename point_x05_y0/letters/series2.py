"""
series2.py -- exact truncated double power series in (t, y) with mpq
coefficients:  S[i][j] = coefficient of t^i y^j,  0 <= i <= NT, 0 <= j <= NY.
"""
import math
from gmpy2 import mpq


class Ser:
    def __init__(self, NT, NY, c=None):
        self.NT, self.NY = NT, NY
        self.c = c if c is not None else [[mpq(0)] * (NY + 1) for _ in range(NT + 1)]

    @classmethod
    def from_poly(cls, p, x0, NT, NY):
        """polynomial {(i, j): c} in (x, y), with x = x0 + t."""
        s = cls(NT, NY)
        for (i, j), c in p.items():
            if j > NY:
                continue
            for k in range(min(i, NT) + 1):          # (x0 + t)^i = sum C(i,k) x0^(i-k) t^k
                s.c[k][j] += c * math.comb(i, k) * x0 ** (i - k)
        return s

    def __add__(self, o):
        return Ser(self.NT, self.NY, [[a + b for a, b in zip(r1, r2)] for r1, r2 in zip(self.c, o.c)])

    def __sub__(self, o):
        return Ser(self.NT, self.NY, [[a - b for a, b in zip(r1, r2)] for r1, r2 in zip(self.c, o.c)])

    def scale(self, q):
        return Ser(self.NT, self.NY, [[a * q for a in r] for r in self.c])

    def __mul__(self, o):
        NT, NY = self.NT, self.NY
        A, B = self.c, o.c
        nzA = [(i, j, A[i][j]) for i in range(NT + 1) for j in range(NY + 1) if A[i][j]]
        out = [[mpq(0)] * (NY + 1) for _ in range(NT + 1)]
        for i, j, a in nzA:
            for k in range(NT + 1 - i):
                Bk, ok = B[k], out[i + k]
                for l in range(NY + 1 - j):
                    b = Bk[l]
                    if b:
                        ok[j + l] += a * b
        return Ser(NT, NY, out)

    def inv(self):
        """1/S, needs S(0,0) != 0 (coefficient recursion, exact)."""
        NT, NY, A = self.NT, self.NY, self.c
        a0 = A[0][0]
        assert a0 != 0, "inverse of a series vanishing at the origin"
        ia0 = 1 / a0
        C = [[mpq(0)] * (NY + 1) for _ in range(NT + 1)]
        nzA = [(k, l, A[k][l]) for k in range(NT + 1) for l in range(NY + 1) if A[k][l] and (k, l) != (0, 0)]
        for i in range(NT + 1):
            for j in range(NY + 1):
                s = mpq(1) if (i, j) == (0, 0) else mpq(0)
                for k, l, a in nzA:
                    if k <= i and l <= j:
                        s -= a * C[i - k][j - l]
                C[i][j] = s * ia0
        return Ser(NT, NY, C)

    def sqrt(self):
        """Sqrt[S] for S(0,0) = a perfect rational square > 0 (principal root)."""
        NT, NY, A = self.NT, self.NY, self.c
        a0 = A[0][0]
        n, d = a0.numerator, a0.denominator
        import gmpy2
        rn, en = gmpy2.iroot(n, 2); rd, ed = gmpy2.iroot(d, 2)
        assert en and ed and a0 > 0, "constant term %s is not a rational square" % a0
        s0 = mpq(rn, rd)
        C = [[mpq(0)] * (NY + 1) for _ in range(NT + 1)]
        C[0][0] = s0
        i2s0 = 1 / (2 * s0)
        for i in range(NT + 1):
            for j in range(NY + 1):
                if (i, j) == (0, 0):
                    continue
                s = A[i][j]
                for k in range(i + 1):
                    for l in range(j + 1):
                        if (k, l) in ((0, 0), (i, j)):
                            continue
                        s -= C[k][l] * C[i - k][j - l]
                C[i][j] = s * i2s0
        return Ser(NT, NY, C)

    def value(self, t, y):
        return sum(self.c[i][j] * t ** i * y ** j
                   for i in range(self.NT + 1) for j in range(self.NY + 1) if self.c[i][j])
