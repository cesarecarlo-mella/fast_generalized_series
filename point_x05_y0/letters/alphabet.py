"""
alphabet.py -- the active H letters (letters.dict.m) as exact polynomials in
(x, y), plus the point we expand around.

    x = X0 + t ,  y = y          (X0 = 1/2, y0 = 0)

Plain letters are  l_i = Log[P_i(x, y)]  with P_i a polynomial.
l7, l8 are  Log[(a - I Sqrt[b]) / (a + I Sqrt[b])]  with a, b polynomials:
    l7 : a = x (1-x-y) ,  b = x (1-x-y) y
    l8 : a = x y       ,  b = x (1-x-y) y
Polynomials are dicts {(i, j): mpq}  meaning  sum c x^i y^j.
"""
from gmpy2 import mpq

X0 = mpq(1, 2)
ACTIVE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16]


def P(*terms):
    """P((c, i, j), ...) -> {(i, j): c}"""
    out = {}
    for c, i, j in terms:
        out[(i, j)] = out.get((i, j), 0) + mpq(c)
    return {k: v for k, v in out.items() if v}


def pmul(a, b):
    out = {}
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            k = (i1 + i2, j1 + j2)
            out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}


def padd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def pd(a, var):
    """d/dx (var=0) or d/dy (var=1)"""
    out = {}
    for (i, j), c in a.items():
        e = (i, j)[var]
        if e:
            k = (i - 1, j) if var == 0 else (i, j - 1)
            out[k] = out.get(k, 0) + c * e
    return out


def peval(a, x, y):
    return sum(c * x ** i * y ** j for (i, j), c in a.items())


X, Y, ONE = P((1, 1, 0)), P((1, 0, 1)), P((1, 0, 0))
Z = P((1, 0, 0), (-1, 1, 0), (-1, 0, 1))                  # 1 - x - y

PLAIN = {
    1: X,
    2: Y,
    3: Z,
    4: P((1, 0, 0), (-1, 1, 0)),                           # 1 - x
    5: P((1, 0, 0), (-1, 0, 1)),                           # 1 - y
    6: P((1, 1, 0), (1, 0, 1)),                            # x + y
    9: P((1, 0, 0), (-1, 0, 1), (-1, 1, 1)),               # 1 - (1+x) y
    10: P((1, 0, 0), (-1, 1, 0), (-1, 1, 1)),              # 1 - x (1+y)
    11: P((1, 1, 0), (-1, 0, 1), (1, 0, 2)),               # x - y + y^2
    12: P((-1, 1, 0), (1, 2, 0), (1, 0, 1)),               # -x + x^2 + y
    13: P((1, 0, 0), (-2, 1, 0), (1, 2, 0), (-1, 0, 1)),   # (x-1)^2 - y
    16: P((-1, 1, 0), (1, 0, 0), (-2, 0, 1), (1, 0, 2)),   # -x + (1-y)^2
}
B_SQ = pmul(pmul(X, Z), Y)                                 # b = x (1-x-y) y
ODD = {7: (pmul(X, Z), B_SQ),                              # (a, b)
       8: (pmul(X, Y), B_SQ)}


def letter_value(i, x, y):
    """numeric letter (mpmath), same branch as letters.dict.m"""
    import mpmath as m
    if i in PLAIN:
        return m.log(peval(PLAIN[i], x, y))
    a, b = ODD[i]
    av, r = peval(a, x, y), m.sqrt(peval(b, x, y))
    return m.log((av - 1j * r) / (av + 1j * r))
