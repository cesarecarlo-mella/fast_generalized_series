#!/usr/bin/env python3
"""
build_letters.py -- dlog series of the active H letters around (x, y) = (1/2, 0)

    x = 1/2 + t ,  y = y          NO blow-up (see check_point.py)

    dlog_dt_c0.m : { l[i] -> d/dt Log[l_i] }      series in t, y
    dlog_dy_c0.m : { l[i] -> d/dy Log[l_i] }
    info_c0.m    : metadata

Truncation: t^i with i <= ORDT and y^e with e <= ORDY (e is a half-integer for
l7, l8), i.e. the same rectangular truncation the recursion uses.

  plain letter  P = y^m Q, Q(1/2, 0) != 0:
        d/dt log P = Q_x / Q ,   d/dy log P = m / y + Q_y / Q
  l7, l8        f = log((a - I r)/(a + I r)),  r = Sqrt[b],  b = y * x (1-x-y):
        d_v f = I (2 a_v b - a b_v) / (r (a^2 + b))
              = I y^(n - d - 1/2) N'/(D' Sqrt[x (1-x-y)])
        (N = y^n N', a^2 + b = y^d D', N'(1/2,0), D'(1/2,0) != 0; Sqrt[x(1-x-y)]
         -> 1/2 at the point, so everything is a plain power series times
         y^(half-integer)).  Branch: principal Sqrt, as in letters.dict.m.

Exact rationals throughout (gmpy2), then a numeric self-check against the
letters themselves (mpmath, 40 digits).

    python3 build_letters.py [ORDT] [ORDY] [--out ../pipeline/inputs]
"""
import argparse
import os
import sys
import time

import mpmath as mpm
from gmpy2 import mpq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))
import alphabet as al
from series2 import Ser
import mparse as mp

ap = argparse.ArgumentParser()
ap.add_argument('ordt', type=int, nargs='?', default=20)
ap.add_argument('ordy', type=int, nargs='?', default=20)
ap.add_argument('--out', default=os.path.join(HERE, '..', 'pipeline', 'inputs'))
a = ap.parse_args()
OT, OY = a.ordt, a.ordy
NYA = OY + 2                          # analytic part: 2 extra y orders, cut at the end
X0 = al.X0


def yval(p):
    return min(j for (_, j) in p)


def ydiv(p, m):
    return {(i, j - m): c for (i, j), c in p.items()}


def ser(p):
    return Ser.from_poly(p, X0, OT, NYA)


def to_poly(s, e2shift, imag, coef=1):
    """Ser * y^(e2shift/2) [* I] -> mparse Poly, truncated to t^OT, y^OY."""
    out = {}
    for i in range(OT + 1):
        for j in range(NYA + 1):
            c = s.c[i][j]
            if not c:
                continue
            e2y = 2 * j + e2shift
            if e2y > 2 * OY:
                continue
            assert e2y >= -2, "pole below y^-1"
            k = []
            if imag:
                k.append(('I', 2))
            if i:
                k.append(('t', 2 * i))
            if e2y:
                k.append(('y', e2y))
            out[tuple(sorted(k))] = c * coef
    return out


t0 = time.time()
DT, DY = {}, {}
for li in al.ACTIVE:
    if li in al.PLAIN:
        Pp = al.PLAIN[li]
        m = yval(Pp)
        Q = ydiv(Pp, m)
        q0 = al.peval(Q, X0, 0)
        assert q0 != 0, "l%d: Q(1/2,0) = 0 -> singular at the point" % li
        iQ = ser(Q).inv()
        DT[li] = to_poly(ser(al.pd(Q, 0)) * iQ, 0, False)
        dy = to_poly(ser(al.pd(Q, 1)) * iQ, 0, False)
        if m:
            dy[(('y', -2),)] = dy.get((('y', -2),), 0) + m
        DY[li] = {k: v for k, v in dy.items() if v}
    else:
        A, B = al.ODD[li]
        mb = yval(B)
        Bq = ydiv(B, mb)
        assert mb == 1, "l%d: b ~ y^%d" % (li, mb)
        D = al.padd(al.pmul(A, A), B)
        d = yval(D)
        Dq = ydiv(D, d)
        assert al.peval(Dq, X0, 0) != 0, "l%d: a^2 + b vanishes beyond y^%d at the point" % (li, d)
        iden = (ser(Dq) * ser(Bq).sqrt()).inv()
        for var, store in ((0, DT), (1, DY)):
            N = al.padd(al.pmul(al.pmul({(0, 0): mpq(2)}, al.pd(A, var)), B), al.pmul(A, al.pd(B, var)), -1)
            n = yval(N) if N else 0
            Nq = ydiv(N, n) if N else {}
            e2 = 2 * (n - d) - mb                    # y^(n - d - 1/2)
            store[li] = to_poly(ser(Nq) * iden, e2, True) if N else {}
    print("  l%-2d  dt: %4d terms   dy: %4d terms" % (li, len(DT[li]), len(DY[li])), flush=True)
print("built in %.1fs" % (time.time() - t0))

# ---------------------------------------------------------------- self-check
mpm.mp.dps = 40
pt = (mpq(1, 50), mpq(3, 100))       # (t, y): |y| well inside the radius 1/4


def evalpoly(p, t, y):
    s = mpm.mpc(0)
    for k, c in p.items():
        term = mpm.mpf(c.numerator) / c.denominator
        for at, e in k:
            if at == 'I':
                term *= 1j
            elif at == 't':
                term *= mpm.mpf(t) ** (mpm.mpf(e) / 2)
            elif at == 'y':
                term *= mpm.mpf(y) ** (mpm.mpf(e) / 2)
        s += term
    return s


worst = 0
tt, yy = mpm.mpf(pt[0].numerator) / pt[0].denominator, mpm.mpf(pt[1].numerator) / pt[1].denominator
for li in al.ACTIVE:
    f = lambda T, Yv: al.letter_value(li, mpm.mpf(1) / 2 + T, Yv)
    rt = mpm.diff(lambda T: f(T, yy), tt)
    ry = mpm.diff(lambda Yv: f(tt, Yv), yy)
    for name, ser_, ref in (('dt', DT[li], rt), ('dy', DY[li], ry)):
        v = evalpoly(ser_, tt, yy)
        rel = abs(v - ref) / max(abs(ref), mpm.mpf(10) ** -30)
        worst = max(worst, rel)
        if rel > 1e-12:
            print("  !! l%d %s: series %s  vs  letter %s  (rel %.1e)" % (li, name, mpm.nstr(v, 15), mpm.nstr(ref, 15), rel))
print("numeric check at (t, y) = (1/50, 3/100): max rel. difference %.1e" % worst)
if worst > 1e-12:
    sys.exit("numeric check FAILED")

# ---------------------------------------------------------------- write
os.makedirs(a.out, exist_ok=True)
for name, S in (('dlog_dt_c0.m', DT), ('dlog_dy_c0.m', DY)):
    with open(os.path.join(a.out, name), 'w') as fh:
        fh.write("{" + ",\n ".join("l[%d] -> %s" % (li, mp.fmt_poly(S[li]) if S[li] else "0")
                                   for li in al.ACTIVE) + "}\n")
with open(os.path.join(a.out, 'info_c0.m'), 'w') as fh:
    fh.write('<|"point" -> {x -> 1/2, y -> 0}, "expansion" -> "x = 1/2 + t, y = y (no blow-up)", '
             '"corner" -> 0, "ORDt" -> %d, "ORDy" -> %d, "active" -> {%s}, "maxReldiff" -> %s|>\n'
             % (OT, OY, ", ".join(map(str, al.ACTIVE)), mpm.nstr(worst, 3)))
print("wrote dlog_dt_c0.m, dlog_dy_c0.m, info_c0.m -> %s" % os.path.abspath(a.out))
