#!/usr/bin/env python3
"""
build_cc_letters.py -- dlog series of the H letters at the corners CC2, CC3
of the continued region, in blow-up variables (exact rationals).

Chart (analytic continuation, x,y < 0, z > 1):
    x = -w/v ,  y = -u/v ,  z = 1/v ,     u + w + v = 1
In (w,u,v) the alphabet closes on itself (see ../../alphabet_map/), with
(x,y,z) -> (w,u,v) relabelled, so the corners are blown up exactly as the
old corners 2 and 3 (00_build_letters.wl, W = 2):

    CC2 : u, v -> 0 (w -> 1)    u = t*v^2 ,  BIG = v ,  w = 1 - v - t v^2
          (u << v, mirror of CC3; NOT the old-c2 shape v = t u^2, see SPEC)
    CC3 : w, v -> 0 (u -> 1)    w = t*v^2 ,  BIG = v ,  u = 1 - v - t v^2
          (old c3: x = t*z^2, BIG = z, y = 1 - x - z)

The letters are the ORIGINAL letters l_k(x, y) of letters.dict.m, composed
with the chart (so the ORIGINAL atilde is used unchanged):

    dlog_dt_<CC>.m : { l[k] -> d/dt   Log[l_k(x(t,B), y(t,B))] }
    dlog_dv_<CC>.m : { l[k] -> d/dv   Log[l_k(x(t,v), y(t,v))] }   (BIG = v for both)

  plain letter : arg = N/D, N = t^a B^b Q, Q(0,0) != 0 (asserted -- this is
                 the statement that the blow-up resolves the letter)
                 dlog = monomial poles + Q'/Q  (same for D)
  l7, l8       : f = Log[(a - I r)/(a + I r)], r = Sqrt[b], b = x y z,
                 a = x z (l7), a = x y (l8)
                 d f = I (2 a' b - a b') / (r (a^2 + b))
                 in the continued region x y z > 0; r = -Sqrt|xyz| by default
                 (see BRANCH below).

Truncation: t^e with e <= ORDT, B^e with e <= ORDB (e half-integer for
l7, l8); the same rectangular truncation the recursion uses.

A numeric self-check (mpmath, 40 digits) compares every series with the
derivative of the true composed letter at one point (t, B) = (3/100, 7/100).

    python3 build_cc_letters.py [ORDT] [ORDB] [--out ../inputs]

BRANCH (l7, l8): the letters are odd under r -> -r.  The only H masters
that are odd (prefactor Sqrt) are 156, 157; they have zero boundary
constant everywhere, and even masters only see products (odd letter) x
(odd master), so the choice r = +Sqrt|xyz| vs -Sqrt|xyz| only flips the
overall sign of masters 156, 157.  DEFAULT: --odd-sign -1, r = -Sqrt|xyz|,
the choice consistent with the continuation Log x -> log|x| + I Pi,
Log y -> log|y| + I Pi used for the CC boundaries:
    Sqrt[x] Sqrt[y] Sqrt[z] = (I Sqrt|x|)(I Sqrt|y|) Sqrt[z] = -Sqrt|xyz|.
--odd-sign +1 gives the principal root.  (Note: ../../notes/galois_parity.tex)
"""
import argparse
import os
import sys
import time

import mpmath as mpm
import sympy as sp
from gmpy2 import mpq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from series2 import Ser
import mparse as mp

ap = argparse.ArgumentParser()
ap.add_argument('ordt', type=int, nargs='?', default=20)
ap.add_argument('ordb', type=int, nargs='?', default=60)
ap.add_argument('--out', default=os.path.join(HERE, '..', 'inputs'))
ap.add_argument('--corners', default='CC2,CC3')
ap.add_argument('--odd-sign', type=int, default=-1, choices=[1, -1],
                help='sign of r = Sqrt[xyz] in the continued region: -1 (default) r = -Sqrt|xyz|, '
                     'consistent with Log x, Log y -> log|.| + I Pi; +1 principal root')
a = ap.parse_args()
OT, OB = a.ordt, a.ordb
XT, XB = 4, 8                      # extra orders kept internally, cut at the end

x, y, t, B = sp.symbols('x y t B')
z = 1 - x - y
I = sp.I

# --------------------------------------------------------------- alphabet
PLAIN = {1: x, 2: y, 3: 1 - x - y, 4: 1 - x, 5: 1 - y, 6: x + y,
         9: 1 - (1 + x) * y, 10: 1 - x * (1 + y), 11: x - y + y**2, 12: -x + x**2 + y,
         13: (x - 1)**2 - y, 14: -y + (x + y)**2, 15: -x + (x + y)**2, 16: -x + (1 - y)**2,
         17: x + x * y + y**2, 18: x**2 + y + x * y, 19: (y - 1)**2 + x * y, 20: (x - 1)**2 + x * y}
ODD = {7: x * z, 8: x * y}          # a; b = x y z
ACTIVE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16]   # letters in H_atilde


def letter_value(k, X, Y):
    """true letter l_k(X, Y) (mpmath), principal branches as in letters.dict.m."""
    Z = 1 - X - Y
    if k in PLAIN:
        return mpm.log(sp.lambdify((x, y), PLAIN[k], 'mpmath')(X, Y))
    r = mpm.sqrt(X * Z * Y)
    A = X * Z if k == 7 else X * Y
    return mpm.log((A - 1j * r) / (A + 1j * r))


def _from_xy(xe, ye):
    """chart given by x(t,B), y(t,B) -> u, w, v (v = 1/z, w = -x v, u = -y v)."""
    ve = 1 / (1 - xe - ye)
    return dict(w=sp.cancel(-xe * ve), u=sp.cancel(-ye * ve), v=sp.cancel(ve))


SPEC = {
    # CC2: u << v (same ordering as the transport of bdy_CC2: u -> 0 first).
    # The old-c2-shaped chart v = t u^2 (v << u) is WRONG here: from weight 2 the
    # masters depend on the ratio u/v (y = -u/v) and the constant term depends on
    # the ordering.  See claude/analytic_continuation.md.
    'CC2': dict(BIG='v', SMALL='u', blowup='u = t*v^2',
                w=1 - B - t * B**2, u=t * B**2, v=B, old=0),
    # edge-side charts (v << u at CC2, v << w at CC3): they reach the v = 0 edge,
    # where |x|, |y| > 1.  Same shape as the old c2 chart (CC2e) and its mirror (CC3e).
    'CC2e': dict(BIG='u', SMALL='v', blowup='v = t*u^2',
                 w=1 - B - t * B**2, u=B, v=t * B**2, old=2),
    'CC3e': dict(BIG='w', SMALL='v', blowup='v = t*w^2',
                 w=B, u=1 - B - t * B**2, v=t * B**2, old=0),
    # ratio-line charts: the bands |y| = 1 (near CC2) and |x| = 1 (near CC3) are where the
    # ratio towers of CC2/CC2e (CC3/CC3e) are at their radius (singularity at ratio -1).
    # Expand around the point x = -inf, y = -1 (resp. y = -inf, x = -1) in coordinates
    # that make the letter 1 - x(1+y) (resp. 1 - y(1+x)) a coordinate line:
    #   CC2r:  x = -1/p ,  y = -1 - p + t     (1 - x(1+y) = t/p)
    #   CC3r:  y = -1/p ,  x = -1 - p + t     (1 - y(1+x) = t/p)
    # t = 0 is a spurious singularity of the masters: no Log[t] may survive.
    'CC2r': dict(BIG='p', SMALL='t', blowup='x = -1/p, y = -1 - p + t (no blow-up)', old=0,
                 **_from_xy(-1 / B, -1 - B + t)),
    'CC3r': dict(BIG='p', SMALL='t', blowup='y = -1/p, x = -1 - p + t (no blow-up)', old=0,
                 **_from_xy(-1 - B + t, -1 / B)),
    'CC3': dict(BIG='v', SMALL='w', blowup='w = t*v^2',
                w=t * B**2, u=1 - B - t * B**2, v=B, old=3),
}


# --------------------------------------------------------------- helpers
def monomial_split(p):
    """polynomial p(t,B) -> (a, b, Q) with p = t^a B^b Q, Q(0,0) = const term."""
    P = sp.Poly(sp.expand(p), t, B)
    a_ = min(m[0] for m in P.monoms())
    b_ = min(m[1] for m in P.monoms())
    Q = sp.expand(P.as_expr() / (t**a_ * B**b_))
    return a_, b_, Q


def to_ser(p, NT, NB):
    P = sp.Poly(sp.expand(p), t, B)
    s = Ser(NT, NB)
    for (i, j), c in zip(P.monoms(), P.coeffs()):
        if i <= NT and j <= NB:
            c = sp.Rational(c)
            s.c[i][j] += mpq(int(c.p), int(c.q))
    return s


def to_poly(s, e2t, e2b, imag, BIGNAME, coef=1):
    """Ser * t^(e2t/2) B^(e2b/2) [* I] -> mparse Poly, truncated to t^OT, B^OB."""
    out = {}
    for i in range(s.NT + 1):
        for j in range(s.NY + 1):
            c = s.c[i][j]
            if not c:
                continue
            et, eb = 2 * i + e2t, 2 * j + e2b
            if et > 2 * OT or eb > 2 * OB:
                continue
            assert et >= -2 and eb >= -2, "pole beyond 1/t or 1/B"
            k = []
            if imag:
                k.append(('I', 2))
            if et:
                k.append(('t', et))
            if eb:
                k.append((BIGNAME, eb))
            key = tuple(sorted(k))
            out[key] = out.get(key, 0) + c * coef
    return {k: v for k, v in out.items() if v}


def padd(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}


def plain_dlog(arg, BIGNAME):
    """dlog_t, dlog_B of a rational arg(t, B)."""
    num, den = sp.fraction(sp.cancel(sp.together(arg)))
    DT, DB = {}, {}
    for part, sgn in ((num, 1), (den, -1)):
        if not part.free_symbols:
            continue
        at, bb, Q = monomial_split(part)
        q0 = Q.subs({t: 0, B: 0})
        assert q0 != 0, "blow-up does not resolve this letter: Q(0,0) = 0  (%s)" % Q
        NT, NB = OT + XT, OB + XB
        iQ = to_ser(Q, NT, NB).inv()
        dt = to_poly(to_ser(sp.diff(Q, t), NT, NB) * iQ, 0, 0, False, BIGNAME, sgn)
        db = to_poly(to_ser(sp.diff(Q, B), NT, NB) * iQ, 0, 0, False, BIGNAME, sgn)
        if at:
            dt = padd(dt, {(('t', -2),): mpq(sgn * at)})
        if bb:
            db = padd(db, {((BIGNAME, -2),): mpq(sgn * bb)})
        DT, DB = padd(DT, dt), padd(DB, db)
    return DT, DB


def odd_dlog(A, Bq, BIGNAME, sign):
    """d f, f = Log[(A - I r)/(A + I r)], r = sign*Sqrt[Bq] (principal Sqrt > 0)."""
    nb, db_ = sp.fraction(sp.cancel(sp.together(Bq)))
    pn, qn, hn = monomial_split(nb)
    pd_, qd, hd = monomial_split(db_)
    c_n, c_d = hn.subs({t: 0, B: 0}), hd.subs({t: 0, B: 0})
    if c_n < 0:                                 # b > 0: flip both signs
        hn, hd, c_n, c_d = -hn, -hd, -c_n, -c_d
    assert c_n > 0 and c_d > 0, "b does not have a positive limit shape"
    # r = t^((pn-pd)/2) B^((qn-qd)/2) Sqrt[hn]/Sqrt[hd]
    NT, NB = OT + XT, OB + XB
    sq_n, sq_d = to_ser(hn, NT, NB).sqrt(), to_ser(hd, NT, NB).sqrt()
    out = []
    for var in (t, B):
        g = sp.cancel(sp.together((2 * sp.diff(A, var) * Bq - A * sp.diff(Bq, var)) / (A**2 + Bq)))
        gn, gd = sp.fraction(g)
        if gn == 0:
            out.append({}); continue
        an, bn, Gn = monomial_split(gn)
        ad, bd, Gd = monomial_split(gd)
        assert Gd.subs({t: 0, B: 0}) != 0, "a^2 + b: non-monomial zero at the corner"
        # d f = I * g / r = I * t^e2t/2 B^e2b/2 * Gn * Sqrt[hd] / (Gd * Sqrt[hn])
        e2t = 2 * (an - ad) - (pn - pd_)
        e2b = 2 * (bn - bd) - (qn - qd)
        S = to_ser(Gn, NT, NB) * sq_d * (to_ser(Gd, NT, NB) * sq_n).inv()
        out.append(to_poly(S, e2t, e2b, True, BIGNAME, sign))
    return out


def evalpoly(p, tv, bv, BIGNAME):
    s = mpm.mpc(0)
    for k, c in p.items():
        term = mpm.mpf(c.numerator) / c.denominator
        for at, e in k:
            if at == 'I':
                term *= 1j
            elif at == 't':
                term *= tv ** (mpm.mpf(e) / 2)
            elif at == BIGNAME:
                term *= bv ** (mpm.mpf(e) / 2)
            else:
                raise ValueError(at)
        s += term
    return s


# --------------------------------------------------------------- build
os.makedirs(a.out, exist_ok=True)
mpm.mp.dps = 40
for CC in a.corners.split(','):
    S = SPEC[CC]
    BIGNAME = S['BIG']
    t0 = time.time()
    sub = {x: -S['w'] / S['v'], y: -S['u'] / S['v']}
    print("=== %s : %s, BIG = %s, x = -w/v, y = -u/v,  order t^%d %s^%d" % (CC, S['blowup'], BIGNAME, OT, BIGNAME, OB), flush=True)
    DT, DB = {}, {}
    for k in ACTIVE:
        if k in PLAIN:
            DT[k], DB[k] = plain_dlog(PLAIN[k].subs(sub, simultaneous=True), BIGNAME)
        else:
            A = ODD[k].subs(sub, simultaneous=True)
            Bq = (x * y * z).subs(sub, simultaneous=True)
            DT[k], DB[k] = odd_dlog(A, Bq, BIGNAME, a.odd_sign)
        print("  l%-2d  dt: %5d terms   d%s: %5d terms" % (k, len(DT[k]), BIGNAME, len(DB[k])), flush=True)
    print("  built in %.1fs" % (time.time() - t0), flush=True)

    # ---- numeric self-check vs the true composed letter
    tv, bv = mpm.mpf(3) / 100, mpm.mpf(7) / 100
    chart = sp.lambdify((t, B), (-S['w'] / S['v'], -S['u'] / S['v']), 'mpmath')

    def Ltrue(k, T, Bv):
        X, Y = chart(T, Bv)
        return letter_value(k, X, Y)

    worst = 0
    for k in ACTIVE:
        rt = mpm.diff(lambda T: Ltrue(k, T, bv), tv)
        rb = mpm.diff(lambda Bv: Ltrue(k, tv, Bv), bv)
        if k in ODD:
            rt, rb = rt * a.odd_sign, rb * a.odd_sign
        for name, ser_, ref in (('dt', DT[k], rt), ('d' + BIGNAME, DB[k], rb)):
            v = evalpoly(ser_, tv, bv, BIGNAME)
            rel = abs(v - ref) / max(abs(ref), mpm.mpf(10) ** -30)
            worst = max(worst, rel)
            if rel > 1e-12:
                print("  !! l%d %s: series %s  vs  letter %s  (rel %.1e)" % (k, name, mpm.nstr(v, 15), mpm.nstr(ref, 15), rel))
    print("  numeric check at (t, %s) = (3/100, 7/100): max rel. difference %s" % (BIGNAME, mpm.nstr(worst, 3)), flush=True)
    if worst > 1e-12:
        sys.exit("numeric check FAILED")

    for name, D in (('dlog_dt_%s.m' % CC, DT), ('dlog_d%s_%s.m' % (BIGNAME, CC), DB)):
        with open(os.path.join(a.out, name), 'w') as fh:
            fh.write("{" + ",\n ".join("l[%d] -> %s" % (k, mp.fmt_poly(D[k]) if D[k] else "0")
                                       for k in ACTIVE) + "}\n")
    with open(os.path.join(a.out, 'info_%s.m' % CC), 'w') as fh:
        fh.write('<|"region" -> "%s", "chart" -> "x = -w/v, y = -u/v, z = 1/v, u + w + v = 1", '
                 '"blowup" -> "%s", "BIG" -> %s, "SMALL" -> %s, "oldCornerShape" -> %d, '
                 '"ORDt" -> %d, "ORDy" -> %d, "active" -> {%s}, "oddSign" -> %d, "maxReldiff" -> %s|>\n'
                 % (CC, S['blowup'], BIGNAME, S['SMALL'], S['old'], OT, OB,
                    ", ".join(map(str, ACTIVE)), a.odd_sign, mpm.nstr(worst, 3)))
    print("  wrote dlog_dt_%s.m, dlog_d%s_%s.m, info_%s.m -> %s" % (CC, BIGNAME, CC, CC, os.path.abspath(a.out)), flush=True)
