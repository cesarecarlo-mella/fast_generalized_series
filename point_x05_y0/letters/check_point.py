#!/usr/bin/env python3
"""
check_point.py -- TASK 1: is a blow-up needed at (x, y) = (1/2, 0)?

A blow-up is needed when two or more singular curves of the alphabet meet at
the expansion point (or one of them is singular/tangent there): then no
double series in (x - 1/2, y) exists. Here we list, for every active letter,
the curves where its dlog is singular (zeros of P_i; for l7, l8 the zeros of
b and of a^2 + b), and check which of them pass through the point.

Also: (1) the alphabet in alphabet.py is checked numerically against
letters.dict.m, (2) the nearest singularities along the two axes are printed
(they bound where the series converge).

    python3 check_point.py [path/to/letters.dict.m]
"""
import os
import re
import sys

import numpy as np
import mpmath as mpm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import alphabet as al

X0, Y0 = al.X0, 0

# ---- (1) alphabet.py == letters.dict.m (numerically, mod 2 pi i) ------------
ld = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'letters.dict.m')
if os.path.isfile(ld):
    s = open(ld).read()
    s = re.sub(r'\(\*.*?\*\)', '', s, flags=re.S)
    rules = dict((int(m.group(1)), m.group(2)) for m in
                 re.finditer(r'l\[(\d+)\]\s*->\s*(Log\[.*?\])\s*(?=,\s*l\[|\}\s*$)', s.strip(), re.S))
    def pyexpr(e):
        e = e.replace('Log[', 'log(').replace('Sqrt[', 'sqrt(').replace(']', ')').replace('^', '**')
        return re.sub(r'\bI\b', '1j', e)
    worst = 0
    for (xv, yv) in ((0.31, 0.12), (0.45, 0.05), (0.2, 0.6)):
        for i in al.ACTIVE:
            ref = eval(pyexpr(rules[i]), {'log': mpm.log, 'sqrt': mpm.sqrt, 'x': mpm.mpf(xv), 'y': mpm.mpf(yv)})
            v = al.letter_value(i, mpm.mpf(xv), mpm.mpf(yv))
            d = abs(mpm.exp(v - ref) - 1)          # equal up to 2 pi i
            worst = max(worst, d)
    print("alphabet.py vs letters.dict.m: max |exp(l - l_ref) - 1| = %.1e  %s"
          % (worst, "OK" if worst < 1e-12 else "MISMATCH"))
else:
    print("(letters.dict.m not found next to this script -- skipping the alphabet check)")

# ---- (2) which singular curves pass through the point? --------------------
print("\nexpansion point (x, y) = (1/2, 0):")
curves = []                                    # (letter, name, poly)
for i in al.ACTIVE:
    if i in al.PLAIN:
        curves.append((i, "P", al.PLAIN[i]))
    else:
        a, b = al.ODD[i]
        curves.append((i, "b", b))
        curves.append((i, "a^2+b", al.padd(al.pmul(a, a), b)))
through = []
for i, name, p in curves:
    v = al.peval(p, X0, 0)
    # y-adic valuation: which power of y divides p, and is the rest regular?
    m = min(j for (_, j) in p)
    rest = {(k, j - m): c for (k, j), c in p.items()}
    r0 = al.peval(rest, X0, 0)
    flag = "" if v != 0 else "   <-- vanishes: = y^%d * (%s at the point)" % (m, r0)
    if v == 0:
        through.append((i, name, m, r0))
    print("  l%-2d %-6s at point = %-8s%s" % (i, name, v, flag))

print("\nsingular curves through the point:")
for i, name, m, r0 in through:
    print("  l%-2d %-6s : y^%d times a factor that is %s there" % (i, name, m, "NONZERO" if r0 else "ZERO"))
only_y = all(r0 != 0 for *_, r0 in through)
print("\n=> " + ("the ONLY singular curve through (1/2, 0) is y = 0 (log y, and Sqrt[y] in l7, l8):\n"
                 "   normal crossing trivially holds, NO blow-up is needed. The dlogs are\n"
                 "   y^(-1) or y^(-1/2) times double power series in (x - 1/2, Sqrt[y])."
                 if only_y else "another curve passes through the point: a blow-up IS needed."))

# ---- (3) nearest singularities along the axes ------------------------------
def roots_1d(p, var):
    """roots of p(1/2 + t, 0) in t (var=0) or of p(1/2, y) in y (var=1)."""
    deg = max(k[var] for k in p)
    c = np.zeros(deg + 1, complex)
    for (i, j), q in p.items():
        if var == 0:
            if j:
                continue
            # (1/2 + t)^i
            for k in range(i + 1):
                import math
                c[k] += float(q) * math.comb(i, k) * 0.5 ** (i - k)
        else:
            c[j] += float(q) * 0.5 ** i
    c = c[::-1]
    while len(c) > 1 and abs(c[0]) < 1e-15:
        c = c[1:]
    return np.roots(c) if len(c) > 1 else np.array([])

print("\nnearest singularities (|root|), excluding y = 0 itself:")
for var, lab in ((0, "in t = x - 1/2 at y = 0"), (1, "in y at x = 1/2")):
    best = []
    for i, name, p in curves:
        for r in roots_1d(p, var):
            if abs(r) > 1e-12:
                best.append((abs(r), i, name, r))
    best.sort()
    rad = best[0][0]
    who = sorted({"l%d(%s)" % (i, n) for R, i, n, _ in best if abs(R - rad) < 1e-6})
    print("  %-22s radius %.4f  from %s" % (lab, rad, ", ".join(who)))
