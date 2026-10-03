#!/usr/bin/env python3
"""
cc1_from_c1.py -- the CC1 expansion (continued region, near x = y = 0 with
x, y < 0) from the OLD corner-1 physical series, by exact analytic
continuation with the convention used for all CC boundaries:

    x = -X, y = -Y  (X, Y > 0;  X = w/v, Y = u/v)
    Log[x] -> Log[X] + I Pi ,  Log[y] -> Log[Y] + I Pi
    x^(n/2) -> I^n X^(n/2)    (so Sqrt[x] Sqrt[y] Sqrt[z] = -Sqrt|xyz|)

    python3 cc1_from_c1.py phys_c1_w2_ord20.m phys_CC1_w2_ord20.m [--ord N]

--ord N keeps only X^e, Y^e with e <= N (default: everything in the file).
"""
import argparse
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline'))
import mparse as mp
from phys import iter_top_list

ap = argparse.ArgumentParser()
ap.add_argument('src'); ap.add_argument('dst')
ap.add_argument('--ord', type=int, default=None)
a = ap.parse_args()



def ipow(n):
    """I^n as a Poly (I^2 = -1 reduced)."""
    sign = (-1) ** ((n % 4) // 2)
    return {(('I', 2),): sign} if n % 2 else {(): sign}


def convert(p):
    out = {}
    for key, c in p.items():
        d = dict(key)
        ex, ey = d.pop('x', 0), d.pop('y', 0)
        if a.ord is not None and (ex > 2 * a.ord or ey > 2 * a.ord):
            continue
        kx, ky = d.pop('Log[x]', 0) // 2, d.pop('Log[y]', 0) // 2
        base = {tuple(sorted(d.items())): c}
        # phases: x^(ex/2) = I^ex X^(ex/2), same for y
        mono = []
        if ex: mono.append(('X', ex))
        if ey: mono.append(('Y', ey))
        term = mp.pmul(mp.pmul(base, {tuple(sorted(mono)): 1}), ipow(ex + ey))   # x^(ex/2) -> I^ex X^(ex/2)
        for kk, L in ((kx, 'Log[X]'), (ky, 'Log[Y]')):     # (Log + I Pi)^k
            if kk:
                s = {}
                for j in range(kk + 1):
                    piece = {tuple(sorted(([(L, 2 * j)] if j else []) + ([('Pi', 2 * (kk - j))] if kk - j else []))): math.comb(kk, j)}
                    s = mp.padd(s, mp.pmul(piece, ipow(kk - j)))
                term = mp.pmul(term, s)
        out = mp.padd(out, term)
    return out


s = open(a.src).read()
s = re.sub(r'\\\r?\n\s*', '', s)
comps = [convert(mp.parse_string(c.strip())) for c in iter_top_list(s)]
with open(a.dst, 'w') as f:
    f.write("{" + ",\n ".join(mp.fmt_poly(c) for c in comps) + "}\n")
print("wrote", a.dst, "(%d components)" % len(comps))
