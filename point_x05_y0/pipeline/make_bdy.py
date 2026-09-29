#!/usr/bin/env python3
"""
make_bdy.py -- boundary constants at (x, y) = (1/2, 0) from the solutions on
the line y = 0:

    eps_sep/w<w>/sol<w>_1.m   J_w on y = 0 (regular part, Log[y]^0), as HPLs
                              G[a1, ..., an, t] with t = x, a_i in {0, 1}

bdy_c0_w<w>.m = sol<w>_1.m at t = 1/2. The values G[a..., 1/2] are kept EXACT,
as symbolic atoms "G[a1, ..., an, 1/2]" (together with Pi, Zeta[n]); the
recursion carries every constant monomial as a label, so nothing has to be
evaluated here. gpl_num.py gives their numerical values (40 digits); a
reduction to a basis (Log[2], Zeta, PolyLog[n, 1/2], ...) with PolyLogTools
is optional and only changes how the constants look.

    python3 make_bdy.py --sols <.../H/solutions/eps_sep> [--wmax 6] [--out inputs]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mparse as mp

ap = argparse.ArgumentParser()
ap.add_argument('--sols', required=True, help='eps_sep dir with w<w>/sol<w>_1.m')
ap.add_argument('--wmax', type=int, default=6)
ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inputs'))
ap.add_argument('--x0', default='1/2', help='the point on y = 0 (default 1/2)')
a = ap.parse_args()

GT = re.compile(r'^G\[(.*),\s*t\]$')
for w in range(1, a.wmax + 1):
    f = os.path.join(a.sols, 'w%d' % w, 'sol%d_1.m' % w)
    vec = mp.parse_file(f)
    out, monos, gs = [], set(), set()
    for p in vec:
        q = {}
        for k, v in p.items():
            nk = []
            for at, e in k:
                m = GT.match(at)
                if m:
                    at = 'G[%s, %s]' % (m.group(1), a.x0)
                    gs.add(at)
                elif at == 't' or re.search(r'[\[, ]t[\],]', at):
                    sys.exit("%s: unexpected t-dependence in atom %r" % (f, at))
                nk.append((at, e))
            nk = tuple(sorted(nk))
            q[nk] = q.get(nk, 0) + v
            monos.add(nk)
        out.append({k: v for k, v in q.items() if v})
    path = os.path.join(a.out, 'bdy_c0_w%d.m' % w)
    with open(path, 'w') as fh:
        fh.write("{" + ",\n ".join(mp.fmt_poly(c) for c in out) + "}\n")
    print("weight %d: %d components, %d constant monomials, %d distinct G[..., %s] -> %s"
          % (w, len(out), len(monos), len(gs), a.x0, path))
