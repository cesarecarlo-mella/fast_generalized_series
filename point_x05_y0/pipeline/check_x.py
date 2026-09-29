#!/usr/bin/env python3
"""
check_x.py -- end-to-end check along y = 0: the y^0 Log[y]^0 part of our
series J_w(t, y) (from rec_c0_w<w>.pkl) must equal the known line solution
sol<w>_1.m at x = 1/2 + t.  Everything is evaluated numerically (gpl_num.py,
40 digits) at a few t; the difference should be the truncation error
~ (2|t|)^(ORDt+1).

    python3 check_x.py --workdir dist_10_10 --sols <.../eps_sep> -w 1,2,3,4 [--t 0.05,-0.05,0.1]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mpm
import mparse as mp
import pipe_common as pc
import gpl_num as gn

ap = argparse.ArgumentParser()
ap.add_argument('--workdir', required=True)
ap.add_argument('--sols', required=True)
ap.add_argument('-w', '--weights', default='1,2,3,4')
ap.add_argument('--t', default='0.05,-0.05,0.1')
a = ap.parse_args()

C = pc.load_pickle(pc.setup_path(a.workdir, 0))['C']
cache = {}


def mono(key, t=None):
    """numerical value of a constant monomial; G[..., t] atoms at x = 1/2 + t."""
    v = mpm.mpf(1)
    for at, e in key:
        if at.startswith('G[') and at.endswith(', t]'):
            word = tuple(int(s) for s in at[2:-4].split(','))
            val = gn.value(word, mpm.mpf(1) / 2 + t)
        else:
            val = cache.get(at)
            if val is None:
                val = cache[at] = gn.atom_value(at)
        v *= val ** (mpm.mpf(e) / 2)
    return v


ok = True
for w in [int(s) for s in a.weights.split(',')]:
    slots = pc.load_pickle(pc.rec_path(a.workdir, 0, w))
    ref = mp.parse_file(os.path.join(a.sols, 'w%d' % w, 'sol%d_1.m' % w))
    # y^0 Log[y]^0 slice:  tag = (pt, py=0, la, lb=0, ck), iy = 1 (y^0)
    part = {}
    for (tag, j, it, iy), q in slots.items():
        pt, py, la, lb, ck = tag
        if py == 0 and lb == 0 and iy == 1:
            assert pt == 0 and la == 0
            part.setdefault(j, []).append((it - 1, ck, q))
    worst = 0
    for tv in [mpm.mpf(s) for s in a.t.split(',')]:
        for j in range(C.ncomp):
            mine = sum((mpm.mpf(q.numerator) / q.denominator) * tv ** e * mono(ck)
                       for e, ck, q in part.get(j, []))
            r = sum((mpm.mpf(q.numerator) / q.denominator) * mono(k, tv) for k, q in ref[j].items())
            d = abs(mine - r) / max(1, abs(r))
            worst = max(worst, d)
    print("weight %d: max |series - sol%d_1| (relative, over %d components, t = %s) = %.1e"
          % (w, w, C.ncomp, a.t, worst), flush=True)
    ok &= worst < 1e-6
print("OK" if ok else "!! MISMATCH")
sys.exit(0 if ok else 1)
