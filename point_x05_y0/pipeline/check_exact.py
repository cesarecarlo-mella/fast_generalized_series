#!/usr/bin/env python3
"""
check_exact.py -- compare the series J_w(t, y) around (1/2, 0) with the exact
full solution (H/solutions/full_sol/eps_sep/eps_<w>.m, GPLs G[a..., x] with
letters depending on y and G[a..., y]), numerically at points (x, y) near the
expansion point. Weights 1 and 2 (the GPLs are evaluated by direct numerical
integration, mpmath).

    python3 check_exact.py --workdir dist_10_10 --full <.../full_sol/eps_sep> -w 1,2
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mpm
import mparse as mp
import pipe_common as pc
import gpl_num as gn

mpm.mp.dps = 30
ap = argparse.ArgumentParser()
ap.add_argument('--workdir', required=True)
ap.add_argument('--full', required=True)
ap.add_argument('-w', '--weights', default='1,2')
ap.add_argument('--points', default='0.53,0.02;0.47,0.03;0.55,0.05;0.45,0.01')
a = ap.parse_args()


def G(args, z):
    """G(a1, ..., an; z) numerically, n <= 2 (n = 3 also works, slowly)."""
    if all(ai == 0 for ai in args):
        return mpm.log(z) ** len(args) / mpm.factorial(len(args))
    if len(args) == 1:
        return mpm.log(1 - z / args[0])
    a1, rest = args[0], args[1:]
    return mpm.quad(lambda s: G(rest, s) / (s - a1), [0, z])


GATOM = re.compile(r'^G\[(.*)\]$')


def atom_val(at, x, y, cache):
    v = cache.get(at)
    if v is not None:
        return v
    m = GATOM.match(at)
    if m and m.group(1).split(',')[-1].strip() in ('x', 'y'):
        parts = [p.strip() for p in m.group(1).split(',')]
        env = {'x': x, 'y': y}
        args = [mpm.mpf(eval(p, {}, env)) for p in parts[:-1]]
        v = G(args, env[parts[-1]])
    else:
        v = gn.atom_value(at)
    cache[at] = v
    return v


def series_value(slots, j, t, y, cache):
    s = mpm.mpc(0)
    lt, ly = (mpm.log(t) if t > 0 else None), mpm.log(y)
    for (tag, jj, it, iy), q in slots.get(j, []):
        pt, py, la, lb, ck = tag
        term = mpm.mpf(q.numerator) / q.denominator
        term *= mpm.mpf(t) ** (it - 1 + mpm.mpf(pt) / 2) if (it - 1 or pt) else 1
        term *= mpm.mpf(y) ** (iy - 1 + mpm.mpf(py) / 2)
        if la:
            term *= lt ** la
        if lb:
            term *= ly ** lb
        for at, e in ck:
            term *= atom_val(at, None, None, cache) ** (mpm.mpf(e) / 2)
        s += term
    return s


ok = True
C = pc.load_pickle(pc.setup_path(a.workdir, 0))['C']
pts = [tuple(mpm.mpf(v) for v in p.split(',')) for p in a.points.split(';')]
for w in [int(s) for s in a.weights.split(',')]:
    rec = pc.load_pickle(pc.rec_path(a.workdir, 0, w))
    byc = {}
    for k, q in rec.items():
        byc.setdefault(k[1], []).append((k, q))
    ref = mp.parse_file(os.path.join(a.full, 'eps_%d.m' % w))
    worst = 0
    for (x, y) in pts:
        cache, ccache = {}, {}
        for j in range(C.ncomp):
            r = mpm.mpc(0)
            for k, q in ref[j].items():
                term = mpm.mpf(q.numerator) / q.denominator
                for at, e in k:
                    term *= atom_val(at, x, y, cache) ** (mpm.mpf(e) / 2)
                r += term
            mine = series_value(byc, j, x - mpm.mpf(1) / 2, y, ccache)
            d = abs(mine - r) / max(1, abs(r))
            worst = max(worst, d)
        print("  weight %d  (x, y) = (%s, %s): max rel. diff so far %.1e"
              % (w, mpm.nstr(x, 3), mpm.nstr(y, 3), worst), flush=True)
    print("weight %d: max |series - exact| (relative, %d components, %d points) = %.1e"
          % (w, C.ncomp, len(pts), worst), flush=True)
    ok &= worst < 1e-6
print("OK" if ok else "!! MISMATCH")
sys.exit(0 if ok else 1)
