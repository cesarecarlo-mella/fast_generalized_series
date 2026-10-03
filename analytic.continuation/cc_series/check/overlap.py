#!/usr/bin/env python3
"""
overlap.py -- independent check of the continued-region series: evaluate the
SAME 371 masters from different expansions at common points of the
continued region (x, y < 0, z > 1) and count the digits of agreement.

    CC1 : the OLD corner-1 physical series phys_c1_w<w>_ord<N>.m (in x, y),
          analytically continued with
              Log[x] -> log|x| + I Pi ,  Log[y] -> log|y| + I Pi ,
              x^(n/2) = |x|^(n/2) e^(I Pi n/2)  (same for y, Sqrt[x], Sqrt[y])
          -- the continuation used for the CC boundaries, and r = -Sqrt|xyz|
    CC2 : phys_CC2_w<w>_ord<N>.m  in (u, v), chart u = t v^2 (u << v)
    CC3 : phys_CC3_w<w>_ord<N>.m  in (w, v)   (this pipeline)

with x = -w/v, y = -u/v, z = 1/v, u + w + v = 1. The CC2/CC3 series know the
continuation only through their boundary constants (transported from CC1 along
u = 0 / w = 0 and Galois-fixed sign of r), so agreement with the continued
corner-1 series away from those lines tests constants, letters and recursion.

Digits per component = min(16, -log10 |a - b| / |b|), components with |b| < CHOP
(default 1e-4) skipped; MIN and MEDIAN over components per point.

    python3 overlap.py --c1 <dir with phys_c1_w*_ord20.m> --cc <workdir with phys_CC*> \
        --ordc1 20 --ordcc 6 --weights 1,2,3 --pairs CC1:CC2,CC1:CC3,CC2:CC3
"""
import argparse
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))
import mparse as mp
from phys import iter_top_list

ZETA = {3: 1.2020569031595942854, 5: 1.0369277551433699263}
REGVARS = {'CC1': ('x', 'y'), 'CC2': ('u', 'v'), 'CC3': ('w', 'v')}


def load_terms(path):
    """phys file -> list (per component) of [(coef float, key)]"""
    s = open(path).read()
    s = re.sub(r'\\\r?\n\s*', '', s)
    comps = []
    for comp in iter_top_list(s):
        p = mp.parse_string(comp.strip())
        comps.append([(float(v), k) for k, v in p.items()])
    return comps


def atom_values(region, pts):
    """symbol -> complex vector over points (pts: array (P,3) of u, w, v)."""
    u, w, v = pts[:, 0], pts[:, 1], pts[:, 2]
    A = {'Pi': np.full(len(u), math.pi + 0j), 'I': np.full(len(u), 1j)}
    for n, z in ZETA.items():
        A['Zeta[%d]' % n] = np.full(len(u), z + 0j)
    if region == 'CC1':
        for s, val in (('x', w / v), ('y', u / v)):          # |x|, |y|
            A[('Log', s)] = np.log(val) + 1j * math.pi      # Log of the negative number
            A[('abs', s)] = val
    else:
        for s, val in (('u', u), ('w', w), ('v', v)):
            A[('Log', s)] = np.log(val) + 0j
            A[('abs', s)] = val
    return A


def evaluate(comps, region, pts):
    A = atom_values(region, pts)
    cont = region == 'CC1'
    cache = {}

    def atom(at, e2):
        key = (at, e2)
        r = cache.get(key)
        if r is not None:
            return r
        m = re.fullmatch(r'Log\[(\w+)\]', at)
        if m:
            r = A[('Log', m.group(1))] ** (e2 // 2)
        else:
            m = re.fullmatch(r'Sqrt\[(\w+)\]', at)
            if m:
                at, e2 = m.group(1), e2 // 2           # Sqrt[x]^k = x^(k/2): exp2 halves
            if at in ('x', 'y', 'u', 'v', 'w'):
                r = A[('abs', at)] ** (e2 / 2) + 0j
                if cont:                                  # |x|^(e/2) * e^(I Pi e/2)
                    r = r * np.exp(1j * math.pi * e2 / 2)
            elif at == 'I':
                r = A['I'] ** (e2 // 2)
            else:
                r = A[at] ** (e2 // 2)
        cache[key] = r
        return r

    out = np.zeros((len(comps), len(pts)), complex)
    for j, terms in enumerate(comps):
        acc = np.zeros(len(pts), complex)
        for c, key in terms:
            t = c
            for at, e2 in key:
                t = t * atom(at, e2)
            acc = acc + t
        out[j] = acc
    return out


def digits(a, b, chop):
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.minimum(16.0, -np.log10(np.abs(a - b) / np.abs(b) + 1e-300))
    d = np.where(np.abs(a - b) < 1e-300, 16.0, d)
    return np.where(np.abs(b) >= chop, d, np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--c1', required=True)
    ap.add_argument('--cc', required=True)
    ap.add_argument('--ordc1', type=int, default=20)
    ap.add_argument('--ordcc', type=int, required=True)
    ap.add_argument('--weights', default='1,2,3')
    ap.add_argument('--pairs', default='CC1:CC2,CC1:CC3,CC2:CC3')
    ap.add_argument('--chop', type=float, default=1e-4)
    ap.add_argument('--grid', type=int, default=20, help='points per side of the (u,w) triangle')
    ap.add_argument('--out', default=None, help='write per-point results (npz)')
    a = ap.parse_args()

    n = a.grid
    pts = np.array([(i / n, j / n, 1 - (i + j) / n) for i in range(1, n) for j in range(1, n) if i + j < n])
    print("grid: %d interior points of the (u, w) triangle (u, w, v > 0)" % len(pts))
    pairs = [tuple(p.split(':')) for p in a.pairs.split(',')]
    regions = sorted({r for p in pairs for r in p})
    res = {}
    for wgt in [int(s) for s in a.weights.split(',')]:
        vals = {}
        for r in regions:
            f = (os.path.join(a.c1, 'phys_c1_w%d_ord%d.m' % (wgt, a.ordc1)) if r == 'CC1'
                 else os.path.join(a.cc, 'phys_%s_w%d_ord%d.m' % (r, wgt, a.ordcc)))
            t0 = time.time()
            comps = load_terms(f)
            vals[r] = evaluate(comps, r, pts)
            print("  w%d %s: %s  (%d terms, %.0fs)" % (wgt, r, os.path.basename(f), sum(map(len, comps)), time.time() - t0), flush=True)
        for r1, r2 in pairs:
            d = digits(vals[r1], vals[r2], a.chop)
            mn = np.nanmin(d, axis=0); md = np.nanmedian(d, axis=0)
            res[(wgt, r1, r2)] = (mn, md)
            k = np.argsort(-mn)[:5]
            print("  w%d %s vs %s: best MIN digits %s at (u,w,v) %s ; points with MIN >= 4: %d / %d ; max MEDIAN %.1f"
                  % (wgt, r1, r2, np.round(mn[k], 1).tolist(),
                     [tuple(np.round(pts[i], 2)) for i in k[:3]], int(np.sum(mn >= 4)), len(pts), np.nanmax(md)), flush=True)
    if a.out:
        np.savez(a.out, pts=pts, **{'%d_%s_%s_min' % k: v[0] for k, v in res.items()},
                 **{'%d_%s_%s_med' % k: v[1] for k, v in res.items()})


if __name__ == '__main__':
    main()
