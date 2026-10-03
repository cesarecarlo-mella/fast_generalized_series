#!/usr/bin/env python3
"""
coverage.py -- whole continued triangle (x, y < 0, z > 1) at weights 1, 2:
every expansion vs the EXACT solution continued (exact_continued.py).

Charts (order N in both variables):
  CC1  : old corner-1 series continued (Log -> log|.| + i pi), in (x, y)   |x|,|y| < 1
  CC2  : u = t v^2 (u << v), in (u, v)                                      |y| < 1, near w = 1
  CC3  : w = t v^2 (w << v), in (w, v)                                      |x| < 1, near u = 1
  CC2e : v = t u^2 (v << u), in (u, v)   edge side                          |y| > 1, near w = 1
  CC3e : v = t w^2 (v << w), in (w, v)   edge side                          |x| > 1, near u = 1
Rule: |x| = w/v, |y| = u/v ; |x|+|y| < 1 -> CC1 ; both < 1 -> CC2 if w >= u else CC3 ;
      |y| < 1 -> CC2 ; |x| < 1 -> CC3 ;
      else (v smallest) CC2e if w >= u else CC3e.
"""
import argparse, os, pickle, sys
import numpy as np
from overlap import load_terms, evaluate, digits

ap = argparse.ArgumentParser()
ap.add_argument('--gt', required=True)
ap.add_argument('--N', type=int, default=12)
ap.add_argument('--c1', default='/mnt/user-data/uploads/series.all/cluster_results/blow_up_form_2_order_20/out')
ap.add_argument('--files', nargs='+', required=True, help='CHART=dir ...')
ap.add_argument('--out', default='coverage')
a = ap.parse_args()
G = pickle.load(open(a.gt, 'rb'))
pts = G['pts']; u, w, v = pts.T
ax_, ay_ = w / v, u / v
# CC1 is a series in (x, y): log(1 - x - y) limits it to |x| + |y| < 1 (u + w < v)
rule = np.where(ax_ + ay_ < 1, 'CC1',
       np.where((ax_ < 1) & (ay_ < 1), np.where(w >= u, 'CC2', 'CC3'),
       np.where(ay_ < 1, 'CC2', np.where(ax_ < 1, 'CC3', np.where(w >= u, 'CC2e', 'CC3e')))))
dirs = dict(f.split('=') for f in a.files)
N = a.N


def trunc(comps, names):
    return [[(c, k) for c, k in t if all(e <= 2 * N for at, e in k if at in names)] for t in comps]


res = {}
for wt in (1, 2):
    ref = G['f%d' % wt]
    vals = {}
    vals['CC1'] = evaluate(trunc(load_terms('%s/phys_c1_w%d_ord20.m' % (a.c1, wt)), ('x', 'y')), 'CC1', pts)
    for ch, d in dirs.items():
        vals[ch] = evaluate(load_terms('%s/phys_%s_w%d_ord%d.m' % (d, ch, wt, N)), ch, pts)
    D = {ch: digits(vals[ch], ref, 1e-4) for ch in vals}
    mn = {ch: np.nanmin(D[ch], axis=0) for ch in D}
    md = {ch: np.nanmedian(D[ch], axis=0) for ch in D}
    sel_min = np.array([mn[r][i] for i, r in enumerate(rule)])
    sel_med = np.array([md[r][i] for i, r in enumerate(rule)])
    best = np.max(np.array([mn[ch] for ch in mn]), axis=0)
    res[wt] = dict(sel_min=sel_min, sel_med=sel_med, best=best, mn=mn, md=md)
    print("\n=== weight %d, order %d, %d points ===" % (wt, N, len(pts)))
    for ch in ('CC1', 'CC2', 'CC3', 'CC2e', 'CC3e'):
        m = rule == ch
        if m.any():
            print("  %-4s %3d pts: min digits  worst %5.1f  median %5.1f  best %5.1f | median-over-components: median %5.1f"
                  % (ch, m.sum(), sel_min[m].min(), np.median(sel_min[m]), sel_min[m].max(), np.median(sel_med[m])))
    for k in (2, 4, 6, 8):
        print("  points with MIN digits >= %d : rule %3d/%d   best-chart %3d/%d" % (k, (sel_min >= k).sum(), len(pts), (best >= k).sum(), len(pts)))
pickle.dump(dict(pts=pts, rule=rule, res=res), open(a.out + '.pkl', 'wb'))

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
cmap = plt.get_cmap('Spectral'); bounds = [-2, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16.01]
norm = BoundaryNorm(bounds, cmap.N)
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.6))
for r, wt in enumerate((1, 2)):
    for c, (key, ttl) in enumerate((('sel_min', 'min over components'), ('sel_med', 'median over components'))):
        ax = axs[r, c]
        sc = ax.scatter(u, w, c=np.clip(res[wt][key], -2, 16), cmap=cmap, norm=norm, marker='s', s=40, linewidths=0)
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color='0.4', lw=0.6)
        ax.plot([0, 0.5], [0.5, 0], 'k:', lw=0.5); ax.plot([0, 1/3, 1], [1, 1/3, 0], 'k:', lw=0.5); ax.plot([0.5, 1/3, 0], [0, 1/3, 0.5], 'k:', lw=0.5); ax.plot([1/3, 0.5], [1/3, 0.5], 'k:', lw=0.5)
        ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('$u$'); ax.set_ylabel('$w$')
        ax.set_title('weight %d, %s' % (wt, ttl), fontsize=9)
fig.colorbar(sc, ax=axs, ticks=[1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16], label='correct digits vs continued exact')
fig.savefig(a.out + '.png', dpi=130, bbox_inches='tight')
print("wrote", a.out + '.png')
