#!/usr/bin/env python3
"""
coverage2.py -- whole continued triangle, seven expansions, weights 1-2, vs the
EXACT solution continued (exact_continued.py).

  CC1                     X + Y < 1                       (X = |x| = w/v, Y = |y| = u/v)
  CC2 / CC2e / CC2r       w >= u side: Y < 1 / Y > 1 / band around Y = 1
  CC3 / CC3e / CC3r       u >  w side: X < 1 / X > 1 / band around X = 1
Rule: X+Y<1 -> CC1; otherwise on the CC2 side (w >= u) use CC2r if |t| < TMAX and p < PMAX
(its own variables: p = 1/X, t = 1 - Y + p), else CC2 (Y < 1) or CC2e (Y > 1); mirror on the CC3 side.
"""
import argparse, pickle
import numpy as np
from overlap import load_terms, evaluate, digits
from rcharts import eval_exact, tp_of

ap = argparse.ArgumentParser()
ap.add_argument('--gt', required=True)
ap.add_argument('--N', type=int, default=20)
ap.add_argument('--c1', default='/mnt/user-data/uploads/series.all/cluster_results/blow_up_form_2_order_20/out')
ap.add_argument('--blowup', required=True, help='dir with phys_CC{2,3,2e,3e}_w<w>_ord<N>.m')
ap.add_argument('--ratio', required=True, help='dir with exact_CC{2r,3r}_w<w>.m')
ap.add_argument('--tmax', type=float, default=0.5)
ap.add_argument('--pmax', type=float, default=0.35)
ap.add_argument('--out', default='coverage7')
a = ap.parse_args()
G = pickle.load(open(a.gt, 'rb')); pts = G['pts']; u, w, v = pts.T
X, Y = w / v, u / v
N = a.N


def trunc(comps, names):
    return [[(c, k) for c, k in t if all(e <= 2 * N for at, e in k if at in names)] for t in comps]


side2 = w >= u
t2, p2 = tp_of('CC2r', pts); t3, p3 = tp_of('CC3r', pts)
rule = np.empty(len(pts), dtype=object)
for i in range(len(pts)):
    if X[i] + Y[i] < 1:
        rule[i] = 'CC1'
    elif side2[i]:
        rule[i] = 'CC2r' if (abs(t2[i]) < a.tmax and p2[i] < a.pmax) else ('CC2' if Y[i] < 1 else 'CC2e')
    else:
        rule[i] = 'CC3r' if (abs(t3[i]) < a.tmax and p3[i] < a.pmax) else ('CC3' if X[i] < 1 else 'CC3e')
CH = ['CC1', 'CC2', 'CC2e', 'CC2r', 'CC3', 'CC3e', 'CC3r']
res = {}
for wt in (1, 2):
    ref = G['f%d' % wt]
    vals = {'CC1': evaluate(trunc(load_terms('%s/phys_c1_w%d_ord20.m' % (a.c1, wt)), ('x', 'y')), 'CC1', pts)}
    for ch in ('CC2', 'CC3', 'CC2e', 'CC3e'):
        vals[ch] = evaluate(load_terms('%s/phys_%s_w%d_ord%d.m' % (a.blowup, ch, wt, N)), ch, pts)
    for ch in ('CC2r', 'CC3r'):
        vals[ch] = eval_exact('%s/exact_%s_w%d.m' % (a.ratio, ch, wt), ch, pts)
    D = {ch: digits(vals[ch], ref, 1e-4) for ch in CH}
    mn = {ch: np.nanmin(D[ch], axis=0) for ch in CH}
    md = {ch: np.nanmedian(D[ch], axis=0) for ch in CH}
    sel_min = np.array([mn[r][i] for i, r in enumerate(rule)])
    sel_med = np.array([md[r][i] for i, r in enumerate(rule)])
    best = np.max(np.array([mn[ch] for ch in CH]), axis=0)
    res[wt] = dict(sel_min=sel_min, sel_med=sel_med, best=best)
    print("\n=== weight %d, order %d, %d points ===" % (wt, N, len(pts)))
    for ch in CH:
        m = rule == ch
        if m.any():
            print("  %-4s %3d pts: MIN digits worst %5.1f median %5.1f | MEDIAN digits median %5.1f"
                  % (ch, m.sum(), sel_min[m].min(), np.median(sel_min[m]), np.median(sel_med[m])))
    for k in (2, 4, 6, 8):
        print("  MIN digits >= %d : rule %3d/%d   best chart %3d/%d" % (k, (sel_min >= k).sum(), len(pts), (best >= k).sum(), len(pts)))
pickle.dump(dict(pts=pts, rule=rule, res=res), open(a.out + '.pkl', 'wb'))

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
cmap = plt.get_cmap('Spectral'); bounds = [-2, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16.01]
norm = BoundaryNorm(bounds, cmap.N)
fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.6))
cols = {'CC1': 0, 'CC2': 1, 'CC2e': 2, 'CC2r': 3, 'CC3': 4, 'CC3e': 5, 'CC3r': 6}
for r, wt in enumerate((1, 2)):
    for c, (key, ttl) in enumerate((('sel_min', 'min over components'), ('sel_med', 'median over components'))):
        ax = axs[r, c]
        sc = ax.scatter(u, w, c=np.clip(res[wt][key], -2, 16), cmap=cmap, norm=norm, marker='s', s=34, linewidths=0)
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color='0.4', lw=0.6)
        ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel('$u$'); ax.set_ylabel('$w$')
        ax.set_title('weight %d, %s' % (wt, ttl), fontsize=9)
    ax = axs[r, 2]
    if r == 1:
        ax.axis('off'); continue
    ax.scatter(u, w, c=[cols[x] for x in rule], cmap='tab10', vmin=0, vmax=9, marker='s', s=34, linewidths=0)
    ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color='0.4', lw=0.6); ax.set_aspect('equal'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title('chart used' if r == 0 else '', fontsize=9); ax.set_xlabel('$u$')
    for ch, k in cols.items():
        ax.scatter([], [], c=[plt.get_cmap('tab10')(k / 9)], marker='s', label=ch)
    ax.legend(fontsize=7, loc='upper right', frameon=False)
fig.colorbar(sc, ax=axs[:, :2], ticks=[1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16], label='correct digits vs continued exact')
fig.savefig(a.out + '.png', dpi=130, bbox_inches='tight')
print("wrote", a.out + '.png')
