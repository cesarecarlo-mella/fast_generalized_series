#!/usr/bin/env python3
"""
plot_grid.py -- one representation (native or Bernoulli), one mode, a grid of
panels: rows = MIN / MEDIAN over components, columns = truncation orders
(vs exact: ordp == ordb == N) or order pairs (self-convergence: N -> N+1).
Same style as plot_w6_paper.py (flat cells, no interpolation, cut as in the
compare files).

    python3 plot_grid.py --indir results_w2_consistent --weight 2 --mode c1only \
        --rep native --pairs 12,12 16,16 20,20 --outdir figures/paper
    python3 plot_grid.py --indir results_w6_fast --weight 6 --mode c1only \
        --rep bernoulli --pairs 18,19 19,20 --outdir figures/paper
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_w6_paper as P

ap = argparse.ArgumentParser()
ap.add_argument("--indir", required=True)
ap.add_argument("--outdir", default="figures/paper")
ap.add_argument("--weight", type=int, required=True)
ap.add_argument("--mode", default="c1only")
ap.add_argument("--rep", choices=["native", "bernoulli"], required=True)
ap.add_argument("--pairs", nargs="+", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

pairs = [tuple(int(v) for v in p.split(",")) for p in a.pairs]
exact = all(L == H for L, H in pairs)
col0 = 4 if a.rep == "native" else 6                  # min; median = col0 + 1
ds = [P.load(os.path.join(a.indir, "compare_%s_w%d_ordp%d_ordb%d.m" % (a.mode, a.weight, L, H)))
      for L, H in pairs]
n = len(ds)
fig, ax = plt.subplots(2, n, figsize=(1.47 * n, 2.75), squeeze=False,
                       gridspec_kw=dict(wspace=0.12, hspace=0.18))
for c, ((L, H), d) in enumerate(zip(pairs, ds)):
    for r, stat in enumerate(("min over\ncomponents", "median over\ncomponents")):
        P.panel(ax[r, c], d, col0 + r, a.mode, show_x=(r == 1), show_y=(c == 0))
        if r == 0:
            ax[r, c].set_title((r"$N=%d$" % L) if exact else (r"$N=%d \rightarrow %d$" % (L, H)),
                               color=P.INK, pad=3, fontsize=9)
        if c == 0:
            ax[r, 0].annotate(stat, xy=(0, 0.5), xycoords="axes fraction", xytext=(-32, 0),
                              textcoords="offset points", rotation=90,
                              ha="center", va="center", color=P.INK, fontsize=9)
P.colorbar(fig, ax)
rep = r"native $(x,y,z)$" if a.rep == "native" else "Bernoulli"
what = "truncation order vs exact" if exact else "self-convergence"
# no title: the caption carries it
tag = "exact" if exact else "selfconv"
base = os.path.join(a.outdir, "w%d_%s_%s_%s_grid" % (a.weight, tag, a.mode, a.rep))
fig.savefig(base + ".pdf", bbox_inches="tight")
print("  wrote", base + ".pdf")
