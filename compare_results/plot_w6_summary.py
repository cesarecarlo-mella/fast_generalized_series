#!/usr/bin/env python3
"""
plot_w6_summary.py -- one 2x2 figure for weight 6 (nearest-corner coverage):

  top row    : Bernoulli self-convergence N=19 -> 20, MIN (left) and MEDIAN (right)
               over components (digits, same colour scale as the other paper figures)
  bottom row : improvement of the MIN digits from (18 -> 19) to (19 -> 20),
               native (left) and Bernoulli (right)

Flat cells, no interpolation, cut as in the compare files (1e-4).

    python3 plot_w6_summary.py --indir results_w6_fast --outdir figures/paper
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_w6_paper as P
import plot_w6_delta as Dl

ap = argparse.ArgumentParser()
ap.add_argument("--indir", required=True)
ap.add_argument("--outdir", default="figures/paper")
ap.add_argument("--mode", default="coverage")
ap.add_argument("--cells-bottom", dest="smooth_bottom", action="store_false",
                help="bottom row as flat cells (no interpolation) instead of the smooth map")
ap.add_argument("--top-only", action="store_true",
                help="only the top row (Bernoulli 19->20, min and median); writes *_bernoulli.pdf")
a = ap.parse_args()

# smooth bottom row: the per-point values are LINEARLY INTERPOLATED on the triangulation
# of the grid onto a fine 600x600 raster and shown with a continuous diverging colour
# map (same colours as the binned one). Only the picture is smoothed; the percentages
# printed in the panels are computed from the raw grid values.
from matplotlib.colors import Normalize
# no point gets worse (checked below), so the whole colour range is spent on 0 .. +1.6:
# light yellow = no gain, green, dark blue = +1.6 digits per order
SCMAP = matplotlib.colormaps["YlGnBu"]
SLO, SHI = 0.0, 1.6
SNORM = Normalize(vmin=SLO, vmax=SHI)


def smooth_fill(ax, d, v):
    import matplotlib.tri as mtri
    tri = P.triangulation(d)
    assert (v > -0.1).all(), "some points got worse: use --cells-bottom (diverging scale)"
    interp = mtri.LinearTriInterpolator(tri, np.clip(v, SLO, SHI))
    g = np.linspace(0, 1, 600)
    X, Y = np.meshgrid(g, g)
    Z = interp(X, Y)                                   # masked outside the data hull
    Z = np.ma.masked_where((X + Y > 1) | np.ma.getmaskarray(Z), Z)
    ax.imshow(Z, origin="lower", extent=(0, 1, 0, 1), cmap=SCMAP,
              norm=SNORM, interpolation="bilinear",
              rasterized=True, zorder=1)
os.makedirs(a.outdir, exist_ok=True)
m = a.mode
f = lambda L, H: os.path.join(a.indir, "compare_%s_w6_ordp%d_ordb%d.m" % (m, L, H))
d1, d2 = P.load(f(18, 19)), P.load(f(19, 20))
assert d1.shape == d2.shape and np.allclose(d1[:, :2], d2[:, :2]), "grids differ"

if a.top_only:
    fig, ax1 = plt.subplots(1, 2, figsize=(4.4, 2.0), gridspec_kw=dict(wspace=0.14))
    ax = np.array([ax1])
else:
    fig, ax = plt.subplots(2, 2, figsize=(4.4, 4.4), gridspec_kw=dict(wspace=0.14, hspace=0.42))
# --- top: Bernoulli self-convergence 19 -> 20 (digits)
for c, (col, stat) in enumerate(((6, "min over components"), (7, "median over components"))):
    P.panel(ax[0, c], d2, col, m, show_x=True, show_y=(c == 0))
    ax[0, c].set_title(stat, color=P.INK, pad=4, fontsize=9)
ax[0, 0].annotate(r"Bernoulli, $N=19\rightarrow20$", xy=(0, 0.5), xycoords="axes fraction", xytext=(-28, 0),
                  textcoords="offset points", rotation=90, ha="center", va="center", color=P.INK, fontsize=9)
sm = plt.cm.ScalarMappable(norm=P.NORM, cmap=P.CMAP)
cb = fig.colorbar(sm, ax=ax[0, :], fraction=0.035, pad=0.03, ticks=[1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16])
cb.set_label("correct digits", color=P.INK)
cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(P.EDGE)
cb.ax.tick_params(width=0.5, length=2, colors=P.INK2, labelsize=7.5)

if a.top_only:
    base = os.path.join(a.outdir, "w6_selfconv_%s_bernoulli" % m)
    fig.savefig(base + ".pdf", bbox_inches="tight")
    print("  wrote", base + ".pdf")
    raise SystemExit

# --- bottom: change of MIN digits (19->20) - (18->19), native and Bernoulli
for c, (col, rep) in enumerate(((4, r"native $(x,y,z)$"), (6, "Bernoulli"))):
    dv = np.nan_to_num(np.clip(d2[:, col], -2, 16) - np.clip(d1[:, col], -2, 16))
    if a.smooth_bottom:
        smooth_fill(ax[1, c], d2, dv)
    else:
        P.fill(ax[1, c], d2, np.clip(dv, -7.99, 7.99), Dl.DCMAP, Dl.DNORM, Dl.DBOUNDS)
    Dl.decorate(ax[1, c], m, show_y=(c == 0))
    ax[1, c].set_title(rep, color=P.INK, pad=4, fontsize=9)
    # third line: strictly negative changes (none: the script asserts below that dv >= 0)
    better, worse = np.mean(dv >= 0.1), np.mean(dv < 0)
    ax[1, c].text(0.98, 0.98, "$\\Delta>+0.1$: %d%%\n$|\\Delta|<0.1$: %d%%\n$\\Delta<0$: %d%%"
                  % (round(100 * better), round(100 * (1 - better - worse)), round(100 * worse)),
                  transform=ax[1, c].transAxes, ha="right", va="top", fontsize=6.5, color=P.INK2)
    print("  %-9s improved %.1f%%  unchanged %.1f%%  worse %.1f%%  median %+.2f"
          % (rep.split()[0], 100 * better, 100 * (1 - better - worse), 100 * worse, np.median(dv)))
ax[1, 0].annotate(r"$\Delta$ min digits", xy=(0, 0.5), xycoords="axes fraction", xytext=(-28, 0),
                  textcoords="offset points", rotation=90, ha="center", va="center", color=P.INK, fontsize=9)
if a.smooth_bottom:
    sm = plt.cm.ScalarMappable(norm=SNORM, cmap=SCMAP)
    cb = fig.colorbar(sm, ax=ax[1, :], fraction=0.035, pad=0.03, ticks=[0, 0.4, 0.8, 1.2, 1.6])
    cb.ax.set_yticklabels(["$0$", "$+0.4$", "$+0.8$", "$+1.2$", "$+1.6$"])
else:
    sm = plt.cm.ScalarMappable(norm=Dl.DNORM, cmap=Dl.DCMAP)
    cb = fig.colorbar(sm, ax=ax[1, :], fraction=0.035, pad=0.03, ticks=Dl.DBOUNDS[1:-1])
    cb.ax.set_yticklabels(["$-1$", "$-0.5$", "$-0.1$", "$+0.1$", "$+0.5$", "$+1$"])
cb.set_label(r"$(19\!\rightarrow\!20)-(18\!\rightarrow\!19)$", color=P.INK)
cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(P.EDGE)
cb.ax.tick_params(width=0.5, length=2, colors=P.INK2, labelsize=7.5)

# no title: the caption carries it
base = os.path.join(a.outdir, "w6_selfconv_%s_summary" % m)
fig.savefig(base + ".pdf", bbox_inches="tight")
print("  wrote", base + ".pdf")
