#!/usr/bin/env python3
"""
plot_w6_paper.py -- publication heat maps of the digits in the
compare_<mode>_w<W>_ordp<P>_ordb<B>.m files (fast_selfconv.py / 03_analyze.wl).
Works for any weight (--weight) and for both kinds of comparison:
  * self-convergence   (ordp < ordb: order L vs order H of the same series)
  * vs exact solution  (ordp == ordb == N: truncated series vs exact result)

Each panel is the physical triangle x, y > 0, x + y < 1, coloured by the
estimated number of correct digits (low vs high truncation order); MIN =
worst component, MEDIAN = typical component.  Sequential single-hue scale in
bins of 2 digits (prints in grayscale, colour-blind safe), with 4/8/12-digit
contours.  Vector PDF.

    python3 plot_w6_paper.py --indir results_w6_fast --outdir figures/paper
        [--weight 6] [--pair 19,20] [--modes c1only,c2only,c3only,coverage] [--compare-pairs 18,19 19,20]
    python3 plot_w6_paper.py --indir results_w2 --outdir figures/paper --weight 2 \
        --pair 20,20 --compare-pairs 12,12 16,16 20,20              # vs exact

Figures written (per mode):
    w<W>_selfconv_<mode>_N<L>-<H>  or  w<W>_exact_<mode>_N<N>   (.pdf)
                        2x2: rows native (x,y,z) / Bernoulli, columns MIN / MEDIAN
    w<W>_pairs_<mode>_min  or  w<W>_exact_orders_<mode>_min
                        MIN digits, one column per order (pair), native & Bernoulli
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

# sequential ramp, light -> dark (validated single-hue blue ramp)
RAMP = ["#e8f1fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
# 1-digit bins up to 8 digits (thin dips, e.g. where a component crosses zero, stay
# visible), 2-digit bins above
BOUNDS = [-2, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16.01]
from matplotlib.colors import LinearSegmentedColormap as _LSC
CMAP = ListedColormap(_LSC.from_list("ramp", RAMP)(np.linspace(0, 1, len(BOUNDS) - 1)))
NORM = BoundaryNorm(BOUNDS, CMAP.N)
INK, INK2, EDGE = "#1f1f1f", "#5a5a5a", "#8c8c8c"

EXPANSION = {"c1only": ((0, 0), r"$(x,y)\to 0$"),
             "c2only": ((1, 0), r"$(y,z)\to 0$"),
             "c3only": ((0, 1), r"$(x,z)\to 0$")}
MODE_TITLE = {"c1only": r"corner 1, $(x,y)\to 0$", "c2only": r"corner 2, $(y,z)\to 0$",
              "c3only": r"corner 3, $(x,z)\to 0$", "coverage": "nearest-corner coverage"}

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9,
    "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": EDGE, "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "pdf.fonttype": 42, "savefig.dpi": 300})


def load(path):
    s = open(path).read()
    s = s[s.index('"data" -> ') + 10:s.rindex('|>')]
    s = re.sub(r'`[0-9.]*', '', s).replace('*^', 'e').replace('{', '[').replace('}', ']')
    s = s.replace('Indeterminate', 'float("nan")')
    return np.array(eval(s), float)          # cols: x y z corner dNmin dNmed dBmin dBmed


def triangulation(d):
    x, y = d[:, 0], d[:, 1]
    tri = mtri.Triangulation(x, y)
    step = np.min(np.diff(np.unique(np.round(x, 9))))
    t = tri.triangles
    ex = np.stack([np.hypot(x[t[:, i]] - x[t[:, (i + 1) % 3]], y[t[:, i]] - y[t[:, (i + 1) % 3]])
                   for i in range(3)], 1)
    tri.set_mask(ex.max(1) > 1.6 * step)       # only true lattice triangles
    return tri


STYLE = "cells"
CONTOURS = False    # 4/8/12 contour lines (interpolated) -- off by default     # "cells": one flat cell per grid point (no interpolation); "smooth": tricontourf


def fill(ax, d, v, cmap, norm, levels):
    """colour the triangle by v at the grid points d[:, 0:2].
    cells  : every grid point is its own flat square (exact data, thin features such as
             the lines where a component crosses zero stay visible)
    smooth : filled contours on the triangulation (linear interpolation between points)"""
    if STYLE == "smooth":
        ax.tricontourf(triangulation(d), v, levels=levels, cmap=cmap, norm=norm)
        return
    xs, ys = np.unique(np.round(d[:, 0], 9)), np.unique(np.round(d[:, 1], 9))
    V = np.full((len(ys), len(xs)), np.nan)
    iy, ix = np.searchsorted(ys, np.round(d[:, 1], 9)), np.searchsorted(xs, np.round(d[:, 0], 9))
    assert len(set(zip(iy.tolist(), ix.tolist()))) == len(d), "two grid points fall in one cell"
    V[iy, ix] = v
    # regular lattice check: cell width = grid spacing everywhere (no stretched cells)
    for c in (xs, ys):
        dc = np.diff(c)
        assert np.allclose(dc, dc[0], rtol=1e-6), "grid is not a regular lattice"
    edges = lambda c: np.concatenate(([c[0] - (c[1] - c[0]) / 2], (c[1:] + c[:-1]) / 2, [c[-1] + (c[-1] - c[-2]) / 2]))
    ax.pcolormesh(edges(xs), edges(ys), np.ma.masked_invalid(V), cmap=cmap, norm=norm,
                  shading="flat", rasterized=True)


def panel(ax, d, col, mode, show_x=True, show_y=True):
    tri = triangulation(d)
    v = np.clip(d[:, col], -2, 16)
    fill(ax, d, v, CMAP, NORM, BOUNDS)
    if CONTOURS:      # optional guide lines: these ARE interpolated (linear, on the triangulation)
        cs = ax.tricontour(tri, v, levels=[4, 8, 12], colors=INK, linewidths=0.45, alpha=0.75)
        ax.clabel(cs, fmt="%d", fontsize=6.5, inline=True, inline_spacing=2)
    ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color=EDGE, lw=0.6, zorder=3)
    if mode == "coverage":            # regions of the nearest-corner switch (largest of x, y, z)
        c = (1 / 3, 1 / 3)
        for e in ((0.5, 0.5), (0.5, 0.0), (0.0, 0.5)):
            # thin dashed line only (no white halo: it would hide the data along the switch)
            ax.plot([c[0], e[0]], [c[1], e[1]], color=INK, lw=0.4, ls=(0, (3, 3)), alpha=0.6, zorder=5)
    elif mode in EXPANSION:
        (px, py), _ = EXPANSION[mode]
        ax.plot(px, py, marker="o", ms=4.5, mfc="white", mec=INK, mew=0.8, zorder=6, clip_on=False)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"] if show_x else [])
    ax.set_yticklabels(["0", "0.5", "1"] if show_y else [])
    if show_x: ax.set_xlabel(r"$x$", labelpad=1)
    if show_y: ax.set_ylabel(r"$y$", labelpad=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def colorbar(fig, axes):
    sm = plt.cm.ScalarMappable(norm=NORM, cmap=CMAP)
    cb = fig.colorbar(sm, ax=axes, fraction=0.035, pad=0.03, ticks=[1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16])
    cb.ax.set_yticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "10", "12", "14", "16"])
    cb.set_label("correct digits", color=INK)
    cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(EDGE)
    cb.ax.tick_params(width=0.5, length=2, colors=INK2)
    return cb


def fig_2x2(indir, outdir, mode, L, H, W):
    f = os.path.join(indir, "compare_%s_w%d_ordp%d_ordb%d.m" % (mode, W, L, H))
    if not os.path.isfile(f):
        print("  skip (missing)", f); return
    d = load(f)
    fig, ax = plt.subplots(2, 2, figsize=(4.9, 4.55), gridspec_kw=dict(wspace=0.12, hspace=0.18))
    for r, (cmin, cmed, rep) in enumerate(((4, 5, r"native $(x,y,z)$"), (6, 7, "Bernoulli"))):
        for c, (col, stat) in enumerate(((cmin, "min over components"), (cmed, "median over components"))):
            panel(ax[r, c], d, col, mode, show_x=(r == 1), show_y=(c == 0))
            if r == 0:
                ax[r, c].set_title(stat, color=INK, pad=4)
        ax[r, 0].annotate(rep, xy=(-0.36, 0.5), xycoords="axes fraction", rotation=90,
                          ha="center", va="center", color=INK, fontsize=9)
    colorbar(fig, ax)
    if L == H:
        title = r"weight %d, %s: truncation $N=%d$ vs exact result" % (W, MODE_TITLE[mode], L)
        base = os.path.join(outdir, "w%d_exact_%s_N%d" % (W, mode, L))
    else:
        title = r"weight %d, %s: agreement of orders $N=%d$ and $N=%d$" % (W, MODE_TITLE[mode], L, H)
        base = os.path.join(outdir, "w%d_selfconv_%s_N%d-%d" % (W, mode, L, H))
    fig.suptitle(title, color=INK, fontsize=9.5, y=0.985)
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig); print("  wrote", base + ".pdf")


def fig_pairs(indir, outdir, mode, pairs, W):
    ds = []
    for L, H in pairs:
        f = os.path.join(indir, "compare_%s_w%d_ordp%d_ordb%d.m" % (mode, W, L, H))
        if os.path.isfile(f):
            ds.append(((L, H), load(f)))
    if len(ds) < 2:
        return
    n = len(ds)
    fig, ax = plt.subplots(2, n, figsize=(2.45 * n, 4.55), squeeze=False,
                           gridspec_kw=dict(wspace=0.12, hspace=0.18))
    for c, ((L, H), d) in enumerate(ds):
        for r, (col, rep) in enumerate(((4, r"native $(x,y,z)$"), (6, "Bernoulli"))):
            panel(ax[r, c], d, col, mode, show_x=(r == 1), show_y=(c == 0))
            if r == 0:
                ax[r, c].set_title((r"$N=%d$" % L) if L == H else (r"$N=%d \rightarrow %d$" % (L, H)),
                                   color=INK, pad=4)
            if c == 0:
                ax[r, 0].annotate(rep, xy=(-0.36, 0.5), xycoords="axes fraction", rotation=90,
                                  ha="center", va="center", color=INK, fontsize=9)
    colorbar(fig, ax)
    exact = all(L == H for (L, H), _ in ds)
    fig.suptitle(r"weight %d, %s: min digits vs truncation order%s"
                 % (W, MODE_TITLE[mode], " (vs exact)" if exact else ""),
                 color=INK, fontsize=9.5, y=0.985)
    base = os.path.join(outdir, ("w%d_exact_orders_%s_min" if exact else "w%d_pairs_%s_min") % (W, mode))
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig); print("  wrote", base + ".pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", default="figures/paper")
    ap.add_argument("--weight", type=int, default=6)
    ap.add_argument("--style", choices=["cells", "smooth"], default="cells")
    ap.add_argument("--contours", action="store_true", help="draw 4/8/12 contour lines (interpolated)")
    ap.add_argument("--pair", default="19,20")
    ap.add_argument("--modes", default="c1only,coverage")
    ap.add_argument("--compare-pairs", nargs="+", default=["18,19", "19,20"])
    a = ap.parse_args()
    STYLE = a.style
    CONTOURS = a.contours
    os.makedirs(a.outdir, exist_ok=True)
    L, H = (int(v) for v in a.pair.split(","))
    pairs = [tuple(int(v) for v in p.split(",")) for p in a.compare_pairs]
    for m in a.modes.split(","):
        fig_2x2(a.indir, a.outdir, m, L, H, a.weight)
        fig_pairs(a.indir, a.outdir, m, pairs, a.weight)
