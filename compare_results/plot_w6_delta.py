#!/usr/bin/env python3
"""
plot_w6_delta.py -- change in correct digits between two order pairs, e.g.
    delta = digits(19 -> 20) - digits(18 -> 19)
for one mode and one statistic (default: coverage, MIN over components).
  --stat min|median : change of the MIN / MEDIAN digits  (reads compare_*.m)
  --stat mean       : MEAN over components of each component's own change
                      (reads rawdata_*.m; components with |ref| < chop are dropped,
                      as in 03_analyze.wl, and so are components at 16 digits in both
                      pairs, which cannot change; a point with none left counts as 0)
Blue = more digits (improved), orange = fewer digits (got worse), gray = no change
(|delta| < 0.1 digit, which includes points saturated at 16 digits in both).
Two panels: native (x,y,z) and Bernoulli.  Vector PDF.

    python3 plot_w6_delta.py --indir results_w6_fast --outdir figures/paper
        [--weight 6] [--mode coverage] [--stat min|median|mean] [--from 18,19] [--to 19,20]
"""
import argparse
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

import plot_w6_paper as P          # style, loader, triangulation

# diverging: orange (worse) <- gray -> blue (better)
DRAMP = ["#a3490f", "#d9782b", "#f2b27a", "#e3e3e3", "#9ec5f4", "#3987e5", "#184f95"]
DBOUNDS = [-8, -1, -0.5, -0.1, 0.1, 0.5, 1, 8]
DCMAP, DNORM = ListedColormap(DRAMP), BoundaryNorm(DBOUNDS, len(DRAMP))


def load_raw(path):
    """rawdata_*.m -> list of (x, y, diffN[], diffB[], ref[])"""
    s = open(path).read()
    s = s[s.index('"data" -> ') + 10:s.rindex('|>')]
    s = re.sub(r'`[0-9.]*', '', s).replace('*^', 'e').replace('{', '[').replace('}', ']')
    return eval(s.replace('Indeterminate', 'float("nan")'))


def digits(d, r):                   # same formula as fast_selfconv.digit_of, vectorised
    d = np.abs(d)
    g = np.minimum(16.0, -np.log10(d / np.abs(r) + 1e-30))
    return np.where(d < 1e-16, 16.0, g)


def mean_delta(raw1, raw2, chop):
    """per point: (x, y, mean dNative, mean dBernoulli, frac worse N, frac worse B)"""
    out = []
    for p1, p2 in zip(raw1, raw2):
        assert abs(p1[0] - p2[0]) < 1e-12 and abs(p1[1] - p2[1]) < 1e-12, "grids differ"
        ref = np.array(p2[6], float)
        keep = np.abs(ref) >= chop
        row = [p2[0], p2[1]]
        fr = []
        for k in (4, 5):
            g1 = digits(np.array(p1[k], float)[keep], ref[keep])
            g2 = digits(np.array(p2[k], float)[keep], ref[keep])
            live = ~((g1 >= 16) & (g2 >= 16))
            dv = (g2 - g1)[live]
            row.append(dv.mean() if dv.size else 0.0)
            fr.append(np.mean(dv < -0.1) if dv.size else 0.0)
        out.append(row + fr)
    return np.array(out, float)


def decorate(ax, mode, show_y):
    ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color=P.EDGE, lw=0.6, zorder=3)
    if mode == "coverage":
        c = (1 / 3, 1 / 3)
        for e in ((0.5, 0.5), (0.5, 0.0), (0.0, 0.5)):
            ax.plot([c[0], e[0]], [c[1], e[1]], color=P.INK, lw=0.4, ls=(0, (3, 3)), alpha=0.6, zorder=5)
    elif mode in P.EXPANSION:
        (px, py), _ = P.EXPANSION[mode]
        ax.plot(px, py, marker="o", ms=4.5, mfc="white", mec=P.INK, mew=0.8, zorder=6, clip_on=False)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"]); ax.set_xlabel(r"$x$", labelpad=1)
    ax.set_yticklabels(["0", "0.5", "1"] if show_y else [])
    if show_y: ax.set_ylabel(r"$y$", labelpad=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", default="figures/paper")
    ap.add_argument("--weight", type=int, default=6)
    ap.add_argument("--mode", default="coverage")
    ap.add_argument("--stat", choices=["min", "median", "mean"], default="min")
    ap.add_argument("--chop", type=float, default=1e-4)
    ap.add_argument("--style", choices=["cells", "smooth"], default="cells")
    ap.add_argument("--from", dest="a", default="18,19")
    ap.add_argument("--to", dest="b", default="19,20")
    a = ap.parse_args()
    P.STYLE = a.style
    os.makedirs(a.outdir, exist_ok=True)
    (L1, H1), (L2, H2) = (tuple(int(v) for v in s.split(",")) for s in (a.a, a.b))
    W, m = a.weight, a.mode
    if a.stat == "mean":
        md = mean_delta(load_raw(os.path.join(a.indir, "rawdata_%s_w%d_ordp%d_ordb%d.m" % (m, W, L1, H1))),
                        load_raw(os.path.join(a.indir, "rawdata_%s_w%d_ordp%d_ordb%d.m" % (m, W, L2, H2))),
                        a.chop)
        tri = P.triangulation(md)
        deltas = [md[:, 2], md[:, 3]]
        worse_comp = [md[:, 4], md[:, 5]]
    else:
        d1 = P.load(os.path.join(a.indir, "compare_%s_w%d_ordp%d_ordb%d.m" % (m, W, L1, H1)))
        d2 = P.load(os.path.join(a.indir, "compare_%s_w%d_ordp%d_ordb%d.m" % (m, W, L2, H2)))
        assert d1.shape == d2.shape and np.allclose(d1[:, :2], d2[:, :2]), "grids differ"
        tri = P.triangulation(d2)
        off = 0 if a.stat == "min" else 1
        deltas = [np.nan_to_num(np.clip(d2[:, c], -2, 16) - np.clip(d1[:, c], -2, 16))
                  for c in (4 + off, 6 + off)]
        worse_comp = None

    fig, ax = plt.subplots(1, 2, figsize=(4.9, 2.55), gridspec_kw=dict(wspace=0.12))
    for c, rep in enumerate((r"native $(x,y,z)$", "Bernoulli")):
        dv = deltas[c]
        P.fill(ax[c], np.c_[md[:, :2]] if a.stat == "mean" else d2, np.clip(dv, -7.99, 7.99), DCMAP, DNORM, DBOUNDS)
        decorate(ax[c], m, show_y=(c == 0))
        ax[c].set_title(rep, color=P.INK, pad=4)
        better, worse = np.mean(dv >= 0.1), np.mean(dv <= -0.1)
        ax[c].text(0.98, 0.98, "grid points with\n"
                   "$\\Delta>+0.1$: %d%%\n$|\\Delta|<0.1$: %d%%\n$\\Delta<-0.1$: %d%%\nmedian $\\Delta$ %+.1f"
                   % (round(100 * better), round(100 * (1 - better - worse)), round(100 * worse), np.median(dv)),
                   transform=ax[c].transAxes, ha="right", va="top", fontsize=7, color=P.INK2)
        print("  %-9s improved %.1f%%  worse %.1f%%  median %+.2f  min %+.2f  max %+.2f"
              % (rep.split()[0], 100 * better, 100 * worse, np.median(dv), dv.min(), dv.max()))
        if worse_comp is not None:
            print("            components that got worse (>0.1 digit), averaged over points: %.2f%%"
                  "  (max at one point %.1f%%)" % (100 * worse_comp[c].mean(), 100 * worse_comp[c].max()))

    sm = plt.cm.ScalarMappable(norm=DNORM, cmap=DCMAP)
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.03, ticks=DBOUNDS[1:-1])
    cb.ax.set_yticklabels(["$-1$", "$-0.5$", "$-0.1$", "$+0.1$", "$+0.5$", "$+1$"])
    cb.set_label(r"$\Delta$ digits", color=P.INK)
    cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(P.EDGE)
    cb.ax.tick_params(width=0.5, length=2, colors=P.INK2)
    fig.suptitle(r"weight %d, %s, %s digits: $(N=%d\!\rightarrow\!%d)-(N=%d\!\rightarrow\!%d)$"
                 % (W, P.MODE_TITLE[m], "mean component" if a.stat == "mean" else a.stat, L2, H2, L1, H1), color=P.INK, fontsize=9.5, y=1.02)
    base = os.path.join(a.outdir, "w%d_delta_%s_%s_N%d-%d_vs_N%d-%d" % (W, m, a.stat, L2, H2, L1, H1))
    fig.savefig(base + ".pdf", bbox_inches="tight")
    print("  wrote", base + ".pdf")
