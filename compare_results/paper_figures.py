#!/usr/bin/env python3
"""
paper_figures.py -- the three paper figures, typeset with LaTeX (Computer Modern, the
paper's font) at their final printed size, so that \\includegraphics WITHOUT a width
option reproduces the font size of the text exactly.

    python3 paper_figures.py --w2 results_w2_consistent --w6 results_w6_fast --outdir figures/paper \\
        [--width 6.0] [--fontsize 11]

--width    : printed width of the 3-column figures in inches (default 6.0 = the width at
             which the figure is included; get \\textwidth from LaTeX with \\the\\textwidth,
             72.27 pt = 1 in). The 2x2 figure is 0.72 of it.
--fontsize : font size in pt of labels (default = the paper's body size); tick labels are
             1 pt smaller.
Include with  \\includegraphics{figures/<name>.pdf}   (no width=...), or with
width=<the same --width>.

Figures (no titles: the caption carries them):
  w2_exact_orders_c1only_min.pdf      corner 1, min digits vs exact; rows native / Bernoulli
  w2_exact_coverage_bernoulli_grid.pdf Bernoulli coverage vs exact; rows min / median
  w6_selfconv_coverage_summary.pdf     w6: Bernoulli self-convergence 19->20 min | median,
                                       Delta min digits native | Bernoulli
Data handling is exactly that of plot_w6_paper.py (min panels: one flat cell per grid
point; median panels: slightly smoothed; cut |ref| >= 1e-4 from the compare files).
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon

import plot_w6_paper as P

ap = argparse.ArgumentParser()
ap.add_argument("--w2", required=True, help="dir with compare_{c1only,coverage}_w2_ordp{12,16,20}_*.m")
ap.add_argument("--w6", required=True, help="dir with compare_coverage_w6_ordp1{8,9}_*.m")
ap.add_argument("--outdir", default="figures/paper")
ap.add_argument("--width", type=float, default=6.0)
ap.add_argument("--fontsize", type=float, default=11)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
FS = a.fontsize
INK = "#000000"
GRID = "#6f6f6f"

plt.rcParams.update({
    # Latin Modern = the OpenType version of Computer Modern (the LaTeX default font);
    # maths in matplotlib's Computer Modern set -> matches the paper's text and formulae
    "font.family": "serif", "font.serif": ["Latin Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.unicode_minus": False,
    "font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS,
    "xtick.labelsize": FS - 1, "ytick.labelsize": FS - 1,
    "axes.edgecolor": INK, "axes.linewidth": 0.5, "xtick.color": INK, "ytick.color": INK,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "xtick.direction": "out", "ytick.direction": "out",
    "savefig.dpi": 400, "pdf.fonttype": 42})
P.CLIP = True
TRI = [(0, 0), (1, 0), (0, 1)]


def triangle_axes(ax, show_x, show_y, mode, soft):
    """frame, ticks and decorations of one panel."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.plot([1, 0], [0, 1], color=INK, lw=0.5, zorder=4, solid_capstyle="butt")
    if mode == "coverage":
        c = (1 / 3, 1 / 3)
        for e in ((0.5, 0.5), (0.5, 0.0), (0.0, 0.5)):
            ax.plot([c[0], e[0]], [c[1], e[1]], color="black", lw=0.45, ls=(0, (2.5, 2)),
                    alpha=0.55, zorder=5)
    elif mode in P.EXPANSION:
        (px, py), _ = P.EXPANSION[mode]
        ax.plot(px, py, marker="o", ms=4, mfc="white", mec=INK, mew=0.6, zorder=6, clip_on=False)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xticklabels([r"$0$", r"$0.5$", r"$1$"] if show_x else [])
    ax.set_yticklabels([r"$0$", r"$0.5$", r"$1$"] if show_y else [])
    if show_x:
        ax.set_xlabel(r"$x$", labelpad=1.5)
    if show_y:
        ax.set_ylabel(r"$y$", labelpad=1.5)


def digits_panel(ax, d, col, mode, show_x, show_y):
    v = np.clip(d[:, col], -2, 16)
    soft = col in (5, 7)                       # medians slightly smoothed, minima exact cells
    P.fill(ax, d, v, P.CMAP, P.NORM, P.BOUNDS, soft=soft)
    triangle_axes(ax, show_x, show_y, mode, soft)


def digits_colorbar(fig, cax, ticks=(1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16)):
    sm = plt.cm.ScalarMappable(norm=P.NORM, cmap=P.CMAP)
    cb = fig.colorbar(sm, cax=cax, ticks=list(ticks))
    cb.ax.set_yticklabels([r"$%d$" % t for t in ticks])
    cb.set_label(r"correct digits", labelpad=4)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2, labelsize=FS - 1)
    return cb


def grid_figure(nrow, ncol, width, colorbars=1, vgap=0.19):
    """square panels, fixed gaps (in inches), colour bar(s) on the right."""
    left, right, gap, cbw, cbgap, top, bottom = 0.92, 0.62, 0.20, 0.11, 0.14, 0.26, 0.42
    pw = (width - left - right - cbw - cbgap - gap * (ncol - 1)) / ncol
    height = top + bottom + nrow * pw + vgap * (nrow - 1)
    fig = plt.figure(figsize=(width, height))
    axes = np.empty((nrow, ncol), object)
    for r in range(nrow):
        for c in range(ncol):
            x0 = left + c * (pw + gap)
            y0 = height - top - (r + 1) * pw - r * vgap
            axes[r, c] = fig.add_axes([x0 / width, y0 / height, pw / width, pw / height])
    xcb = (left + ncol * pw + (ncol - 1) * gap + cbgap) / width
    if colorbars == 1:
        y0 = (height - top - nrow * pw - (nrow - 1) * vgap) / height
        cax = [fig.add_axes([xcb, y0, cbw / width, (nrow * pw + (nrow - 1) * vgap) / height])]
    else:
        cax = [fig.add_axes([xcb, (height - top - (r + 1) * pw - r * vgap) / height,
                             cbw / width, pw / height]) for r in range(nrow)]
    return fig, axes, cax


def row_label(ax, text):
    ax.annotate(text, xy=(0, 0.5), xycoords="axes fraction", xytext=(-3.7 * FS, 0),
                textcoords="offset points", rotation=90, ha="center", va="center")


# -------------------------------------------------------------------- figure 1
def fig1():
    ds = [P.load(os.path.join(a.w2, "compare_c1only_w2_ordp%d_ordb%d.m" % (N, N))) for N in (12, 16, 20)]
    fig, ax, cax = grid_figure(2, 3, a.width)
    for c, (N, d) in enumerate(zip((12, 16, 20), ds)):
        for r, col in enumerate((4, 6)):
            digits_panel(ax[r, c], d, col, "c1only", show_x=(r == 1), show_y=(c == 0))
        ax[0, c].set_title(r"$N=%d$" % N, pad=5)
    row_label(ax[0, 0], r"native $(x,y,z)$")
    row_label(ax[1, 0], r"Bernoulli")
    digits_colorbar(fig, cax[0])
    return fig, "w2_exact_orders_c1only_min"


# -------------------------------------------------------------------- figure 2
def fig2():
    ds = [P.load(os.path.join(a.w2, "compare_coverage_w2_ordp%d_ordb%d.m" % (N, N))) for N in (12, 16, 20)]
    fig, ax, cax = grid_figure(2, 3, a.width)
    for c, (N, d) in enumerate(zip((12, 16, 20), ds)):
        for r, col in enumerate((6, 7)):
            digits_panel(ax[r, c], d, col, "coverage", show_x=(r == 1), show_y=(c == 0))
        ax[0, c].set_title(r"$N=%d$" % N, pad=5)
    row_label(ax[0, 0], r"min")
    row_label(ax[1, 0], r"median")
    digits_colorbar(fig, cax[0])
    return fig, "w2_exact_coverage_bernoulli_grid"


# -------------------------------------------------------------------- figure 3
def fig3():
    f = lambda L, H: os.path.join(a.w6, "compare_coverage_w6_ordp%d_ordb%d.m" % (L, H))
    d1, d2 = P.load(f(18, 19)), P.load(f(19, 20))
    assert d1.shape == d2.shape and np.allclose(d1[:, :2], d2[:, :2])
    fig, ax, cax = grid_figure(2, 2, 0.72 * a.width, colorbars=2, vgap=0.62)
    # top: Bernoulli self-convergence 19 -> 20
    for c, (col, t) in enumerate(((6, r"min"), (7, r"median"))):
        digits_panel(ax[0, c], d2, col, "coverage", show_x=True, show_y=(c == 0))
        ax[0, c].set_title(t, pad=5)
    row_label(ax[0, 0], r"Bernoulli, $N=19\to 20$")
    cb = digits_colorbar(fig, cax[0], ticks=(2, 4, 6, 8, 12, 16))
    cb.set_label(r"correct digits", labelpad=4)
    # bottom: Delta min digits, linearly interpolated, sequential scale 0 .. 1.6
    import matplotlib.tri as mtri
    cmap, norm = matplotlib.colormaps["YlGnBu"], Normalize(0.0, 1.6)
    for c, (col, t) in enumerate(((4, r"native $(x,y,z)$"), (6, r"Bernoulli"))):
        dv = np.nan_to_num(np.clip(d2[:, col], -2, 16) - np.clip(d1[:, col], -2, 16))
        assert (dv > -0.1).all(), "a point got worse: a sequential scale would hide it"
        interp = mtri.LinearTriInterpolator(P.triangulation(d2), np.clip(dv, 0, 1.6))
        g = np.linspace(0, 1, 700)
        X, Y = np.meshgrid(g, g)
        Z = interp(X, Y)
        Z = np.ma.masked_where((X + Y > 1) | np.ma.getmaskarray(Z), Z)
        im = ax[1, c].imshow(Z, origin="lower", extent=(0, 1, 0, 1), cmap=cmap, norm=norm,
                             interpolation="bilinear", rasterized=True, zorder=1)
        im.set_clip_path(Polygon(TRI, closed=True, transform=ax[1, c].transData))
        triangle_axes(ax[1, c], True, c == 0, "coverage", True)
        ax[1, c].set_title(t, pad=5)
    row_label(ax[1, 0], r"$\Delta$ min digits")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=cax[1], ticks=[0, 0.4, 0.8, 1.2, 1.6])
    cb.ax.set_yticklabels([r"$0$", r"$0.4$", r"$0.8$", r"$1.2$", r"$1.6$"])
    cb.set_label(r"$\Delta$ digits", labelpad=4)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2, labelsize=FS - 1)
    return fig, "w6_selfconv_coverage_summary"


for make in (fig1, fig2, fig3):
    fig, name = make()
    path = os.path.join(a.outdir, name + ".pdf")
    fig.savefig(path)                     # NO bbox_inches="tight": the size must stay exact
    w, h = fig.get_size_inches()
    plt.close(fig)
    print("  wrote %s  (%.2f x %.2f in)" % (path, w, h))
