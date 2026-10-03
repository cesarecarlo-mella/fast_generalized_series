"""precision map of the direct DE solver (corner-1 start) vs the exact GPL solution, weights 1-2."""
import sys, pickle; sys.path.insert(0, '/mnt/user-data/outputs/compare_results')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plot_w6_paper as P

res = pickle.load(open('pg_A.pkl', 'rb')) + pickle.load(open('pg_B.pkl', 'rb'))


def dig(a, b, chop=1e-4):
    d = np.abs(a - b); r = np.abs(b)
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.minimum(16.0, -np.log10(d / r + 1e-300))
    g = np.where(d < 1e-300, 16.0, g)
    return g[r >= chop]

rows = []
for p in res:
    g1, g2 = dig(p['J1'], p['ref1']), dig(p['J2'], p['ref2'])
    rows.append([p['x'], p['y'], 1 - p['x'] - p['y'], 1, g1.min(), np.median(g1), g2.min(), np.median(g2)])
d = np.array(rows)
print("points", len(d))
for c, name in ((4, 'w1 min'), (5, 'w1 median'), (6, 'w2 min'), (7, 'w2 median')):
    v = d[:, c]
    print("%-10s worst %5.1f  10%% %5.1f  median %5.1f  best %5.1f" % (name, v.min(), np.percentile(v, 10), np.median(v), v.max()))
pan = np.array([p['panels'] for p in res]); print("panels used: median %d, max %d" % (np.median(pan), pan.max()))
np.save('prec_w2_grid.npy', d)

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})
fig, ax = plt.subplots(2, 2, figsize=(4.4, 4.3), gridspec_kw=dict(wspace=0.12, hspace=0.22))
for r, (w, c0) in enumerate(((1, 4), (2, 6))):
    for c, ttl in enumerate(("min over components", "median over components")):
        P.panel(ax[r, c], d, c0 + c, 'c1only', show_x=(r == 1), show_y=(c == 0), soft=False)
        ax[r, c].plot(0.1, 0.07, marker='x', ms=4, color=P.INK, mew=0.8, zorder=7)
        if r == 0:
            ax[r, c].set_title(ttl, color=P.INK, pad=3, fontsize=9)
    ax[r, 0].annotate("weight %d" % w, xy=(0, 0.5), xycoords="axes fraction", xytext=(-28, 0),
                      textcoords="offset points", rotation=90, ha="center", va="center", color=P.INK, fontsize=9)
P.colorbar(fig, ax)
fig.savefig('direct_de_precision_w2.pdf', bbox_inches='tight'); fig.savefig('direct_de_precision_w2.png', dpi=160, bbox_inches='tight')
print("wrote direct_de_precision_w2.{pdf,png}")
