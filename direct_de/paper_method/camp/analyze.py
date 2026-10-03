"""Analysis of the campaign (paper method, compiled connection).
  w6: 500 random points per region; double (1e-10, 1e-12) against the dd (1e-14) reference, timings.
  w2: 2016-point lattice per region; double (1e-12) and dd (1e-16) against the exact GPLs.
usage: python3 analyze.py [w6] [w2]"""
import sys, os, pickle, glob
import numpy as np
from decimal import Decimal, getcontext
getcontext().prec = 40
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/mnt/user-data/outputs/compare_results')
import plot_w6_paper as P

H = os.path.dirname(os.path.abspath(__file__))
n = 371
REG = {'e': 'Euclidean region', 's': 'scattering region'}
ACC = {'dbl10': '#3987e5', 'dbl12': '#184f95', 'dd14': '#d95f02', 'dd16': '#d95f02'}
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5, "axes.edgecolor": P.EDGE, "axes.linewidth": 0.6,
                     "xtick.color": P.INK2, "ytick.color": P.INK2})


def hl(s):
    """decimal string -> (hi, lo) float pair (double-double)."""
    d = Decimal(s); hi = float(d); return hi, float(d - Decimal(hi))


def read(tag, hp=False):
    """-> list of dicts (x, y, steps, evals, rej, time, re, im[, re_lo, im_lo])"""
    out = []
    for h in 'AB':
        fn = os.path.join(H, '%s_%s.out' % (tag, h))
        if not os.path.exists(fn): continue
        cur = None
        for line in open(fn):
            if line.startswith('#'):
                p = line.split()
                cur = dict(x=p[1], y=p[2], failed='FAILED' in line)
                if not cur['failed']:
                    cur.update(steps=int(p[3]), evals=int(p[4]), rej=int(p[5]), time=float(p[6]))
                cur['v'] = []
                out.append(cur)
            else:
                a, b = line.split()
                cur['v'].append((hl(a), hl(b)) if hp else (float(a), float(b)))
    for r in out:
        v = r.pop('v')
        if not v: continue
        if hp:
            r['hi'] = np.array([a[0] + 1j * b[0] for a, b in v]); r['lo'] = np.array([a[1] + 1j * b[1] for a, b in v])
        else:
            r['hi'] = np.array([a + 1j * b for a, b in v]); r['lo'] = np.zeros(len(v), complex)
    return [r for r in out if 'hi' in r]


def digits(ahi, alo, bhi, blo, chop=1e-4, cap=32.0):
    """per component -log10 |a-b|/|b| (NaN where |b| < chop); double-double differences."""
    d = (ahi - bhi) + (alo - blo)
    with np.errstate(divide='ignore', invalid='ignore'):
        g = -np.log10(np.abs(d) / np.abs(bhi))
    g = np.where(np.abs(d) == 0, cap, np.minimum(g, cap))
    return np.where(np.abs(bhi) >= chop, g, np.nan)


def tstats(rs):
    t = np.array([r['time'] for r in rs]); e = np.array([r['evals'] for r in rs])
    return dict(n=len(t), median=np.median(t), mean=t.mean(), p10=np.percentile(t, 10), p90=np.percentile(t, 90),
                min=t.min(), max=t.max(), total=t.sum(), evals_med=np.median(e), ms_per_eval=1e3 * np.median(t / e))


# ============================================================ w6
def w6(reg):
    ref = {(r['x'], r['y']): r for r in read('%s6_dd14' % reg, hp=True)}
    res = {'dd14': tstats(list(ref.values()))}
    per = {}
    for acc in ('dbl10', 'dbl12'):
        rs = [r for r in read('%s6_%s' % (reg, acc)) if (r['x'], r['y']) in ref]
        res[acc] = tstats(rs)
        G = np.array([digits(r['hi'][:6 * n], 0, ref[(r['x'], r['y'])]['hi'][:6 * n], ref[(r['x'], r['y'])]['lo'][:6 * n], cap=16)
                      for r in rs]).reshape(len(rs), 6, n)
        per[acc] = dict(G=G, pts=[(r['x'], r['y']) for r in rs], time=np.array([r['time'] for r in rs]))
    return res, per, ref


def report_w6(reg, res, per):
    L = ['## %s, weight 6, %d random points\n' % (REG[reg], res['dd14']['n'])]
    L.append('| run | median time [s] | mean | 10%% | 90%% | min | max | total [s] | evals (median) | ms / eval |')
    L.append('|---|---|---|---|---|---|---|---|---|---|')
    for acc, lab in (('dbl10', 'double, err 1e-10'), ('dbl12', 'double, err 1e-12'), ('dd14', 'dd_real, err 1e-14')):
        s = res[acc]
        L.append('| %s | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.0f | %d | %.2f |' % (lab, s['median'], s['mean'], s['p10'], s['p90'],
                 s['min'], s['max'], s['total'], s['evals_med'], s['ms_per_eval']))
    for acc in ('dbl10', 'dbl12'):
        G = per[acc]['G']; mn = np.nanmin(G, axis=2)                       # (points, weight)
        L.append('\n**double, err %s vs dd reference** -- worst component per point, by weight:\n' % ('1e-10' if acc == 'dbl10' else '1e-12'))
        L.append('| weight | worst point | 1% | 10% | median | points >= 12 | >= 10 | >= 8 |')
        L.append('|---|---|---|---|---|---|---|---|')
        for w in range(6):
            v = mn[:, w]
            L.append('| %d | %.1f | %.1f | %.1f | %.1f | %.0f%% | %.0f%% | %.0f%% |' % (w + 1, v.min(), np.percentile(v, 1), np.percentile(v, 10),
                     np.median(v), 100 * np.mean(v >= 12), 100 * np.mean(v >= 10), 100 * np.mean(v >= 8)))
        # which components are bad
        G6 = G[:, 5, :]
        worst_c = np.nanargmin(G6, axis=1)
        cnt = np.bincount(worst_c, minlength=n)
        med = np.nanmedian(G6, axis=0)
        bad = np.argsort(med)[:10]
        L.append('\nweight 6: components with the lowest median digits (index 0-based: median / worst / #points where it is the worst):')
        L.append(', '.join('%d: %.1f / %.1f / %d' % (c, med[c], np.nanmin(G6[:, c]), cnt[c]) for c in bad))
        frac = np.mean(G6 < 12, axis=0)
        L.append('components below 12 digits at >= 1 point: %d of %d; below 12 at > 50%% of points: %d' %
                 (np.sum(np.nanmin(G6, axis=0) < 12), np.sum(~np.all(np.isnan(G6), axis=0)), np.sum(frac > 0.5)))
    return '\n'.join(L) + '\n'


def coords(reg, pts):
    """plot coordinates: (x, y) Euclidean; (w, u) for scattering (x = -w/v, y = -u/v)."""
    a = np.array([(float(x), float(y)) for x, y in pts])
    if reg == 'e': return a
    X, Y = -a[:, 0], -a[:, 1]; v = 1 / (1 + X + Y)
    return np.c_[X * v, Y * v]


def plot_w6(allres):
    fig, ax = plt.subplots(2, 3, figsize=(7.4, 4.6), gridspec_kw=dict(wspace=0.35, hspace=0.45))
    for r_, reg in enumerate('es'):
        res, per, ref = allres[reg]
        # (a) timing distribution
        a = ax[r_, 0]
        tdd = np.array([v['time'] for v in ref.values()])
        bins = np.logspace(np.log10(0.3), np.log10(max(tdd.max(), 100)), 40)
        for acc, lab, col in (('dbl10', 'double 1e-10', '#9ec5f4'), ('dbl12', 'double 1e-12', '#256abf')):
            a.hist(per[acc]['time'], bins=bins, color=col, alpha=0.9, label=lab)
        a.hist(tdd, bins=bins, color='#d95f02', alpha=0.85, label='dd 1e-14')
        a.set_xscale('log'); a.set_xlabel('time per point [s]'); a.set_ylabel('points')
        a.set_title(REG[reg], color=P.INK, loc='left', fontsize=8.5)
        if r_ == 0: a.legend(frameon=False, fontsize=6.5)
        # (b) digits of double (1e-12) vs dd, worst component per point, per weight
        a = ax[r_, 1]
        G = per['dbl12']['G']; mn = np.nanmin(G, axis=2); md = np.nanmedian(G, axis=2)
        a.boxplot([mn[:, w] for w in range(6)], positions=np.arange(1, 7) - 0.17, widths=0.3, showfliers=True,
                  flierprops=dict(ms=1.5, mec='#184f95'), boxprops=dict(color='#184f95'), medianprops=dict(color='#184f95'),
                  whiskerprops=dict(color='#184f95'), capprops=dict(color='#184f95'))
        a.boxplot([md[:, w] for w in range(6)], positions=np.arange(1, 7) + 0.17, widths=0.3, showfliers=False,
                  boxprops=dict(color=P.INK2), medianprops=dict(color=P.INK2), whiskerprops=dict(color=P.INK2), capprops=dict(color=P.INK2))
        a.set_xticks(range(1, 7)); a.set_xticklabels(range(1, 7)); a.set_xlabel('weight')
        a.set_ylabel('digits (double 1e-12 vs dd)'); a.axhline(12, color=P.INK, lw=0.5, ls=(0, (3, 3)))
        a.set_ylim(4, 16.5)
        if r_ == 0: a.text(0.02, 0.03, 'blue: worst component\ngrey: median component', transform=a.transAxes, fontsize=6.5, color=P.INK2)
        # (c) per component, weight 6
        a = ax[r_, 2]
        G6 = G[:, 5, :]
        a.plot(np.arange(n), np.nanmedian(G6, axis=0), '.', ms=2, color=P.INK2, label='median over points')
        a.plot(np.arange(n), np.nanmin(G6, axis=0), '.', ms=2, color='#184f95', label='worst point')
        a.set_xlabel('component'); a.set_ylabel('digits at weight 6'); a.axhline(12, color=P.INK, lw=0.5, ls=(0, (3, 3)))
        a.set_ylim(4, 16.5)
        if r_ == 0: a.legend(frameon=False, fontsize=6.5, loc='lower left', markerscale=3)
        for a in ax[r_]:
            for s in ('top', 'right'): a.spines[s].set_visible(False)
    fig.savefig(os.path.join(H, 'w6_timing_digits.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(H, 'w6_timing_digits.png'), dpi=170, bbox_inches='tight')
    # map: worst digits over components at weight 6 (double 1e-12) and dd time, per point
    fig, ax = plt.subplots(2, 2, figsize=(6.0, 5.4), gridspec_kw=dict(wspace=0.3, hspace=0.35))
    for r_, reg in enumerate('es'):
        res, per, ref = allres[reg]
        c = coords(reg, per['dbl12']['pts']); mn = np.nanmin(per['dbl12']['G'][:, 5, :], axis=1)
        sc = ax[r_, 0].scatter(c[:, 0], c[:, 1], c=np.clip(mn, -2, 16), cmap=P.CMAP, norm=P.NORM, s=6, lw=0)
        ax[r_, 0].set_title('%s: worst digits, w6, double 1e-12' % REG[reg].split()[0], fontsize=8, color=P.INK)
        cd = coords(reg, list(ref.keys())); td = np.array([v['time'] for v in ref.values()])
        sc2 = ax[r_, 1].scatter(cd[:, 0], cd[:, 1], c=td, cmap='viridis', s=6, lw=0)
        ax[r_, 1].set_title('%s: dd time per point [s]' % REG[reg].split()[0], fontsize=8, color=P.INK)
        fig.colorbar(sc2, ax=ax[r_, 1], shrink=0.8)
        fig.colorbar(sc, ax=ax[r_, 0], shrink=0.8, ticks=[1, 4, 8, 12, 16])
        for a in ax[r_]:
            a.plot([0, 1, 0, 0], [0, 0, 1, 0], color=P.EDGE, lw=0.6); a.set_aspect('equal')
            a.set_xlabel('$x$' if reg == 'e' else '$w$  ($x=-w/v$)'); a.set_ylabel('$y$' if reg == 'e' else '$u$  ($y=-u/v$)')
    fig.savefig(os.path.join(H, 'w6_maps.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(H, 'w6_maps.png'), dpi=170, bbox_inches='tight')


# ============================================================ w2 lattice maps
def w2(reg, acc):
    refs = pickle.load(open(os.path.join(H, 'ref_%s2000.pkl' % ('eucl' if reg == 'e' else 'scat')), 'rb'))
    hp = acc.startswith('dd')
    rows, comp, pts, times = [], [], [], []
    for r in read('%s2_%s' % (reg, acc), hp=hp):
        R = refs[(r['x'], r['y'])]
        g = []
        for w in (0, 1):
            bh = np.array([hl(a)[0] + 1j * hl(b)[0] for a, b in R[w]]); bl = np.array([hl(a)[1] + 1j * hl(b)[1] for a, b in R[w]])
            g.append(digits(r['hi'][w * n:(w + 1) * n], r['lo'][w * n:(w + 1) * n], bh, bl, cap=16 if not hp else 32))
        pts.append((r['x'], r['y'])); times.append(r['time'])
        rows.append([np.nanmin(g[0]), np.nanmedian(g[0]), np.nanmin(g[1]), np.nanmedian(g[1])]); comp.append(g[1])
    return np.array(rows), np.array(comp), pts, np.array(times)


def plot_w2(data):
    fig, ax = plt.subplots(2, 4, figsize=(8.6, 4.6), gridspec_kw=dict(wspace=0.12, hspace=0.3))
    for r_, reg in enumerate('es'):
        for c_, acc in enumerate(('dbl12', 'dd16')):
            if (reg, acc) not in data: continue
            rows, comp, pts, t = data[(reg, acc)]
            xy = coords(reg, pts)
            d = np.c_[xy, 1 - xy.sum(1), np.ones(len(xy)), rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]]
            # snap lattice coordinates (scattering: (w, u) = (j, i)/65)
            d[:, 0] = np.round(d[:, 0] * 65) / 65; d[:, 1] = np.round(d[:, 1] * 65) / 65
            for k, col in enumerate((6, 7)):
                a = ax[r_, 2 * c_ + k]
                P.panel(a, d, col, 'none', show_x=(r_ == 1), show_y=(c_ == 0 and k == 0), soft=False)
                if reg == 's':
                    a.set_xlabel('$w$' if r_ == 1 else ''); a.set_ylabel('$u$' if (c_ == 0 and k == 0) else '')
                if r_ == 0:
                    a.set_title('%s, %s' % ('double 1e-12' if acc == 'dbl12' else 'dd 1e-16', 'worst comp.' if k == 0 else 'median comp.'),
                                fontsize=7.5, color=P.INK, pad=3)
        ax[r_, 0].annotate(REG[reg], xy=(0, 0.5), xycoords='axes fraction', xytext=(-30, 0), textcoords='offset points',
                           rotation=90, ha='center', va='center', color=P.INK, fontsize=8.5)
    P.colorbar(fig, ax)
    fig.savefig(os.path.join(H, 'w2_digits_maps.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(H, 'w2_digits_maps.png'), dpi=170, bbox_inches='tight')
    # per function (all 371 weight-2 functions)
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 4.2), sharex=True, gridspec_kw=dict(hspace=0.15))
    for r_, reg in enumerate('es'):
        a = ax[r_]
        for acc, col in (('dbl12', '#184f95'), ('dd16', '#d95f02')):
            if (reg, acc) not in data: continue
            comp = data[(reg, acc)][1]
            ok = ~np.all(np.isnan(comp), axis=0)
            idx = np.arange(n)[ok]
            a.plot(idx, np.nanmedian(comp[:, ok], axis=0), '.', ms=2.2, color=col, alpha=0.5)
            a.plot(idx, np.nanmin(comp[:, ok], axis=0), 'v', ms=1.8, color=col, label=('double 1e-12' if acc == 'dbl12' else 'dd 1e-16') + ': worst point (median: dots)')
        a.set_ylabel('digits, weight 2'); a.set_title(REG[reg], loc='left', fontsize=8, color=P.INK)
        a.axhline(12, color=P.INK, lw=0.5, ls=(0, (3, 3)))
        for s in ('top', 'right'): a.spines[s].set_visible(False)
    ax[0].legend(frameon=False, fontsize=6.5, markerscale=2.5, loc='lower left')
    ax[1].set_xlabel('weight-2 function (component index)')
    fig.savefig(os.path.join(H, 'w2_digits_per_function.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(H, 'w2_digits_per_function.png'), dpi=170, bbox_inches='tight')


def report_w2(data):
    L = ['## Weight 2, 2016-point lattice, vs exact GPLs (components with |J| >= 1e-4)\n',
         '| region | run | worst comp.: worst pt | 1% | median | median comp.: median | points >= 12 (worst comp.) | time median [s] | max |',
         '|---|---|---|---|---|---|---|---|---|']
    for (reg, acc), (rows, comp, pts, t) in sorted(data.items()):
        v = rows[:, 2]
        L.append('| %s | %s | %.1f | %.1f | %.1f | %.1f | %.0f%% | %.3f | %.2f |' % (REG[reg], acc, v.min(), np.percentile(v, 1), np.median(v),
                 np.median(rows[:, 3]), 100 * np.mean(v >= 12), np.median(t), t.max()))
    return '\n'.join(L) + '\n'


if __name__ == '__main__':
    what = sys.argv[1:] or ['w6', 'w2']
    rep = ''
    if 'w6' in what:
        allres = {}
        for reg in 'es':
            if glob.glob(os.path.join(H, '%s6_dd14_*.out' % reg)):
                res, per, ref = w6(reg); allres[reg] = (res, per, ref); rep += report_w6(reg, res, per)
        if len(allres) == 2: plot_w6(allres)
        pickle.dump({k: (v[0], {a: dict(G=p['G'], pts=p['pts'], time=p['time']) for a, p in v[1].items()}) for k, v in allres.items()},
                    open(os.path.join(H, 'w6_results.pkl'), 'wb'))
    if 'w2' in what:
        data = {}
        for reg in 'es':
            for acc in ('dbl12', 'dd16'):
                if glob.glob(os.path.join(H, '%s2_%s_*.out' % (reg, acc))) and os.path.exists(os.path.join(H, 'ref_%s2000.pkl' % ('eucl' if reg == 'e' else 'scat'))):
                    data[(reg, acc)] = w2(reg, acc)
        if data:
            plot_w2(data); rep += report_w2(data)
            pickle.dump({k: (v[0], v[2], v[3]) for k, v in data.items()}, open(os.path.join(H, 'w2_results.pkl'), 'wb'))
    open(os.path.join(H, 'REPORT.md'), 'a').write(rep)
    print(rep)
