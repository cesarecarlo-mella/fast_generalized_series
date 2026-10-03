#!/usr/bin/env python3
"""Collect the chunks and analyse.

  * double vs dd (reference), every weight 1..W, every component
  * double and dd vs the EXACT weight-1,2 solution (closed-form GPLs; scattering region: x+i0, y+i0,
    weight-2 GPLs by quadrature on [0, x]) -- also certifies the dd reference
  * evaluation time per point

Digits of a component = -log10 |a - b| / |b|, only where |b| >= 1e-4 (the 136 functions that vanish
identically at weights 1-2, and any other zero, are reported separately as absolute errors).

Output (results/): summary.md, digits.npz, ref_<region>.npz (exact cache),
  maps_<region>.pdf/png        double vs dd, rows worst / median component, one column per weight
  exact_<region>.pdf/png       weight 2 vs exact: double and dd
  per_function_<region>.pdf    per component: median and worst point (w2 vs exact; w=W double vs dd)
  timing.pdf/png               time per point histograms
Scattering-region maps are drawn in (u, w)  (x = -w/v, y = -u/v, v = 1 - u - w).

usage: 02_analyze.py [--jobs N]
"""
import argparse, glob, os, sys
from decimal import Decimal, getcontext
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'exact')); sys.path.insert(0, os.path.join(HERE, 'tools'))
getcontext().prec = 40
N_ = 371
REGION = {'eucl': 'Euclidean region', 'scat': 'scattering region'}


def hl(s):
    d = Decimal(s); h = float(d); return h, float(d - Decimal(h))


# ------------------------------------------------------------------ reading
def read_chunk(args):
    fn, hp = args
    pts = {}
    cur = None
    re_, im_ = [], []
    def flush():
        if cur is not None and not cur['failed']:
            if hp:
                a = np.array([hl(s) for s in re_]); b = np.array([hl(s) for s in im_])
                cur['hi'] = a[:, 0] + 1j * b[:, 0]; cur['lo'] = a[:, 1] + 1j * b[:, 1]
            else:
                cur['hi'] = np.array(re_, float) + 1j * np.array(im_, float); cur['lo'] = None
    for line in open(fn):
        if line.startswith('#'):
            flush(); re_, im_ = [], []
            p = line.split()
            cur = dict(failed='FAILED' in line)
            if not cur['failed']:
                cur.update(steps=int(p[3]), evals=int(p[4]), rej=int(p[5]), time=float(p[6]))
            pts[(p[1], p[2])] = cur
        else:
            a, b = line.split()
            if hp: re_.append(a); im_.append(b)
            else: re_.append(float(a)); im_.append(float(b))
    flush()
    return pts


def read_all(region, prec, pool):
    fs = sorted(glob.glob(os.path.join(HERE, 'out', '%s_%s_*.out' % (region, prec))))
    out = {}
    for d in pool.map(read_chunk, [(f, prec == 'dd') for f in fs]):
        out.update(d)
    return out


# ------------------------------------------------------------------ exact reference (weights 1, 2)
def exact_point(args):
    region, sa, sb = args
    import mpmath as mpm
    from exact_w2 import Exact, G, G_quad
    mpm.mp.dps = 30
    a, b = mpm.mpf(sa), mpm.mpf(sb)
    if region == 'eucl':
        x, y, gf = a, b, G
    else:
        u, v = a, b; x, y, gf = -(1 - u - v) / v, -u / v, G_quad
    res = []
    for w in (1, 2):
        vals = Exact(w, 30)(x, y, Gf=gf)
        for v_ in vals:
            v_ = mpm.mpc(v_)
            res.append((hl(mpm.nstr(v_.real, 34)), hl(mpm.nstr(v_.imag, 34))))
    r = np.array(res)
    return (sa, sb), r[:, 0, 0] + 1j * r[:, 1, 0], r[:, 0, 1] + 1j * r[:, 1, 1]


def exact_all(region, keys, pool):
    cache = os.path.join(HERE, 'results', 'ref_%s.npz' % region)
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        d = {tuple(k): (h, l) for k, h, l in zip(z['keys'], z['hi'], z['lo'])}
        if all(k in d for k in keys): return d
    d = {}
    for k, h, l in pool.imap_unordered(exact_point, [(region,) + k for k in keys], chunksize=4):
        d[k] = (h, l)
    ks = list(d)
    np.savez(cache, keys=np.array(ks), hi=np.array([d[k][0] for k in ks]), lo=np.array([d[k][1] for k in ks]))
    return d


# ------------------------------------------------------------------ digits
def digits(ahi, alo, bhi, blo, cap, chop=1e-4):
    d = ahi - bhi
    if alo is not None: d = d + alo
    if blo is not None: d = d - blo
    with np.errstate(divide='ignore', invalid='ignore'):
        g = -np.log10(np.abs(d) / np.abs(bhi))
    g = np.where(np.abs(d) == 0, cap, np.minimum(g, cap))
    g = np.where(np.isfinite(d), g, 0.0)                 # NaN / inf in the result counts as 0 digits
    return np.where(np.abs(bhi) >= chop, g, np.nan), np.where(np.abs(bhi) < chop, np.abs(d), np.nan)


def lattice_xy(region, keys):
    a = np.array([(float(p), float(q)) for p, q in keys])
    N = int(round(1 / a.min()))
    if region == 'scat':                       # (u, v) -> plot in (u, w)
        a = np.c_[a[:, 0], 1 - a[:, 0] - a[:, 1]]
    return np.round(a * N) / N, N


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    args = ap.parse_args()
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_w6_paper as P
    regions = [r for r in ('eucl', 'scat') if glob.glob(os.path.join(HERE, 'out', '%s_*.out' % r))]
    rep = ['# cluster_eval results\n']
    save = {}
    tim = {}
    with Pool(args.jobs) as pool:
        for R in regions:
            dbl, dd = read_all(R, 'double', pool), read_all(R, 'dd', pool)
            keys_all = [tuple(l.split()) for l in open(os.path.join(HERE, 'points', '%s.txt' % R)) if l.strip()]
            ok = [k for k in keys_all if k in dbl and k in dd and not dbl[k]['failed'] and not dd[k]['failed']]
            W = (len(dd[ok[0]]['hi']) - 1) // N_
            print(R, 'points', len(keys_all), 'complete (both precisions)', len(ok), 'weight', W, flush=True)
            ex = exact_all(R, ok, pool)
            # ---- timing
            rep.append('## %s: %d lattice points, weights 1..%d\n' % (REGION[R], len(keys_all), W))
            rep.append('| precision | done | failed | missing | time median [s] | mean | 10% | 90% | max | total [core-h] | evals median | ms / eval |')
            rep.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
            for lab, D in (('double', dbl), ('dd', dd)):
                t = np.array([v['time'] for v in D.values() if not v['failed']]); e = np.array([v['evals'] for v in D.values() if not v['failed']])
                nf = sum(v['failed'] for v in D.values())
                tim[(R, lab)] = t
                rep.append('| %s | %d | %d | %d | %.2f | %.2f | %.2f | %.2f | %.1f | %.2f | %d | %.2f |' % (lab, len(t), nf, len(keys_all) - len(D),
                           np.median(t), t.mean(), np.percentile(t, 10), np.percentile(t, 90), t.max(), t.sum() / 3600, np.median(e), 1e3 * np.median(t / e)))
            # ---- digits
            n = len(ok)
            Gd = np.full((n, W, N_), np.nan)            # double vs dd
            Ad = np.full((n, W), np.nan)                # max abs error on vanishing components
            Xd = np.full((n, 2, N_), np.nan); Xq = np.full((n, 2, N_), np.nan)   # double / dd vs exact
            for i, k in enumerate(ok):
                a, b = dbl[k], dd[k]
                for w in range(W):
                    s = slice(w * N_, (w + 1) * N_)
                    g, z = digits(a['hi'][s], None, b['hi'][s], b['lo'][s], 16)
                    Gd[i, w] = g; Ad[i, w] = np.nanmax(z) if np.any(~np.isnan(z)) else np.nan
                eh, el = ex[k]
                for w in range(min(W, 2)):
                    s = slice(w * N_, (w + 1) * N_)
                    Xd[i, w] = digits(a['hi'][s], None, eh[s], el[s], 16)[0]
                    Xq[i, w] = digits(b['hi'][s], b['lo'][s], eh[s], el[s], 32)[0]
            with np.errstate(all='ignore'):
                import warnings; warnings.simplefilter('ignore')
                mn, md = np.nanmin(Gd, 2), np.nanmedian(Gd, 2)
                xmn, xmd = np.nanmin(Xd, 2), np.nanmedian(Xd, 2)
                qmn = np.nanmin(Xq, 2)
            rep.append('\n**double vs dd reference** (per point: worst component; then the median component)\n')
            rep.append('| weight | worst point | 1% | 10% | median | median comp. (median) | points >= 12 | >= 10 | >= 8 | max abs err, vanishing comps |')
            rep.append('|---|---|---|---|---|---|---|---|---|---|')
            for w in range(W):
                v = mn[:, w]
                rep.append('| %d | %.1f | %.1f | %.1f | %.1f | %.1f | %.0f%% | %.0f%% | %.0f%% | %.1e |' % (w + 1, np.nanmin(v), np.nanpercentile(v, 1),
                           np.nanpercentile(v, 10), np.nanmedian(v), np.nanmedian(md[:, w]), 100 * np.mean(v >= 12), 100 * np.mean(v >= 10),
                           100 * np.mean(v >= 8), np.nanmax(Ad[:, w])))
            rep.append('\n**vs exact** (weights 1, 2; worst component per point)\n')
            rep.append('| run | weight | worst point | 1% | median | median comp. (median) | points >= 12 |')
            rep.append('|---|---|---|---|---|---|---|')
            for lab, M, MD in (('double', xmn, xmd), ('dd', qmn, None)):
                for w in range(min(W, 2)):
                    v = M[:, w]
                    rep.append('| %s | %d | %.1f | %.1f | %.1f | %s | %.0f%% |' % (lab, w + 1, np.nanmin(v), np.nanpercentile(v, 1), np.nanmedian(v),
                               '%.1f' % np.nanmedian(MD[:, w]) if MD is not None else '-', 100 * np.mean(v >= 12)))
            rep.append('')
            xy, NL = lattice_xy(R, ok)
            save.update({'%s_keys' % R: np.array(ok), '%s_xy' % R: xy, '%s_dbl_vs_dd' % R: Gd, '%s_dbl_vs_exact' % R: Xd,
                         '%s_dd_vs_exact' % R: Xq, '%s_abs_vanishing' % R: Ad,
                         '%s_time_double' % R: np.array([dbl[k]['time'] for k in ok]), '%s_time_dd' % R: np.array([dd[k]['time'] for k in ok])})
            # ---- maps: double vs dd, one column per weight
            lx, ly = ('$x$', '$y$') if R == 'eucl' else ('$u$', '$w$')
            fig, ax = plt.subplots(2, W, figsize=(1.75 * W + 0.8, 4.1), gridspec_kw=dict(wspace=0.1, hspace=0.12), squeeze=False)
            for w in range(W):
                d = np.c_[xy, 1 - xy.sum(1), np.ones(n), mn[:, w], md[:, w]]
                for r_, col in enumerate((4, 5)):
                    a_ = ax[r_, w]
                    P.panel(a_, d, col, 'none', show_x=(r_ == 1), show_y=(w == 0), soft=False)
                    a_.set_xlabel(lx if r_ == 1 else ''); a_.set_ylabel(ly if w == 0 else '')
                    if r_ == 0: a_.set_title('weight %d' % (w + 1), color=P.INK, pad=3)
            for r_, t in enumerate(('worst component', 'median component')):
                ax[r_, 0].annotate(t, xy=(0, 0.5), xycoords='axes fraction', xytext=(-32, 0), textcoords='offset points',
                                   rotation=90, ha='center', va='center', color=P.INK, fontsize=8.5)
            fig.suptitle('%s: double vs dd reference' % REGION[R], color=P.INK, fontsize=9, y=1.0)
            P.colorbar(fig, ax)
            for ext in ('pdf', 'png'): fig.savefig(os.path.join(HERE, 'results', 'maps_%s.%s' % (R, ext)), bbox_inches='tight', dpi=170)
            plt.close(fig)
            # ---- maps: weight 2 vs exact
            if W >= 2:
                fig, ax = plt.subplots(1, 4, figsize=(8.0, 2.4), gridspec_kw=dict(wspace=0.1), squeeze=False)
                for c_, (lab, A) in enumerate((('double', Xd), ('dd', Xq))):
                    with np.errstate(all='ignore'):
                        d = np.c_[xy, 1 - xy.sum(1), np.ones(n), np.nanmin(A[:, 1], 1), np.nanmedian(A[:, 1], 1)]
                    for k_, col in enumerate((4, 5)):
                        a_ = ax[0, 2 * c_ + k_]
                        P.panel(a_, d, col, 'none', show_x=True, show_y=(c_ == 0 and k_ == 0), soft=False)
                        a_.set_xlabel(lx); a_.set_ylabel(ly if (c_ == 0 and k_ == 0) else '')
                        a_.set_title('%s vs exact, %s' % (lab, 'worst comp.' if k_ == 0 else 'median comp.'), fontsize=7.5, color=P.INK, pad=3)
                fig.suptitle('%s, weight 2' % REGION[R], color=P.INK, fontsize=9, y=1.04)
                P.colorbar(fig, ax)
                for ext in ('pdf', 'png'): fig.savefig(os.path.join(HERE, 'results', 'exact_%s.%s' % (R, ext)), bbox_inches='tight', dpi=170)
                plt.close(fig)
            # ---- per function
            fig, ax = plt.subplots(2, 1, figsize=(7.4, 4.2), sharex=True, gridspec_kw=dict(hspace=0.2))
            for a_, (A, t) in zip(ax, ((Xd[:, 1] if W >= 2 else Xd[:, 0], 'weight %d, double vs exact' % min(W, 2)),
                                       (Gd[:, W - 1], 'weight %d, double vs dd' % W))):
                okc = ~np.all(np.isnan(A), axis=0); idx = np.arange(N_)[okc]
                with np.errstate(all='ignore'):
                    a_.plot(idx, np.nanmedian(A[:, okc], 0), '.', ms=2.4, color=P.INK2, label='median over points')
                    a_.plot(idx, np.nanmin(A[:, okc], 0), '.', ms=2.4, color='#184f95', label='worst point')
                a_.axhline(12, color=P.INK, lw=0.5, ls=(0, (3, 3))); a_.set_ylabel('digits'); a_.set_ylim(2, 16.5)
                a_.set_title('%s: %s' % (REGION[R], t), loc='left', fontsize=8, color=P.INK)
                for s_ in ('top', 'right'): a_.spines[s_].set_visible(False)
            ax[0].legend(frameon=False, fontsize=6.5, markerscale=2.5, loc='lower left'); ax[1].set_xlabel('component (0-based)')
            for ext in ('pdf', 'png'): fig.savefig(os.path.join(HERE, 'results', 'per_function_%s.%s' % (R, ext)), bbox_inches='tight', dpi=170)
            plt.close(fig)
    # ---- timing histograms
    if tim:
        fig, ax = plt.subplots(1, len(regions), figsize=(3.6 * len(regions), 2.6), squeeze=False)
        for a_, R in zip(ax[0], regions):
            allt = np.concatenate([tim[(R, p)] for p in ('double', 'dd')])
            bins = np.logspace(np.log10(max(allt.min(), 1e-3)) - 0.05, np.log10(allt.max()) + 0.05, 50)
            a_.hist(tim[(R, 'double')], bins=bins, color='#3987e5', label='double')
            a_.hist(tim[(R, 'dd')], bins=bins, color='#d95f02', label='dd_real')
            a_.set_xscale('log'); a_.set_xlabel('time per point [s] (one core)'); a_.set_ylabel('points')
            a_.set_title(REGION[R], loc='left', fontsize=8.5)
            for s_ in ('top', 'right'): a_.spines[s_].set_visible(False)
        ax[0, 0].legend(frameon=False, fontsize=7)
        for ext in ('pdf', 'png'): fig.savefig(os.path.join(HERE, 'results', 'timing.%s' % ext), bbox_inches='tight', dpi=170)
    np.savez_compressed(os.path.join(HERE, 'results', 'digits.npz'), **save)
    open(os.path.join(HERE, 'results', 'summary.md'), 'w').write('\n'.join(rep) + '\n')
    print('\n'.join(rep))


if __name__ == '__main__':
    main()
