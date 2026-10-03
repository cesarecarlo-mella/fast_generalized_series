"""weight-2 lattice maps (Euclidean: (x, y); scattering: (u, v) chart, integrated in u, v) vs exact GPLs,
and weight-6 double vs dd in the scattering region.   usage: python3 analyze_uv.py [w2] [w6]"""
import sys, os, pickle
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'camp'))
import analyze as AN          # read helpers, digits, plot style
import matplotlib.pyplot as plt
P = AN.P
H = os.path.dirname(os.path.abspath(__file__))
n = 371


def readfile(fn, hp=False):
    AN.H, keep = H, AN.H
    tag = os.path.basename(fn)[:-4]
    # AN.read expects <tag>_A.out / _B.out: emulate with a single file
    out = []
    cur = None
    for line in open(fn):
        if line.startswith('#'):
            p = line.split()
            cur = dict(x=p[1], y=p[2], failed='FAILED' in line, v=[])
            if not cur['failed']:
                cur.update(steps=int(p[3]), evals=int(p[4]), rej=int(p[5]), time=float(p[6]))
            out.append(cur)
        else:
            a, b = line.split()
            cur['v'].append((AN.hl(a), AN.hl(b)) if hp else (float(a), float(b)))
    AN.H = keep
    res = []
    for r in out:
        v = r.pop('v')
        if r['failed'] or not v: res.append(r); continue
        r['hi'] = np.array([a[0] + 1j * b[0] for a, b in v]) if hp else np.array([a + 1j * b for a, b in v])
        r['lo'] = np.array([a[1] + 1j * b[1] for a, b in v]) if hp else np.zeros(len(v), complex)
        res.append(r)
    return res


def lattice(fn, refs, cap=16):
    rows, comp, pts, times = [], [], [], []
    for r in readfile(fn):
        if r['failed']: continue
        R = refs[(r['x'], r['y'])]
        g = []
        for w in (0, 1):
            bh = np.array([AN.hl(a)[0] + 1j * AN.hl(b)[0] for a, b in R[w]]); bl = np.array([AN.hl(a)[1] + 1j * AN.hl(b)[1] for a, b in R[w]])
            g.append(AN.digits(r['hi'][w * n:(w + 1) * n], r['lo'][w * n:(w + 1) * n], bh, bl, cap=cap))
        rows.append([np.nanmin(g[0]), np.nanmedian(g[0]), np.nanmin(g[1]), np.nanmedian(g[1])]); comp.append(g[1])
        pts.append((float(r['x']), float(r['y']))); times.append(r['time'])
    return np.array(rows), np.array(comp), np.array(pts), np.array(times)


def w2():
    E = lattice(os.path.join(H, 'e2_dbl12.out'), pickle.load(open(os.path.join(H, '..', 'camp', 'ref_eucl2000.pkl'), 'rb')))
    S = lattice(os.path.join(H, 's2_dbl12.out'), pickle.load(open(os.path.join(H, 'ref_uv_lat.pkl'), 'rb')))
    rep = ['## Weight 2, 2016-point lattice, double, requested error 1e-12, vs exact GPLs (all 371 functions, |J| >= 1e-4)\n',
           '| region (integration variables) | worst comp.: worst point | 1% | 10% | median | median comp.: median | points >= 12 | >= 10 | time median [s] | 90% | max |',
           '|---|---|---|---|---|---|---|---|---|---|---|']
    for lab, (rows, comp, pts, t) in (('Euclidean (x, y)', E), ('scattering (u, v)', S)):
        v = rows[:, 2]
        rep.append('| %s | %.1f | %.1f | %.1f | %.1f | %.1f | %.0f%% | %.0f%% | %.3f | %.3f | %.2f |' % (lab, v.min(), np.percentile(v, 1),
                   np.percentile(v, 10), np.median(v), np.median(rows[:, 3]), 100 * np.mean(v >= 12), 100 * np.mean(v >= 10),
                   np.median(t), np.percentile(t, 90), t.max()))
    # maps
    fig, ax = plt.subplots(1, 4, figsize=(8.6, 2.5), gridspec_kw=dict(wspace=0.12))
    for c_, (lab, (rows, comp, pts, t)) in enumerate((('Euclidean', E), ('scattering', S))):
        xy = np.round(pts * 65) / 65
        d = np.c_[xy, 1 - xy.sum(1), np.ones(len(xy)), rows]
        for k, col in enumerate((6, 7)):
            a = ax[2 * c_ + k]
            P.panel(a, d, col, 'none', show_x=True, show_y=(k == 0), soft=False)
            a.set_title('%s, %s' % (lab, 'worst comp.' if k == 0 else 'median comp.'), fontsize=7.5, color=P.INK, pad=3)
            if c_ == 1:
                a.set_xlabel('$u$', labelpad=1); a.set_ylabel('$v$' if k == 0 else '', labelpad=1)
    P.colorbar(fig, ax[None, :])
    fig.savefig(os.path.join(H, 'w2_digits_maps.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(H, 'w2_digits_maps.png'), dpi=170, bbox_inches='tight')
    # per function
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 4.0), sharex=True, gridspec_kw=dict(hspace=0.18))
    for a, (lab, (rows, comp, pts, t)) in zip(ax, (('Euclidean region, (x, y)', E), ('scattering region, (u, v)', S))):
        ok = ~np.all(np.isnan(comp), axis=0); idx = np.arange(n)[ok]
        a.plot(idx, np.nanmedian(comp[:, ok], axis=0), '.', ms=2.4, color=P.INK2, label='median over the 2016 points')
        a.plot(idx, np.nanmin(comp[:, ok], axis=0), '.', ms=2.4, color='#184f95', label='worst point')
        a.axhline(12, color=P.INK, lw=0.5, ls=(0, (3, 3))); a.set_ylabel('digits, weight 2'); a.set_ylim(4, 16.5)
        a.set_title(lab, loc='left', fontsize=8, color=P.INK)
        for s in ('top', 'right'): a.spines[s].set_visible(False)
    ax[0].legend(frameon=False, fontsize=6.5, markerscale=2.5, loc='lower left'); ax[1].set_xlabel('weight-2 function (component index, 0-based)')
    fig.savefig(os.path.join(H, 'w2_digits_per_function.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(H, 'w2_digits_per_function.png'), dpi=170, bbox_inches='tight')
    pickle.dump(dict(E=E, S=S), open(os.path.join(H, 'w2_results.pkl'), 'wb'))
    return '\n'.join(rep) + '\n'


def w6():
    ref = {}
    for h in 'AB':
        for r in readfile(os.path.join(H, 's6_dd14_%s.out' % h), hp=True): ref[(r['x'], r['y'])] = r
    rep = ['## Scattering region, (u, v) chart, weight 6, %d random points: double vs dd (1e-14) reference\n' % len(ref)]
    rep.append('| run | time median [s] | mean | max | evals median | worst digits w1..w6 (worst point) | median over points of worst comp., w1..w6 | points with w6 >= 10 |')
    rep.append('|---|---|---|---|---|---|---|---|')
    td = np.array([r['time'] for r in ref.values() if not r['failed']])
    ed = np.array([r['evals'] for r in ref.values() if not r['failed']])
    rep.append('| dd 1e-14 | %.2f | %.2f | %.2f | %d | (reference) | | |' % (np.median(td), td.mean(), td.max(), np.median(ed)))
    for acc, lab in (('dbl10', 'double abs 1e-10'), ('rel12', 'double abs+rel 1e-12'), ('rel13', 'double abs+rel 1e-13')):
        rs = [r for h in 'AB' for r in readfile(os.path.join(H, 's6_%s_%s.out' % (acc, h))) if not r['failed']]
        G = np.array([[np.nanmin(AN.digits(r['hi'][w * n:(w + 1) * n], 0, ref[(r['x'], r['y'])]['hi'][w * n:(w + 1) * n],
                                           ref[(r['x'], r['y'])]['lo'][w * n:(w + 1) * n], cap=16)) for w in range(6)] for r in rs])
        t = np.array([r['time'] for r in rs]); e = np.array([r['evals'] for r in rs])
        rep.append('| %s | %.2f | %.2f | %.2f | %d | %s | %s | %.0f%% |' % (lab, np.median(t), t.mean(), t.max(), np.median(e),
                   ' '.join('%.1f' % v for v in G.min(0)), ' '.join('%.1f' % v for v in np.median(G, 0)), 100 * np.mean(G[:, 5] >= 10)))
    return '\n'.join(rep) + '\n'


if __name__ == '__main__':
    what = sys.argv[1:] or ['w2', 'w6']
    rep = ''.join(w2() if k == 'w2' else w6() for k in what)
    open(os.path.join(H, 'REPORT.md'), 'a').write(rep); print(rep)
