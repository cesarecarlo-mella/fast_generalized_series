#!/usr/bin/env python3
# =============================================================================
#  fast_exact.py -- series vs EXACT solution, all truncation orders in one pass.
#
#  The weight-3 version of the comparison behind the paper figures
#  (cluster_results/results_w2_consistent, weight 2).  Same definitions:
#
#    * series at order N  = the ceiling-order (20) phys_/ber_ files truncated
#      to BOTH of the corner's variables <= N  (exactly chop_phys_ber.wl;
#      for phys_ and ber_ alike, as in fast_selfconv.py);
#    * reference          = the exact solution, eps^w coefficient of
#      HW3[i] in H_exact_sol/H_solution_w3.m (weights 0..3), evaluated with
#      exact_eval.py (mpmath, 30 digits, real part, kept as double-double);
#    * diff               = series_N - exact   (signed), formed in extended
#      precision (hp_eval.py: exact coefficients, longdouble sums) so that
#      16 agreeing digits come out as 16 -- like the 20-digit Mathematica
#      evaluation of 01_eval_grid.wl.  Needs x86-64 (longdouble = 80 bit);
#    * digits             = min(16, -log10 |diff|/|ref|), components with
#      |ref| < chop (1e-4) dropped, MIN and MEDIAN over components
#      (same formula as 03_analyze.wl / fast_selfconv.py);
#    * modes c1only/c2only/c3only (one corner everywhere) and coverage
#      (nearest corner, pickC of 01_eval_grid.wl).
#
#  Output files have the SAME names and format as the weight-2 ones:
#     <outdir>/rawdata_<mode>_w<w>_ordp<N>_ordb<N>.m   {x,y,z,corner,diffN,diffB,ref}
#     <outdir>/compare_<mode>_w<w>_ordp<N>_ordb<N>.m   {x,y,z,corner,minN,medN,minB,medB}
#  so plot_w6_paper.py / plot_grid.py / plot_compare.wl work unchanged
#  (pass --weight 3).
#
#  USAGE
#     python3 fast_exact.py --datadir <dir with phys_c*_w3_ord20.m, ber_c*_w3_ordp20_ordb20.m> \
#         --weight 3 --ceiling 20 --orders 12 16 20 --modes c1only,c2only,c3only,coverage \
#         --grid grids/grid_N2000.m --exact exact_solutions/H_solution_w3.m \
#         --outdir ../cluster_results/results_w3_consistent -n 8
#  (python3 + numpy + gmpy2 + mpmath; the blow_up_array venv works)
# =============================================================================
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--datadir', required=True)
    ap.add_argument('--weight', type=int, default=3)
    ap.add_argument('--ceiling', type=int, default=20)
    ap.add_argument('--orders', nargs='+', type=int, default=[12, 16, 20])
    ap.add_argument('--modes', default='c1only,c2only,c3only,coverage')
    ap.add_argument('--grid', required=True)
    ap.add_argument('--exact', required=True, help='H_solution_w3.m (HW3[i] = {eps^0..eps^3})')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--cachedir', default=None)
    ap.add_argument('--chop', type=float, default=1e-4)
    ap.add_argument('-n', '--ncores', type=int, default=os.cpu_count())
    a = ap.parse_args()
    import numpy as np
    import fast_selfconv as fs
    fs.find_mparse()
    import mparse as mp
    import exact_eval as ee
    import hp_eval as hp
    hp.check_longdouble()

    w, C = a.weight, a.ceiling
    if not 0 <= w <= 3:
        sys.exit("weight must be 0..3 (H_solution_w3.m)")
    for N in a.orders:
        if not 1 <= N <= C:
            sys.exit("bad order %d (ceiling %d)" % (N, C))
    os.makedirs(a.outdir, exist_ok=True)
    cache = a.cachedir or os.path.join(a.outdir, 'cache')
    os.makedirs(cache, exist_ok=True)
    modes = a.modes.split(',')
    corners = sorted({int(m[1]) for m in modes if m.endswith('only')} | ({1, 2, 3} if 'coverage' in modes else set()))
    T0 = time.time()

    # 1. parse ceiling series files once (cached)
    jobs = []
    for c in corners:
        for fn in ('phys_c%d_w%d_ord%d.m' % (c, w, C), 'ber_c%d_w%d_ordp%d_ordb%d.m' % (c, w, C, C)):
            src = os.path.join(a.datadir, fn)
            if not os.path.isfile(src):
                sys.exit("missing " + src)
            dst = os.path.join(cache, fn[:-2] + '.hp.npz')
            if not os.path.isfile(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
                jobs.append((src, dst))
    if jobs:
        print("[1] parsing %d series file(s) -> %s" % (len(jobs), cache), flush=True)
        import multiprocessing as mproc
        with mproc.get_context('fork').Pool(min(a.ncores, len(jobs))) as pool:
            for msg in pool.starmap(hp.build_table_hp, jobs):
                print("    " + msg, flush=True)
    else:
        print("[1] series tables already cached in %s" % cache)

    # 2. grid, exact reference (cached), series at all orders
    grid = mp.parse_file(a.grid)
    gq = [[mp.pconst(v) for v in p] for p in grid]                 # exact rationals x, y, z
    for p in gq:
        assert p[0] + p[1] + p[2] == 1, "grid point with x+y+z != 1"
    # x, y, z to extended precision from the EXACT rationals: the series and the exact solution
    # must see the same point (some components move by ~1e-14 for a 1e-17 shift of x)
    import mpmath as mpm
    mpm.mp.dps = 30
    ld = lambda q: np.longdouble(mpm.nstr(mpm.mpf(int(q.numerator)) / int(q.denominator), 25))
    pts_ld = np.array([[ld(v) for v in p] for p in gq], dtype=np.longdouble)
    pts = pts_ld.astype(float)
    print("[2] grid %s: %d points" % (a.grid, len(pts)), flush=True)
    exf = os.path.join(cache, 'exact_w%d_%s.npz' % (w, os.path.splitext(os.path.basename(a.grid))[0]))
    if os.path.isfile(exf) and os.path.getmtime(exf) > os.path.getmtime(a.exact):
        d = np.load(exf)
        hi, lo = d['hi'], d['lo']
        print("    exact reference: cached %s" % exf)
    else:
        t1 = time.time()
        hi, lo = ee.evaluate(a.exact, w, gq, a.ncores)
        np.savez(exf, hi=hi, lo=lo)
        print("    exact reference eps^%d evaluated: %d comps x %d pts  %.0fs" % (w, hi.shape[0], hi.shape[1], time.time() - t1), flush=True)
    ex = hi.astype(np.longdouble) + lo
    t1 = time.time()
    vals = {}
    for c in corners:
        for kind, (A, B), fn in (('phys', fs.PHYS_AB[c], 'phys_c%d_w%d_ord%d' % (c, w, C)),
                                 ('ber', fs.BER_AB[c], 'ber_c%d_w%d_ordp%d_ordb%d' % (c, w, C, C))):
            vals[(c, kind)] = hp.eval_orders(os.path.join(cache, fn + '.hp.npz'), pts_ld, A, B, a.orders)
            print("    corner %d %-4s evaluated (extended precision)  %.0fs" % (c, kind, time.time() - t1), flush=True)
    assert ex.shape[0] == vals[(corners[0], 'phys')][a.orders[0]].shape[0], "component count: exact vs series"

    # 3. rawdata + compare files
    def pick(j):                   # 01_eval_grid.wl pickC = Ordering[{x,y,z}, -1]: on a tie the LAST maximum
        v = gq[j]
        i = max(range(3), key=lambda k: (v[k], k))
        return {0: 2, 1: 3, 2: 1}[i]
    note = "order N = order-%d array result truncated to both variables <= N; reference = exact HW3 eps^%d" % (C, w)
    for mode in modes:
        for N in a.orders:
            tag = "%s_w%d_ordp%d_ordb%d" % (mode, w, N, N)
            rows, crow = [], []
            for j, (xv, yv, zv) in enumerate(pts):
                c = int(mode[1]) if mode.endswith('only') else pick(j)
                refL = ex[:, j]
                dN = (vals[(c, 'phys')][N][:, j] - refL).astype(float)     # difference in extended precision
                dB = (vals[(c, 'ber')][N][:, j] - refL).astype(float)
                ref = refL.astype(float)
                rows.append("{%s, %s, %s, %d, %s, %s, %s}" % (fs.mnum(xv), fs.mnum(yv), fs.mnum(zv), c,
                                                             fs.mlist(dN), fs.mlist(dB), fs.mlist(ref)))
                keep = np.abs(ref) >= a.chop
                gN = [fs.digit_of(d, r) for d, r in zip(dN[keep], ref[keep])]
                gB = [fs.digit_of(d, r) for d, r in zip(dB[keep], ref[keep])]
                if gN:
                    crow.append("{%s, %s, %s, %d, %s, %s, %s, %s}" % (
                        fs.mnum(xv), fs.mnum(yv), fs.mnum(zv), c, fs.mnum(min(gN)), fs.mnum(float(np.median(gN))),
                        fs.mnum(min(gB)), fs.mnum(float(np.median(gB)))))
            with open(os.path.join(a.outdir, "rawdata_%s.m" % tag), 'w') as f:
                f.write('<|"mode" -> "%s", "weight" -> %d, "ordp" -> %d, "ordb" -> %d, "note" -> "%s", '
                        '"header" -> {"x", "y", "z", "corner", "diffN", "diffB", "ref"}, '
                        '"data" -> {%s}|>\n' % (mode, w, N, N, note, ",\n ".join(rows)))
            with open(os.path.join(a.outdir, "compare_%s.m" % tag), 'w') as f:
                f.write('<|"mode" -> "%s", "weight" -> %d, "ordp" -> %d, "ordb" -> %d, "chopThreshold" -> %s, "note" -> "%s", '
                        '"header" -> {"x", "y", "z", "corner", "digits_native_min", "digits_native_median", '
                        '"digits_bernoulli_min", "digits_bernoulli_median"}, "data" -> {%s}|>\n'
                        % (mode, w, N, N, fs.mnum(a.chop), note, ",\n ".join(crow)))
            print("[3] wrote rawdata_%s.m + compare_%s.m" % (tag, tag), flush=True)
    print("DONE in %.1f min" % ((time.time() - T0) / 60))


if __name__ == '__main__':
    main()
