#!/usr/bin/env python3
# =============================================================================
#  fast_selfconv.py -- FAST weight-6 (or any weight) SELF-CONVERGENCE check:
#  replaces  chop_phys_ber.wl + 00_build_coeffcache.wl + 01_eval_grid.wl +
#  02_merge_raw.wl + 03_analyze.wl  for the self-convergence case, with the
#  SAME output files (rawdata_*.m and compare_*.m, identical format), so
#  plot_compare.wl works unchanged.  No Mathematica needed.
# -----------------------------------------------------------------------------
#  Why the old path is slow: every worker re-reads the 250 MB phys_/ber_
#  files, and every (order pair) needs its own chopped files + Expand +
#  CoefficientRules (~30 min per worker, x50 workers, x every pair/mode).
#
#  What this does instead:
#   1. parse each CEILING file (phys_..._ord<C>.m, ber_..._ordp<C>_ordb<C>.m)
#      ONCE into a numeric monomial table (exponents + coefficient), cached as
#      .npz -- lower orders are NOT separate files, just exponent masks:
#      "order N" = both of the corner's variables <= N  (exactly what
#      chop_phys_ber.wl does, for phys_ AND ber_);
#   2. evaluate every monomial on the whole grid once, summed into buckets by
#      its "level" max(deg_A, deg_B)  ->  the value at EVERY order N is a
#      cumulative sum over levels, so all pairs (18/19, 19/20, 12/20, ...)
#      come out of one pass;
#   3. write rawdata_<mode>_w<w>_ordp<L>_ordb<H>.m and the analysed
#      compare_<mode>_w<w>_ordp<L>_ordb<H>.m (same digit formula, chop and
#      MIN/MEDIAN as 03_analyze.wl).
#
#  USAGE
#     python3 fast_selfconv.py --datadir <dir with ceiling phys_/ber_> --weight 6 \
#         --ceiling 20 --pairs 18,19 19,20 --modes c1only,c2only,c3only,coverage \
#         --grid results/grid_N1000.m --outdir results_w6_fast -n 6
#  (needs python3 + numpy + gmpy2 -- the blow_up_array venv works)
# =============================================================================
import argparse
import os
import re
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))

# corner -> (phys variables A,B) and (ber variables A,B)  [same as chop_phys_ber.wl]
PHYS_AB = {1: ('x', 'y'), 2: ('z', 'y'), 3: ('x', 'z')}
BER_AB = {1: ('sx', 'sy'), 2: ('sy', 'sz'), 3: ('sx', 'sz')}
# every variable a table may contain, and how to get its value at (x,y,z)
VARS = ['x', 'y', 'z', 'Log[x]', 'Log[y]', 'Log[z]', 'sx', 'sy', 'sz', 'LogX', 'LogY', 'LogZ']


def find_mparse():
    for d in (HERE, os.path.join(HERE, '..', 'blow_up_array', 'pipeline'),
              os.path.join(HERE, 'lib')):
        if os.path.isfile(os.path.join(d, 'mparse.py')):
            sys.path.insert(0, os.path.abspath(d)); return
    sys.exit("mparse.py not found (expected next to this script or in ../blow_up_array/pipeline)")


# ---------------------------------------------------------------- parsing
_SEP = re.compile(r'[(){}\[\],]')


def iter_components(s):
    """top-level list elements, one at a time."""
    s = re.sub(r'\(\*.*?\*\)', ' ', s, flags=re.S)
    s = re.sub(r'\\\r?\n\s*', '', s).strip()
    assert s[0] == '{' and s[-1] == '}'
    depth, start = 0, 1
    for m in _SEP.finditer(s, 1, len(s) - 1):
        c = m.group(0)
        if c in '({[':
            depth += 1
        elif c in ')}]':
            depth -= 1
        elif depth == 0:
            yield s[start:m.start()]; start = m.end()
    yield s[start:len(s) - 1]


def const_value(atom, e2):
    import math
    base = {'Pi': math.pi, 'E': math.e, 'EulerGamma': 0.5772156649015329, 'Log[2]': math.log(2)}
    if atom in base:
        v = base[atom]
    else:
        m = re.fullmatch(r'Zeta\[(\d+)\]', atom)
        if not m:
            return None
        import mpmath
        v = float(mpmath.zeta(int(m.group(1))))
    return v ** (e2 / 2)


def build_table(path, out_npz):
    """parse one .m file -> npz with comp(int32), exps(int16, M x len(VARS), in
    HALF-units), coef(float64: real part; the monomials are real on the grid)."""
    import numpy as np
    import mparse as mp
    t0 = time.time()
    s = open(path).read()
    comps, exps, coefs = [], [], []
    vidx = {v: i for i, v in enumerate(VARS)}
    for j, cs in enumerate(iter_components(s)):
        poly = mp.parse_string(cs)
        for k, q in poly.items():
            e = [0] * len(VARS)
            c = complex(q)
            for at, e2 in k:
                if at == 'I':
                    c *= 1j
                elif at in vidx:
                    e[vidx[at]] = e2
                elif at.startswith('Sqrt['):
                    e[vidx[at[5:-1]]] += 1
                else:
                    cv = const_value(at, e2)
                    if cv is None:
                        raise ValueError("unknown symbol %r in %s" % (at, path))
                    c *= cv
            comps.append(j); exps.append(e); coefs.append(c.real)
    ncomp = j + 1
    np.savez(out_npz, comp=np.array(comps, np.int32), exps=np.array(exps, np.int16).reshape(-1, len(VARS)),
             coef=np.array(coefs, np.float64), ncomp=ncomp)
    return "%s: %d components, %d monomials, %.0fs" % (os.path.basename(path), ncomp, len(coefs), time.time() - t0)


# ---------------------------------------------------------------- evaluation
def var_values(v, X, Y, Z):
    import numpy as np
    return {'x': X, 'y': Y, 'z': Z,
            'Log[x]': np.log(X), 'Log[y]': np.log(Y), 'Log[z]': np.log(Z),
            'sx': -np.log1p(-X), 'sy': -np.log1p(-Y), 'sz': -np.log1p(-Z),
            'LogX': np.log(X), 'LogY': np.log(Y), 'LogZ': np.log(Z)}[v]


def eval_levels(npz, pts, A, B, ceiling, chunk=4000, only_B=False):
    """-> array (ncomp, 2*ceiling+1, npts): sum of monomials per LEVEL
    (level = max(exp_A, exp_B) in half-units); value at order N =
    cumulative sum over levels 0..2N."""
    import numpy as np
    d = np.load(npz)
    comp, E, coef, ncomp = d['comp'], d['exps'].astype(np.int32), d['coef'], int(d['ncomp'])
    X, Y, Z = pts[:, 0], pts[:, 1], pts[:, 2]
    ia, ib = VARS.index(A), VARS.index(B)
    lev = E[:, ib] if only_B else np.maximum(E[:, ia], E[:, ib])
    keep = lev <= 2 * ceiling
    comp, E, coef, lev = comp[keep], E[keep], coef[keep], lev[keep]
    nlev = 2 * ceiling + 1
    out = np.zeros((ncomp * nlev, len(pts)))
    used = [i for i in range(len(VARS)) if E[:, i].any()]
    # power tables, indexed by the exponent in HALF-units (as stored):
    #   x, y, z, s*:  v^(k/2) = sqrt(v)^k     (v > 0 on the grid)
    #   Log*:         v^(k/2) with k even     (v < 0 is fine: integer powers)
    tabs = {}
    for i in used:
        vals = var_values(VARS[i], X, Y, Z)
        lo, hi = int(E[:, i].min()), int(E[:, i].max())
        ks = range(lo, hi + 1)
        if VARS[i].startswith('Log'):
            T = np.array([vals ** (k // 2) if k % 2 == 0 else np.full_like(vals, np.nan) for k in ks])
        else:
            r = np.sqrt(vals)
            T = np.array([r ** k for k in ks])
        tabs[i] = (lo, T)
    key = comp.astype(np.int64) * nlev + lev
    order = np.argsort(key, kind='stable')
    comp, E, coef, key = comp[order], E[order], coef[order], key[order]
    for s0 in range(0, len(coef), chunk):
        sl = slice(s0, s0 + chunk)
        m = np.repeat(coef[sl][:, None], len(pts), axis=1)
        for i in used:
            lo, T = tabs[i]
            m *= T[E[sl, i] - lo]
        k = key[sl]
        starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
        out[k[starts]] += np.add.reduceat(m, starts, axis=0)
    return out.reshape(ncomp, nlev, len(pts))


# ---------------------------------------------------------------- output (Mathematica format)
def mnum(v):
    if v != v or v in (float('inf'), float('-inf')):
        return "Indeterminate"
    s = repr(float(v))
    return s.replace('e', '*^') if 'e' in s else s


def mlist(vals):
    return "{" + ", ".join(mnum(v) for v in vals) + "}"


def digit_of(d, r):
    import math
    dd = abs(d)
    if dd < 1e-16:
        return 16.0
    rel = dd / abs(r)
    return min(16.0, -math.log10(rel + 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--datadir', required=True)
    ap.add_argument('--weight', type=int, default=6)
    ap.add_argument('--ceiling', type=int, required=True, help='order of the files in --datadir (e.g. 20)')
    ap.add_argument('--pairs', nargs='+', required=True, help='low,high pairs, e.g. 18,19 19,20')
    ap.add_argument('--modes', default='c1only,c2only,c3only,coverage')
    ap.add_argument('--grid', required=True, help='grid file from gen_grid.wl (list of {x,y,z})')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--cachedir', default=None, help='where parsed tables go (default <outdir>/cache)')
    ap.add_argument('--chop', type=float, default=1e-4, help='same as 03_analyze.wl chopThreshold')
    ap.add_argument('--phys-trunc', choices=['both', 'big'], default='both',
                    help="'both' (default): order N = SMALL and BIG <= N, as chop_phys_ber.wl does; "
                         "'big': only BIG <= N (the OLD 03_to_physical_vars.wl files -- blind to "
                         "non-convergence in SMALL; only for reproducing old results)")
    ap.add_argument('-n', '--ncores', type=int, default=os.cpu_count())
    a = ap.parse_args()
    find_mparse()
    import numpy as np
    import mparse as mp
    os.makedirs(a.outdir, exist_ok=True)
    cache = a.cachedir or os.path.join(a.outdir, 'cache')
    os.makedirs(cache, exist_ok=True)
    pairs = [tuple(int(v) for v in p.split(',')) for p in a.pairs]
    for L, H in pairs:
        if not (L < H <= a.ceiling):
            sys.exit("bad pair %d,%d (need low < high <= ceiling %d)" % (L, H, a.ceiling))
    modes = a.modes.split(',')
    corners = sorted({int(m[1]) for m in modes if m.endswith('only')} | ({1, 2, 3} if 'coverage' in modes else set()))
    w, C = a.weight, a.ceiling
    T0 = time.time()

    # 1. parse ceiling files once (cached), in parallel
    jobs = []
    for c in corners:
        for kind, fn in (('phys', 'phys_c%d_w%d_ord%d.m' % (c, w, C)),
                         ('ber', 'ber_c%d_w%d_ordp%d_ordb%d.m' % (c, w, C, C))):
            src = os.path.join(a.datadir, fn)
            if not os.path.isfile(src):
                sys.exit("missing " + src)
            dst = os.path.join(cache, fn[:-2] + '.npz')
            if not os.path.isfile(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
                jobs.append((src, dst))
    if jobs:
        print("[1] parsing %d ceiling file(s) once -> %s" % (len(jobs), cache), flush=True)
        import multiprocessing as mproc
        with mproc.get_context('fork').Pool(min(a.ncores, len(jobs))) as pool:
            for msg in pool.starmap(build_table, jobs):
                print("    " + msg, flush=True)
    else:
        print("[1] parsed tables already cached in %s" % cache)

    # 2. grid + evaluation (all orders at once)
    grid = mp.parse_file(a.grid)
    pts = np.array([[float(mp.pconst(v)) for v in p] for p in grid])
    print("[2] grid %s: %d points" % (a.grid, len(pts)), flush=True)
    t1 = time.time()
    vals = {}
    for c in corners:
        for kind, (A, B), fn in (('phys', PHYS_AB[c], 'phys_c%d_w%d_ord%d' % (c, w, C)),
                                 ('ber', BER_AB[c], 'ber_c%d_w%d_ordp%d_ordb%d' % (c, w, C, C))):
            lv = eval_levels(os.path.join(cache, fn + '.npz'), pts, A, B, C,
                             only_B=(kind == 'phys' and a.phys_trunc == 'big'))
            vals[(c, kind)] = (np.cumsum(lv, axis=1), lv)    # cs[:, 2N, :] = value at order N
            print("    corner %d %-4s evaluated (%d comps x %d orders x %d pts)  %.0fs"
                  % (c, kind, lv.shape[0], C + 1, len(pts), time.time() - t1), flush=True)

    # 3. write rawdata + compare files
    def pick(xv, yv, zv):          # same as 01_eval_grid.wl pickC
        i = int(np.argmax([xv, yv, zv]))
        return {0: 2, 1: 3, 2: 1}[i]
    for mode in modes:
        for L, H in pairs:
            tag = "%s_w%d_ordp%d_ordb%d" % (mode, w, L, H)
            rows, crow = [], []
            for j, (xv, yv, zv) in enumerate(pts):
                c = int(mode[1]) if mode.endswith('only') else pick(xv, yv, zv)
                (cN, lN), (cB, lB) = vals[(c, 'phys')], vals[(c, 'ber')]
                ref = cN[:, 2 * H, j]
                # low - high = -(terms of levels L+1/2 .. H), summed DIRECTLY: no cancellation
                # between two large totals (the old path's diffs had a ~1e-11 round-off floor)
                dN = -lN[:, 2 * L + 1:2 * H + 1, j].sum(axis=1)
                dB = -lB[:, 2 * L + 1:2 * H + 1, j].sum(axis=1)
                rows.append("{%s, %s, %s, %d, %s, %s, %s}" % (mnum(xv), mnum(yv), mnum(zv), c,
                                                             mlist(dN), mlist(dB), mlist(ref)))
                keep = np.abs(ref) >= a.chop
                gN = [digit_of(d, r) for d, r in zip(dN[keep], ref[keep])]
                gB = [digit_of(d, r) for d, r in zip(dB[keep], ref[keep])]
                if gN:
                    crow.append("{%s, %s, %s, %d, %s, %s, %s, %s}" % (
                        mnum(xv), mnum(yv), mnum(zv), c, mnum(min(gN)), mnum(float(np.median(gN))),
                        mnum(min(gB)), mnum(float(np.median(gB)))))
            with open(os.path.join(a.outdir, "rawdata_%s.m" % tag), 'w') as f:
                f.write('<|"mode" -> "%s", "weight" -> %d, "ordp" -> %d, "ordb" -> %d, '
                        '"header" -> {"x", "y", "z", "corner", "diffN", "diffB", "ref"}, '
                        '"data" -> {%s}|>\n' % (mode, w, L, H, ",\n ".join(rows)))
            with open(os.path.join(a.outdir, "compare_%s.m" % tag), 'w') as f:
                f.write('<|"mode" -> "%s", "weight" -> %d, "ordp" -> %d, "ordb" -> %d, "chopThreshold" -> %s, '
                        '"header" -> {"x", "y", "z", "corner", "digits_native_min", "digits_native_median", '
                        '"digits_bernoulli_min", "digits_bernoulli_median"}, "data" -> {%s}|>\n'
                        % (mode, w, L, H, mnum(a.chop), ",\n ".join(crow)))
            print("[3] wrote rawdata_%s.m + compare_%s.m" % (tag, tag), flush=True)
    print("DONE in %.1f min" % ((time.time() - T0) / 60))


if __name__ == '__main__':
    main()
