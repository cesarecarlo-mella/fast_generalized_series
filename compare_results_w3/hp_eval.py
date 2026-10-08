#!/usr/bin/env python3
# =============================================================================
#  hp_eval.py -- extended-precision evaluation of the phys_/ber_ series on a
#  grid, for the comparison against the EXACT solution.
#
#  Why: fast_selfconv.py evaluates the series in float64.  That is fine for
#  the self-convergence check (it sums the dropped levels directly), but for
#  series - exact the full series value is needed, and the monomial sum
#  cancels: in float64 the result is only good to ~1e-12..1e-14 relative,
#  so points that agree to 16 digits would show 12-14.  (The weight-2 paper
#  data got 16 there because 01_eval_grid.wl evaluates in Mathematica with
#  20-digit arithmetic.)
#
#  What this does:
#    * coefficients are kept as a double-double pair (hi, lo) built from the
#      EXACT rational times Pi^a Zeta[3]^b (mpmath, 40 digits);
#    * the monomials are evaluated and summed in numpy longdouble (x86-64:
#      80-bit, 64-bit mantissa, ~19 digits), and the sum is compensated
#      (Kahan) over chunks;
#    * only the requested orders N are accumulated (level <= 2N, level =
#      max(deg_A, deg_B) in half-units), as in chop_phys_ber.wl.
#  On machines where longdouble == float64 (arm64 / Apple Silicon) it stops
#  with an error instead of silently losing digits.
# =============================================================================
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

VARS = ['x', 'y', 'z', 'Log[x]', 'Log[y]', 'Log[z]', 'sx', 'sy', 'sz', 'LogX', 'LogY', 'LogZ']


def check_longdouble():
    import numpy as np
    if np.finfo(np.longdouble).eps > 1e-18:
        sys.exit("numpy longdouble has no extended precision on this machine (eps=%g): run on x86-64 Linux"
                 % np.finfo(np.longdouble).eps)


def _const(atom):
    import mpmath as mpm
    if atom == 'Pi':
        return mpm.pi
    m = re.fullmatch(r'Zeta\[(\d+)\]', atom)
    if m:
        return mpm.zeta(int(m.group(1)))
    if atom == 'Log[2]':
        return mpm.log(2)
    if atom == 'EulerGamma':
        return mpm.euler
    return None


def build_table_hp(path, out_npz):
    """.m series file -> npz: comp, exps (half-units), chi, clo (double-double coefficient)."""
    import numpy as np
    import mpmath as mpm
    import fast_selfconv as fs
    fs.find_mparse()
    import mparse as mp
    mpm.mp.dps = 40
    t0 = time.time()
    s = open(path).read()
    vidx = {v: i for i, v in enumerate(VARS)}
    comps, exps, chi, clo = [], [], [], []
    cache = {}
    for j, cs in enumerate(fs.iter_components(s)):
        poly = mp.parse_string(cs)
        for k, q in poly.items():
            e = [0] * len(VARS)
            cst = []
            for at, e2 in k:
                if at in vidx:
                    e[vidx[at]] = e2
                elif at.startswith('Sqrt['):
                    e[vidx[at[5:-1]]] += 1
                elif at == 'I':
                    raise ValueError("imaginary unit in %s (not supported)" % path)
                else:
                    cst.append((at, e2))
            key = tuple(sorted(cst))
            cv = cache.get(key)
            if cv is None:
                cv = mpm.mpf(1)
                for at, e2 in key:
                    b = _const(at)
                    if b is None:
                        raise ValueError("unknown symbol %r in %s" % (at, path))
                    cv *= b ** (mpm.mpf(e2) / 2)
                cache[key] = cv
            v = mpm.mpf(q.numerator) / q.denominator * cv
            h = float(v)
            comps.append(j); exps.append(e); chi.append(h); clo.append(float(v - h))
    np.savez(out_npz, comp=np.array(comps, np.int32), exps=np.array(exps, np.int16).reshape(-1, len(VARS)),
             chi=np.array(chi), clo=np.array(clo), ncomp=j + 1)
    return "%s: %d components, %d monomials, %.0fs" % (os.path.basename(path), j + 1, len(chi), time.time() - t0)


def _vals(v, X, Y, Z):
    import numpy as np
    one = np.longdouble(1)
    return {'x': X, 'y': Y, 'z': Z, 'Log[x]': np.log(X), 'Log[y]': np.log(Y), 'Log[z]': np.log(Z),
            'sx': -np.log1p(-X), 'sy': -np.log1p(-Y), 'sz': -np.log1p(-Z),
            'LogX': np.log(X), 'LogY': np.log(Y), 'LogZ': np.log(Z)}[v]


def eval_orders(npz, pts, A, B, orders, chunk=2000):
    """-> dict N -> longdouble array (ncomp, npts): series truncated to A, B <= N."""
    import numpy as np
    check_longdouble()
    LD = np.longdouble
    d = np.load(npz)
    comp, E, ncomp = d['comp'], d['exps'].astype(np.int32), int(d['ncomp'])
    coef = d['chi'].astype(LD) + d['clo'].astype(LD)
    X, Y, Z = (pts[:, 0].astype(LD), pts[:, 1].astype(LD), pts[:, 2].astype(LD))
    ia, ib = VARS.index(A), VARS.index(B)
    lev = np.maximum(E[:, ia], E[:, ib])
    Nmax = max(orders)
    keep = lev <= 2 * Nmax
    comp, E, coef, lev = comp[keep], E[keep], coef[keep], lev[keep]
    used = [i for i in range(len(VARS)) if E[:, i].any()]
    tabs = {}
    for i in used:
        vals = _vals(VARS[i], X, Y, Z)
        lo, hi = int(E[:, i].min()), int(E[:, i].max())
        if VARS[i].startswith('Log'):
            T = np.array([vals ** (k // 2) if k % 2 == 0 else np.full_like(vals, np.nan) for k in range(lo, hi + 1)])
        else:
            r = np.sqrt(vals)
            T = np.array([r ** k for k in range(lo, hi + 1)])
        tabs[i] = (lo, T)
    order = np.argsort(comp, kind='stable')
    comp, E, coef, lev = comp[order], E[order], coef[order], lev[order]
    npts = len(pts)
    out = {N: np.zeros((ncomp, npts), LD) for N in orders}
    err = {N: np.zeros((ncomp, npts), LD) for N in orders}     # Kahan compensation
    for s0 in range(0, len(coef), chunk):
        sl = slice(s0, s0 + chunk)
        m = np.repeat(coef[sl][:, None], npts, axis=1)
        for i in used:
            lo, T = tabs[i]
            m *= T[E[sl, i] - lo]
        cc, ll = comp[sl], lev[sl]
        for N in orders:
            sel = ll <= 2 * N
            if not sel.any():
                continue
            mm, ci = m[sel], cc[sel]
            starts = np.flatnonzero(np.r_[True, ci[1:] != ci[:-1]])
            part = np.add.reduceat(mm, starts, axis=0)
            rows = ci[starts]
            # Kahan: out += part
            yv = part - err[N][rows]
            t = out[N][rows] + yv
            err[N][rows] = (t - out[N][rows]) - yv
            out[N][rows] = t
    return out
