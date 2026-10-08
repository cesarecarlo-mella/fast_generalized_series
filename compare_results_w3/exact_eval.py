#!/usr/bin/env python3
# =============================================================================
#  exact_eval.py -- numerical evaluation of the explicit H solution
#  H_exact_sol/H_solution_w3.m  (HW3[i] = {eps^0, eps^1, eps^2, eps^3}
#  coefficients of the canonical master BASISLC[H,"",i], i = 1..371) on a grid.
#
#  The file contains only Log, PolyLog[2|3, rational argument], Pi, Zeta[3],
#  ArcTan and Ti2/Ti3 (Ti_n(u) = Im Li_n(i u), master 157 at eps^3 only).
#  Each entry is translated ONCE into a Python expression and evaluated with
#  mpmath (default 30 digits); all PolyLog/Log/Ti calls are memoised per point,
#  so the ~70 distinct transcendental atoms are computed once per point.
#
#  python3 exact_eval.py --exact H_solution_w3.m --grid grid_N2000.m \
#          --weight 3 --out exact_w3_grid_N2000.npz -n 8
#  -> npz with hi, lo (371, npts): the eps^weight coefficient at each grid point as a
#     double-double hi+lo (real part; the script checks |Im| < 1e-20 * max(1,|Re|)).
# =============================================================================
import argparse
import os
import re
import sys
import time

import mpmath as mpm

HERE = os.path.dirname(os.path.abspath(__file__))

_INT = re.compile(r'(?<![\w.])(\d+)(?![\w.])')


def translate(e):
    """Mathematica expression string -> Python expression string (mpmath)."""
    e = e.replace('\n', ' ')
    e = e.replace('PolyLog[', 'Li(').replace('Log[', 'LOG(').replace('ArcTan[', 'ATAN(')
    e = e.replace('Ti2[', 'TI2(').replace('Ti3[', 'TI3(').replace('Zeta[3]', 'Z3')
    e = e.replace('[', '(').replace(']', ')').replace('^', '**')
    e = re.sub(r'\bPi\b', 'PI', e)
    e = _INT.sub(r'Q(\1)', e)            # every integer literal -> exact mpf
    return e


def load_solution(path):
    """-> dict i -> [4 Mathematica strings] (eps^0..eps^3)."""
    s = open(path).read()
    s = re.sub(r'\(\*.*?\*\)', ' ', s, flags=re.S)
    out = {}
    for m in re.finditer(r'HW3\[(\d+)\]\s*=\s*\{(.*?)\};', s, flags=re.S):
        body, parts, depth, cur = m.group(2), [], 0, ''
        for ch in body:
            if ch in '[({':
                depth += 1
            elif ch in '])}':
                depth -= 1
            if ch == ',' and depth == 0:
                parts.append(cur.strip()); cur = ''
            else:
                cur += ch
        parts.append(cur.strip())
        assert len(parts) == 4, (m.group(1), len(parts))
        out[int(m.group(1))] = parts
    assert sorted(out) == list(range(1, len(out) + 1)), "HW3 indices not 1..n"
    return out


_CODE = None


def compile_weight(path, w):
    sol = load_solution(path)
    return [compile(translate(sol[i][w]), 'HW3[%d][%d]' % (i, w), 'eval') for i in range(1, len(sol) + 1)]


def _namespace(xv, yv):
    memo = {}

    def M(tag, f, *a):
        k = (tag,) + tuple(a)
        v = memo.get(k)
        if v is None:
            v = memo[k] = f(*a)
        return v
    ns = {
        'Q': mpm.mpf, 'PI': mpm.pi, 'Z3': mpm.zeta(3),
        'Li': lambda n, z: M('Li', lambda n_, z_: mpm.polylog(int(n_), z_), n, z),
        'LOG': lambda z: M('log', mpm.log, z),
        'ATAN': lambda z: M('atan', mpm.atan, z),
        'TI2': lambda u: M('ti2', lambda u_: mpm.im(mpm.polylog(2, 1j * u_)), u),
        'TI3': lambda u: M('ti3', lambda u_: mpm.im(mpm.polylog(3, 1j * u_)), u),
        'x': xv, 'y': yv, 'I': mpm.mpc(0, 1),
    }
    return ns


def _mp(v):
    """float, mpf, Fraction / gmpy2 mpq or (num, den) -> mpf at the current precision."""
    if isinstance(v, tuple):
        return mpm.mpf(int(v[0])) / int(v[1])
    if hasattr(v, 'numerator') and hasattr(v, 'denominator') and not isinstance(v, float):
        return mpm.mpf(int(v.numerator)) / int(v.denominator)
    return mpm.mpf(v)


def eval_point(args):
    xv, yv, dps = args
    mpm.mp.dps = dps
    ns = _namespace(_mp(xv), _mp(yv))
    vals, worst_im = [], 0.0
    for code in _CODE:
        v = eval(code, ns)
        v = mpm.mpc(v)
        if abs(v.imag) > 1e-20 * max(1, abs(v.real)):
            worst_im = max(worst_im, float(abs(v.imag)))
        h = float(v.real)
        vals.append((h, float(v.real - h)))       # double-double: hi + lo
    return vals, worst_im


def evaluate(path, w, pts, ncores=1, dps=30):
    """pts: rows with x, y first -- EXACT rationals (Fraction/mpq or (num, den)) are strongly
    preferred: some components are very steep (a 2e-17 shift of x moves them by 6e-15), so the
    series and the exact solution must be evaluated at the very same point.  -> (hi, lo)."""
    import numpy as np
    global _CODE
    _CODE = compile_weight(path, w)
    conv = lambda v: (int(v.numerator), int(v.denominator)) if hasattr(v, 'denominator') and not isinstance(v, float) else float(v)
    jobs = [(conv(p[0]), conv(p[1]), dps) for p in pts]
    if ncores > 1:
        import multiprocessing as mproc
        with mproc.get_context('fork').Pool(ncores) as pool:
            res = pool.map(eval_point, jobs, chunksize=4)
    else:
        res = [eval_point(j) for j in jobs]
    worst = max(r[1] for r in res)
    if worst:
        print("  WARNING: largest non-negligible imaginary part %.2e" % worst)
    hl = np.array([r[0] for r in res])            # (npts, ncomp, 2)
    return hl[:, :, 0].T.copy(), hl[:, :, 1].T.copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exact', required=True)
    ap.add_argument('--grid', required=True)
    ap.add_argument('--weight', type=int, default=3, help='eps power (0..3)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--dps', type=int, default=30)
    ap.add_argument('-n', '--ncores', type=int, default=os.cpu_count())
    a = ap.parse_args()
    sys.path.insert(0, HERE)
    import numpy as np
    import mparse as mp
    grid = mp.parse_file(a.grid)
    pts = [[mp.pconst(v) for v in p] for p in grid]          # exact rationals
    t0 = time.time()
    hi, lo = evaluate(a.exact, a.weight, pts, a.ncores, a.dps)
    np.savez(a.out, hi=hi, lo=lo)
    print("exact eps^%d: %d components x %d points -> %s  (%.0fs)" % (a.weight, hi.shape[0], hi.shape[1], a.out, time.time() - t0))


if __name__ == '__main__':
    main()
