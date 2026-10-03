"""Export the H-family canonical DE for the C++ solver (paper method, arXiv:2606.30354).

  data/atilde.txt   : dimbasis nletters, then per letter "k nnz" + nnz lines "i j re im" (40 digits)
  data/letters.txt  : plain letters as polynomials in x, y: "k nterms" + lines "i j coeff"
                      (odd letters l7, l8 are hard-coded in h_family.cpp: (a - i r)/(a + i r), r^2 = x y z)
  data/bdy_w2_<tag>.txt : boundary stream in the paper's layout, at x0, weights 0..2, 34 digits,
                      J_1, J_2 from the exact GPL solution, J_0 = sol0, special function r = +Sqrt[x y z].

usage: python3 export_data.py  [x0 y0]
"""
import sys, os
sys.path.insert(0, '..')
import mpmath as mpm
import numpy as np
from common import load_atilde_exact, ACTIVE
from letters import PLAIN
from exact_w2 import Exact
import mparse as mp

INP = '../../ccpipe/inputs'
mpm.mp.dps = 45
os.makedirs('data', exist_ok=True)


def s(q):
    return mpm.nstr(mpm.mpf(q.numerator) / q.denominator, 40, min_fixed=-1, max_fixed=-1) if q else '0'


def export_static():
    A = load_atilde_exact(os.path.join(INP, 'atilde.m'))
    n = 371
    with open('data/atilde.txt', 'w') as f:
        f.write('%d %d\n' % (n, len(ACTIVE)))
        for k in ACTIVE:
            ent = [(ij, v) for ij, v in sorted(A[k].items()) if v[0] or v[1]]
            f.write('%d %d\n' % (k, len(ent)))
            for (i, j), (re, im) in ent:
                f.write('%d %d %s %s\n' % (i, j, s(re), s(im)))
    with open('data/letters.txt', 'w') as f:
        f.write('%d\n' % len(PLAIN))
        for k, p in PLAIN.items():
            f.write('%d %d\n' % (k, len(p)))
            for (i, j), c in p.items():
                f.write('%d %d %d\n' % (i, j, c))


def export_boundary(x0, y0, tag):
    x0, y0 = mpm.mpf(x0), mpm.mpf(y0)
    J0 = [sum(mpm.mpf(v.numerator) / v.denominator for v in p.values()) if p else mpm.mpf(0)
          for p in mp.parse_file(os.path.join(INP, 'sol0.m'))]
    for p in mp.parse_file(os.path.join(INP, 'sol0.m')):
        assert all(k == () for k in p), 'sol0 must be rational'
    J1 = Exact(1, dps=40)(x0, y0)
    J2 = Exact(2, dps=40)(x0, y0)
    r = mpm.sqrt(x0 * y0 * (1 - x0 - y0))
    fmt = lambda v: '%s %s' % (mpm.nstr(mpm.re(v), 34, min_fixed=-1, max_fixed=-1),
                               mpm.nstr(mpm.im(v), 34, min_fixed=-1, max_fixed=-1))
    with open('data/bdy_w2_%s.txt' % tag, 'w') as f:
        f.write('2\n%s %s\n' % (mpm.nstr(x0, 34), mpm.nstr(y0, 34)))
        f.write('2 371 1\n')
        for J in (J0, J1, J2):
            for v in J:
                f.write(fmt(mpm.mpc(v)) + '\n')
        f.write(fmt(mpm.mpc(r)) + '\n')


if __name__ == '__main__':
    export_static()
    x0, y0 = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ('0.3', '0.3')
    export_boundary(x0, y0, '%s_%s' % (x0, y0))
    print('ok')


def export_boundary_frobenius(ray=(1.0, 0.7), s0=0.1, N=60, tag=None):
    """weights 0..6 at P0 = s0*ray from the corner-1 Frobenius series (double precision)."""
    from solver import Solver
    S = Solver(INP, '/mnt/user-data/uploads/series.all/blow_up_array/pipeline/inputs', ray=ray, s0=s0, N=N)
    x0, y0 = S.P0
    J = [np.asarray(S.J0[w]) * S.sc for w in range(7)]
    r = np.sqrt(x0 * y0 * (1 - x0 - y0))
    tag = tag or 'frob_%g_%g' % (x0, y0)
    with open('data/bdy_w6_%s.txt' % tag, 'w') as f:
        f.write('2\n%.17e %.17e\n6 371 1\n' % (x0, y0))
        for Jw in J:
            for v in Jw:
                f.write('%.17e %.17e\n' % (v.real, v.imag))
        f.write('%.17e 0\n' % r)
    return 'data/bdy_w6_%s.txt' % tag, (x0, y0), J


def const_mp(key):
    v = mpm.mpf(1)
    for at, e in key:
        if at == 'Pi': v *= mpm.pi ** (e // 2)
        elif at.startswith('Zeta['): v *= mpm.zeta(int(at[5:-1])) ** (e // 2)
        elif at == 'Log[2]': v *= mpm.log(2) ** (e // 2)
        else: raise ValueError(at)
    return v


def export_ray_constants(X, Y, bdy_dir='/mnt/user-data/uploads/series.all/blow_up_array/pipeline/inputs', wmax=6,
                         prefix='bdy_c1', tag=''):
    """c^ray_w = sum_{a+b<=w} Ax^a Ay^b c_{w-a-b} log^a X log^b Y /(a! b!), 40 digits.
    X, Y: strings (exact decimals).  -> data/raycst<tag>_X_Y.txt : wmax, then (wmax+1)*371 lines 're im'.
    Continued region (corner CC1, x = -w/v, y = -u/v): prefix='bdy_CC1', X, Y < 0; the constants are
    regularised in Log[w], Log[u] (w ~ -x, u ~ -y at the corner), so the logs are log|X|, log|Y|."""
    from math import factorial
    from common import load_vector_exact, matvec_exact
    A = load_atilde_exact(os.path.join(INP, 'atilde.m'))
    Ax, Ay = A[1], A[2]
    c = [load_vector_exact(os.path.join(INP, 'sol0.m'))] + \
        [load_vector_exact(os.path.join(bdy_dir, '%s_w%d.m' % (prefix, w))) for w in range(1, wmax + 1)]
    lX, lY = mpm.log(abs(mpm.mpf(X))), mpm.log(abs(mpm.mpf(Y)))
    tonum = lambda vec: _tonum(vec)
    out = []
    for w in range(wmax + 1):
        tot = [mpm.mpc(0)] * 371
        for a in range(w + 1):
            for b in range(w + 1 - a):
                v = c[w - a - b]
                for _ in range(b): v = matvec_exact(Ay, v)
                for _ in range(a): v = matvec_exact(Ax, v)
                fac = lX ** a * lY ** b / (factorial(a) * factorial(b))
                for ck, comp in v.items():
                    cv = const_mp(ck) * fac
                    for j, (re, im) in comp.items():
                        tot[j] += cv * mpm.mpc(mpm.mpf(re.numerator) / re.denominator, mpm.mpf(im.numerator) / im.denominator)
        out.append(tot)
    fn = 'data/raycst%s_%s_%s.txt' % (tag, X, Y)
    with open(fn, 'w') as f:
        f.write('%d\n' % wmax)
        for tot in out:
            for v in tot:
                f.write('%s %s\n' % (mpm.nstr(v.real, 36, min_fixed=-1, max_fixed=-1), mpm.nstr(v.imag, 36, min_fixed=-1, max_fixed=-1)))
    return fn
