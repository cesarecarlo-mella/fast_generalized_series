#!/usr/bin/env python3
"""
exact_continued.py -- the EXACT weight-0,1,2 solution of H, analytically
continued to the region x, y < 0, z > 1 (x = -w/v, y = -u/v, z = 1/v).

Starting point P0 = (rho, rho) in the Euclidean triangle: the exact GPL
solution (compare_results/exact_solutions/eps_{1,2}.m) evaluated there with
mpmath.  It is then continued along a path in C^2 by integrating the canonical
DE exactly at the level of iterated integrals (weight <= 2 needs only single
and double integrals of dlog forms):

    f1(P) = f1(P0) + sum_k A_k f0 dL_k
    f2(P) = f2(P0) + sum_k A_k f1(P0) dL_k + sum_{j,k} A_k A_j f0 I_jk
    dL_k = int dlog l_k ,  I_jk = int (L_j - L_j(P0)) dlog l_k

Path:   (1) x = y = rho e^{i pi s}, s: 0 -> 1     (both phases +pi near the origin:
            Log x, Log y -> log|.| + i pi, Sqrt[xyz] -> -Sqrt|xyz|, exactly the
            continuation used for the CC boundaries)
        (2) straight line (-rho, -rho) -> (x_T, y_T) with an imaginary bump
            +/- i*delta*s(1-s) (BUMP = +1 / -1): equal results for both bumps
            means no singularity of the masters between the two paths.
Logs are tracked continuously along the path (phase unwrapping on a dense
Gauss-Legendre grid); the square root r = Sqrt[xyz] of l7, l8 is tracked
continuously as well.
"""
import os
import re
import sys

import mpmath as mpm
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))
import mparse as mp

ACTIVE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16]
X, Y = sp.symbols('x y')
PLAIN = {1: X, 2: Y, 3: 1 - X - Y, 4: 1 - X, 5: 1 - Y, 6: X + Y, 9: 1 - (1 + X) * Y, 10: 1 - X * (1 + Y),
         11: X - Y + Y**2, 12: -X + X**2 + Y, 13: (X - 1)**2 - Y, 16: -X + (1 - Y)**2}
ODD_A = {7: X * (1 - X - Y), 8: X * Y}
B_RAD = X * Y * (1 - X - Y)
_lam = lambda e: sp.lambdify((X, Y), e, 'numpy')
PL = {k: (_lam(p), _lam(sp.diff(p, X)), _lam(sp.diff(p, Y))) for k, p in PLAIN.items()}
OD = {k: (_lam(a), _lam(sp.diff(a, X)), _lam(sp.diff(a, Y))) for k, a in ODD_A.items()}
BR = (_lam(B_RAD), _lam(sp.diff(B_RAD, X)), _lam(sp.diff(B_RAD, Y)))


def _arr(f, x, y):
    v = f(x, y)
    return np.broadcast_to(np.asarray(v, complex), np.shape(x)).copy()


# ------------------------------------------------------------- DE data
def load_atilde(path):
    A = mp.parse_file(path)
    n = len(A)
    M = {k: np.zeros((n, n), complex) for k in ACTIVE}
    for i, row in enumerate(A):
        for j, e in enumerate(row):
            for key, v in e.items():
                d = dict(key)
                im = d.pop('I', 0)
                (lt, _), = d.items()
                k = int(lt[2:-1])
                M[k][i, j] += float(v) * (1j if im else 1)
    return M


def load_const_vector(path):
    vec = mp.parse_file(path)
    out = []
    for p in vec:
        s = 0
        for key, v in p.items():
            t = mpm.mpf(v.numerator) / v.denominator
            for at, e in key:
                if at == 'Pi':
                    t *= mpm.pi ** (e // 2)
                elif at == 'I':
                    t *= mpm.mpc(0, 1) ** (e // 2)
                else:
                    raise ValueError(at)
            s += t
        out.append(complex(s))
    return np.array(out)


# ------------------------------------------------------------- exact GPLs at P0
def G(args, z):
    if all(a == 0 for a in args):
        return mpm.log(z) ** len(args) / mpm.factorial(len(args))
    if len(args) == 1:
        return mpm.log(1 - z / args[0])
    a1, rest = args[0], args[1:]
    return mpm.quad(lambda s: G(rest, s) / (s - a1), [0, z])


def eval_exact(path, x0, y0):
    mpm.mp.dps = 30
    vec = mp.parse_file(path)
    cache = {}
    env = {'x': mpm.mpf(x0), 'y': mpm.mpf(y0)}

    def atom(at):
        if at in cache:
            return cache[at]
        m = re.fullmatch(r'G\[(.*)\]', at)
        if m:
            parts = [p.strip() for p in m.group(1).split(',')]
            args = [mpm.mpf(eval(p, {}, env)) for p in parts[:-1]]
            v = G(args, env[parts[-1]])
        elif at == 'Pi':
            v = mpm.pi
        else:
            raise ValueError(at)
        cache[at] = v
        return v

    out = []
    for p in vec:
        s = mpm.mpf(0)
        for key, v in p.items():
            t = mpm.mpf(v.numerator) / v.denominator
            for at, e in key:
                if at == 'I':
                    t *= mpm.mpc(0, 1) ** (e // 2)
                else:
                    t *= atom(at) ** (e // 2)
            s += t
        out.append(complex(s))
    return np.array(out)


# ------------------------------------------------------------- path machinery
GL_N = 16
_gx, _gw = np.polynomial.legendre.leggauss(GL_N)


def nodes(M):
    """Gauss-Legendre nodes/weights on [0,1] split into M pieces (sorted),
    plus the two end points with weight 0 (so logs are anchored exactly)."""
    edges = np.linspace(0, 1, M + 1)
    s = ((edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * _gx[None, :]).ravel()
    w = ((edges[1:, None] - edges[:-1, None]) / 2 * _gw[None, :]).ravel()
    return np.concatenate([[0.0], s, [1.0]]), np.concatenate([[0.0], w, [0.0]])


def letters_along(x, y, dx, dy, r_prev):
    """values (for unwrapping) and d/ds log of all letters at the nodes.
    r = Sqrt[xyz] tracked continuously, starting from r_prev."""
    vals, ders = {}, {}
    for k, (f, fx, fy) in PL.items():
        v = _arr(f, x, y)
        vals[k] = v
        ders[k] = (_arr(fx, x, y) * dx + _arr(fy, x, y) * dy) / v
    b, bx, by = (_arr(g, x, y) for g in BR)
    r = np.sqrt(b)
    prev = r_prev
    for i in range(len(r)):                      # continuous branch of the root
        if abs(r[i] - prev) > abs(r[i] + prev):
            r[i] = -r[i]
        prev = r[i]
    db = bx * dx + by * dy
    for k, (fa, fax, fay) in OD.items():
        a = _arr(fa, x, y)
        da = _arr(fax, x, y) * dx + _arr(fay, x, y) * dy
        vals[k] = (a - 1j * r) / (a + 1j * r)
        ders[k] = 1j * (2 * da * b - a * db) / (r * (a * a + b))
    return vals, ders, r[-1]


def continue_logs(vals_seq, L0):
    """continuous log along consecutive nodes, starting from L0 (complex)."""
    ang = np.unwrap(np.angle(vals_seq))
    ang = ang - ang[0] + np.angle(vals_seq[0])
    L = np.log(np.abs(vals_seq)) + 1j * ang
    # align the start with L0 (L0 may differ by 2 pi i from principal)
    return L + (L0 - np.log(vals_seq[0]))


class Continuer:
    def __init__(self, atilde_path, sol0_path, eps1_path, eps2_path, rho=0.1, M=150):
        self.A = load_atilde(atilde_path)
        self.f0 = load_const_vector(sol0_path)
        self.rho, self.M = rho, M
        self.f1P0 = eval_exact(eps1_path, rho, rho)
        self.f2P0 = eval_exact(eps2_path, rho, rho)
        self.V = {k: self.A[k] @ self.f0 for k in ACTIVE}          # A_k f0
        self.W = {k: self.A[k] @ self.f1P0 for k in ACTIVE}        # A_k f1(P0)
        self.AV = {(j, k): self.A[k] @ self.V[j] for j in ACTIVE for k in ACTIVE}
        # segment 1 (common to every target): rotation of both phases by +pi
        s, w = nodes(M)
        x = rho * np.exp(1j * np.pi * s); dx = 1j * np.pi * x
        L0 = {}
        x0 = np.array([rho + 0j]); dummy = np.zeros(1)
        v0, _, r0 = letters_along(x0, x0, dummy, dummy, np.sqrt(complex(rho * rho * (1 - 2 * rho))))
        for k in ACTIVE:
            L0[k] = np.log(v0[k][0])
        self.r_start = np.sqrt(complex(rho * rho * (1 - 2 * rho)))     # principal, > 0 at P0
        vals, ders, r_end = letters_along(x, x, dx, dx, self.r_start)
        self.seg1 = dict(w=w, vals=vals, ders=ders, L0=L0, r_end=r_end)

    def evaluate(self, xT, yT, bump=+1, delta=0.05, tri=True):
        """f0, f1, f2 at (xT, yT) (complex 371-vectors).
        tri=True: segment 2 is a straight line in the TRIANGLE coordinates
        (u, w) (v = 1-u-w, x = -w/v, y = -u/v) with an imaginary bump
        +/- i*delta*s(1-s) on u and w -- well conditioned also for targets
        with v -> 0 (|x|, |y| large).  tri=False: straight line in (x, y)."""
        rho = self.rho
        s, w = nodes(self.M)
        if tri:
            z0 = 1 + 2 * rho; u1 = w1 = rho / z0               # (u, w) of (-rho, -rho)
            zT = 1 - xT - yT; uT, wT = -yT / zT, -xT / zT
            bb = bump * 1j * delta * s * (1 - s); dbb = bump * 1j * delta * (1 - 2 * s)
            uu = u1 + (uT - u1) * s + bb; ww = w1 + (wT - w1) * s + bb
            du = (uT - u1) + dbb; dw = (wT - w1) + dbb
            vv = 1 - uu - ww; dv = -du - dw
            x = -ww / vv; y = -uu / vv
            dx = -(dw * vv - ww * dv) / vv ** 2
            dy = -(du * vv - uu * dv) / vv ** 2
        else:
            bx = bump * 1j * delta * s * (1 - s)
            x = -rho + (xT + rho) * s + bx
            y = -rho + (yT + rho) * s + bx
            dbx = bump * 1j * delta * (1 - 2 * s)
            dx = (xT + rho) + dbx
            dy = (yT + rho) + dbx
        S1 = self.seg1
        vals2, ders2, r_end = letters_along(x, y, dx, dy, S1['r_end'])
        dL, I = {}, {}
        Lnodes = {}
        for k in ACTIVE:
            seq = np.concatenate([S1['vals'][k], vals2[k]])
            L = continue_logs(seq, S1['L0'][k])
            Lnodes[k] = L - S1['L0'][k]                       # L_k - L_k(P0) at all nodes
        wts = np.concatenate([S1['w'], w])
        ders = {k: np.concatenate([S1['ders'][k], ders2[k]]) for k in ACTIVE}
        for k in ACTIVE:
            dL[k] = np.sum(wts * ders[k])
        for j in ACTIVE:
            for k in ACTIVE:
                I[(j, k)] = np.sum(wts * Lnodes[j] * ders[k])
        f1 = self.f1P0 + sum(self.V[k] * dL[k] for k in ACTIVE)
        f2 = self.f2P0 + sum(self.W[k] * dL[k] for k in ACTIVE) + sum(self.AV[(j, k)] * I[(j, k)] for j in ACTIVE for k in ACTIVE)
        # consistency: the endpoint log from unwrapping vs integrated dlog
        return self.f0, f1, f2, dict(dL=dL, Lend={k: Lnodes[k][-1] for k in ACTIVE}, r_end=r_end)


if __name__ == '__main__':
    import time
    base = '/tmp/claude-0/-home-claude/62cee7e5-47c1-5a82-adde-578e1e20317b/scratchpad/ccpipe'
    ex = '/mnt/user-data/uploads/series.all/compare_results/exact_solutions'
    t0 = time.time()
    C = Continuer(base + '/inputs/atilde.m', base + '/inputs/sol0.m', ex + '/eps_1.m', ex + '/eps_2.m')
    print("setup %.1fs" % (time.time() - t0))
    for (xT, yT) in ((-0.3, -0.2), (-2.0, -0.5), (-5.0, -4.0)):
        for b in (+1, -1):
            f0, f1, f2, info = C.evaluate(xT, yT, bump=b)
            dmax = max(abs(info['dL'][k] - info['Lend'][k]) for k in ACTIVE)
            print("(%g,%g) bump %+d: |integrated dlog - unwrapped log| max %.1e ; f2[7] = %s" % (xT, yT, b, dmax, f2[6]))


def evaluate_adaptive(C, xT, yT, tol=1e-12, Mmax=64000, bump=+1):
    """evaluate with increasing node count until the integrated dlogs agree with
    the unwrapped logs (to tol) and f2 is stable between two refinements."""
    M = C.M
    prev = None
    delta = 0.15
    while True:
        C.M = M
        f0, f1, f2, info = C.evaluate(xT, yT, bump=bump, delta=delta)
        mm = max(abs(info['dL'][k] - info['Lend'][k]) for k in ACTIVE)
        stable = prev is not None and np.max(np.abs(f2 - prev) / np.maximum(np.abs(f2), 1e-3)) < 1e-11
        if (mm < tol and stable) or M >= Mmax:
            C.M = 1000
            return f0, f1, f2, dict(M=M, dlog_err=mm, ok=(mm < tol and stable))
        prev = f2
        M *= 2
