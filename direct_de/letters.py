"""
letters.py -- the H alphabet, numerically.

ray_dlog_series(X, Y, N): along x = s X, y = s Y, each letter is
    d/ds log l_k = m_k / s + sum_i B_k[i] s^i
(m_k = order of vanishing at s = 0).  B_k from an FFT of the regular part on a
circle well inside the nearest singularity.

eval_path(x, y, dx, dy, r_prev): letter dlogs (d/dt log l_k) at complex path
nodes, with Sqrt[xyz] tracked continuously (start: principal, > 0 in the
Euclidean triangle).
"""
import numpy as np
from numpy.polynomial import polynomial as Pl

ACTIVE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16]


def plain_poly_xy():
    """plain letters as {k: dict {(i,j): coeff}} polynomials in x, y"""
    return {1: {(1, 0): 1}, 2: {(0, 1): 1}, 3: {(0, 0): 1, (1, 0): -1, (0, 1): -1}, 4: {(0, 0): 1, (1, 0): -1},
            5: {(0, 0): 1, (0, 1): -1}, 6: {(1, 0): 1, (0, 1): 1}, 9: {(0, 0): 1, (0, 1): -1, (1, 1): -1},
            10: {(0, 0): 1, (1, 0): -1, (1, 1): -1}, 11: {(1, 0): 1, (0, 1): -1, (0, 2): 1},
            12: {(1, 0): -1, (2, 0): 1, (0, 1): 1}, 13: {(0, 0): 1, (1, 0): -2, (2, 0): 1, (0, 1): -1},
            16: {(1, 0): -1, (0, 0): 1, (0, 1): -2, (0, 2): 1}}


PLAIN = plain_poly_xy()


def peval(p, x, y):
    return sum(c * x ** i * y ** j for (i, j), c in p.items())


def pdx(p):
    return {(i - 1, j): c * i for (i, j), c in p.items() if i}


def pdy(p):
    return {(i, j - 1): c * j for (i, j), c in p.items() if j}


def ray_poly(p, X, Y):
    """polynomial in s of p(sX, sY): coefficient array (ascending)"""
    deg = max(i + j for i, j in p)
    out = np.zeros(deg + 1, complex)
    for (i, j), c in p.items():
        out[i + j] += c * X ** i * Y ** j
    return out


def odd_parts(k, x, y, dx, dy):
    """a, a', b, b' for l7 (a = x z), l8 (a = x y); b = x y z (' = d/dparam)."""
    z = 1 - x - y; dz = -dx - dy
    if k == 7:
        a = x * z; da = dx * z + x * dz
    else:
        a = x * y; da = dx * y + x * dy
    b = x * y * z; db = dx * y * z + x * dy * z + x * y * dz
    return a, da, b, db


def odd_dlog(a, da, b, db, r):
    return 1j * (2 * da * b - a * db) / (r * (a * a + b))


def ray_radius(X, Y):
    roots = []
    for k, p in PLAIN.items():
        c = ray_poly(p, X, Y)
        c = np.trim_zeros(c, 'b')
        if len(c) > 1:
            roots += [r for r in Pl.polyroots(c) if abs(r) > 1e-12]
    # odd letters: zeros of a^2 + b and of z
    for poly in (Pl.polymul([0, X], [1, -(X + Y)]), ):
        pass
    xs, ys = np.array([0, X]), np.array([0, Y])
    z = np.array([1, -(X + Y)])
    x_ = np.array([0, X], complex); y_ = np.array([0, Y], complex)
    for a in (Pl.polymul(x_, z), Pl.polymul(x_, y_)):
        b = Pl.polymul(Pl.polymul(x_, y_), z)
        q = Pl.polyadd(Pl.polymul(a, a), b)
        q = np.trim_zeros(q, 'b')
        roots += [r for r in Pl.polyroots(q) if abs(r) > 1e-12]
    roots.append(1 / (X + Y))
    return min(abs(r) for r in roots)


def ray_dlog_series(X, Y, N, M=1024):
    Rmin = ray_radius(X, Y)
    rho = 0.5 * Rmin
    th = 2 * np.pi * np.arange(M) / M
    s = rho * np.exp(1j * th)
    m, B = {}, {}
    for k in ACTIVE:
        if k in PLAIN:
            c = ray_poly(PLAIN[k], X, Y)
            mk = 0
            while abs(c[mk]) < 1e-15:
                mk += 1
            f = Pl.polyval(s, Pl.polyder(c)) / Pl.polyval(s, c) - mk / s
        else:
            x, y = s * X, s * Y
            a, da, b, db = odd_parts(k, x, y, X + 0 * s, Y + 0 * s)
            # r = s Sqrt[X Y] Sqrt[1 - s (X + Y)]  (analytic, principal near s = 0)
            r = s * np.sqrt(X * Y) * np.sqrt(1 - s * (X + Y))
            f = odd_dlog(a, da, b, db, r)
            mk = 0
        coef = np.fft.fft(f) / M                          # c_n rho^n e^{i n th}
        n = np.arange(N)
        B[k] = (coef[:N] / rho ** n)
        m[k] = mk
    return m, B


def eval_nodes(x, y, dx, dy, r_prev):
    """d/dt log l_k at path nodes (arrays), Sqrt[xyz] tracked from r_prev."""
    out = {}
    for k, p in PLAIN.items():
        v = peval(p, x, y)
        out[k] = (peval(pdx(p), x, y) * dx + peval(pdy(p), x, y) * dy) / v
    b = x * y * (1 - x - y)
    r = np.sqrt(b + 0j)
    prev = r_prev
    for i in range(len(r)):
        if abs(r[i] - prev) > abs(r[i] + prev):
            r[i] = -r[i]
        prev = r[i]
    for k in (7, 8):
        a, da, bb, db = odd_parts(k, x, y, dx, dy)
        out[k] = odd_dlog(a, da, bb, db, r)
    return out, r[-1]


def letter_values(x, y, r):
    """the letters themselves (for numerical checks)"""
    v = {k: peval(p, x, y) for k, p in PLAIN.items()}
    for k in (7, 8):
        a, _, _, _ = odd_parts(k, x, y, 0, 0)
        v[k] = (a - 1j * r) / (a + 1j * r)
    return v
