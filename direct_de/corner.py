"""
corner.py -- start of the solution at corner 1 (x, y -> 0) from the
shuffle-regularised boundary constants.

Near the corner (all letters vanishing there: x, y, x+y, x-y(+...)):
    J = sum_{a,b} log^a x log^b y  C_ab(x, y),
    C_ab(0, 0) = Ax^a Ay^b c_{w-a-b} / (a! b!)       (weight w)
with Ax = a_1 (letter x), Ay = a_2 (letter y), c_w = regularised boundary
(coefficient of log^0 x log^0 y; Log x, Log y -> 0).  Checked against the
x^0 y^0 log-part of the old corner-1 series (1e-9 in doubles, exact here).
The letters x+y (a_6) and x - y (a_11 + a_12) annihilate every C_ab(0,0):
they are spurious at the corner.

Along a ray x = s X, y = s Y the DE is 1-dimensional, dJ/ds = eps (R/s + B(s)) J,
and the solution is built order by order (Frobenius) with log^k(s) times
power series; the s-regularised constants (log s -> 0) are
    c^ray_w = sum_{a+b<=w} C_ab log^a X log^b Y .
"""
from math import factorial
import numpy as np
from common import *
import letters as L

WMAX = 6


class Corner1:
    def __init__(self, atilde_exact, A_num, sol0_path, bdy_paths):
        self.A = A_num
        Ax, Ay = atilde_exact[1], atilde_exact[2]
        c = [load_vector_exact(sol0_path)] + [load_vector_exact(p) for p in bdy_paths]
        # C[w][(a, b)] numeric, built exactly
        self.C = []
        for w in range(WMAX + 1):
            d = {}
            for a in range(w + 1):
                for b in range(w + 1 - a):
                    v = c[w - a - b]
                    for _ in range(b): v = matvec_exact(Ay, v)
                    for _ in range(a): v = matvec_exact(Ax, v)
                    d[(a, b)] = to_numeric(v) / (factorial(a) * factorial(b))
            self.C.append(d)

    def ray_constants(self, X, Y):
        lX, lY = np.log(X), np.log(Y)
        return [sum(Cab * lX ** a * lY ** b for (a, b), Cab in self.C[w].items()) for w in range(WMAX + 1)]

    def frobenius(self, X, Y, s0, N=60):
        """J_w(s0) along the ray (x, y) = s (X, Y), w = 0..WMAX, by the regularised
        log-power-series recursion.  Returns list of complex 371-vectors."""
        m, B = L.ray_dlog_series(X, Y, N)          # orders of vanishing, Taylor coeffs of the regular part
        n = self.A[1].shape[0]
        R = sum(m[k] * self.A[k] for k in ACTIVE if m[k])
        Bn = [sum(B[k][i] * self.A[k] for k in ACTIVE) for i in range(N)]   # coefficient matrices of s^i
        cr = self.ray_constants(X, Y)
        # J_w(s) = sum_k log^k(s) sum_i P[w][k][i] s^i
        P = [{0: np.zeros((N + 1, n), complex)}]
        P[0][0][0] = cr[0]
        for w in range(1, WMAX + 1):
            prev = P[w - 1]
            # integrand  F(s) = (R/s + sum_i Bn_i s^i) J_{w-1}  = sum_k log^k sum_i G[k][i] s^(i-1)
            G = {}
            for k, Pk in prev.items():
                g = np.zeros((N + 1, n), complex)
                g[:] = Pk @ R.T                                   # R/s * s^i -> s^(i-1)
                for i in range(1, N + 1):                         # B_j s^j * s^l -> s^(j+l) = s^(i-1)
                    acc = np.zeros(n, complex)
                    for j in range(0, i):
                        acc += Bn[j] @ Pk[i - 1 - j]
                    g[i] += acc
                G[k] = g
            # integrate s^(i-1) log^k s
            new = {}
            for k, g in G.items():
                for i in range(N + 1):
                    if not g[i].any():
                        continue
                    if i == 0:                                    # s^-1 log^k -> log^(k+1)/(k+1)
                        new.setdefault(k + 1, np.zeros((N + 1, n), complex))[0] += g[0] / (k + 1)
                    else:                                         # int s^(i-1) log^k = s^i sum_j (-1)^j k!/(k-j)! log^(k-j) / i^(j+1)
                        for j in range(k + 1):
                            coef = (-1) ** j * factorial(k) / factorial(k - j) / i ** (j + 1)
                            new.setdefault(k - j, np.zeros((N + 1, n), complex))[i] += coef * g[i]
            new.setdefault(0, np.zeros((N + 1, n), complex))[0] += cr[w]
            P.append(new)
        ls = np.log(s0)
        pw = s0 ** np.arange(N + 1)
        return [sum((ls ** k) * (pw @ Pk) for k, Pk in P[w].items()) for w in range(WMAX + 1)], P
