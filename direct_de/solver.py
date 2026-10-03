"""
solver.py -- direct numerical solution of the canonical DE of H (Euclidean
region), started at corner 1 from the shuffle-regularised boundary constants.

    S = Solver()                    # builds J(P0), P0 = s0 (X, Y), from corner 1
    J = S.evaluate(x, y)            # list of 371-vectors, weights 0..6

Steps:
 1. corner.Corner1: exact C_ab = Ax^a Ay^b c / (a! b!) (shuffle regularisation:
    the boundary is the log^0 x log^0 y coefficient), translated to the ray.
 2. Frobenius (log-power series) along the ray x = s X, y = s Y up to s0.
 3. From P0 to (x, y): complex-deformed straight path
        z(t) = P0 + (1 + 4 i delta_k t (1-t)) (P - P0)_k
    and the weight recursion J_w(t) = J_w(0) + int_0^t A J_{w-1} computed with
    Gauss-Legendre panels and a spectral cumulative-integration matrix.
    Refined (panels doubled) until two refinements agree to `tol` (relative, floor 1).
"""
import os
import numpy as np
from numpy.polynomial import legendre as Lg
import scipy.sparse as sps
import scipy.linalg as sla
from common import *
import letters as L
from corner import Corner1, WMAX

HERE = os.path.dirname(os.path.abspath(__file__))
NG = 20
_tau, _wts = Lg.leggauss(NG)
_V = Lg.legvander(_tau, NG - 1)                   # V[i, j] = P_j(tau_i)
_Int = np.zeros((NG, NG))
for j in range(NG):
    c = np.zeros(NG); c[j] = 1
    ci = Lg.legint(c, lbnd=-1)                    # antiderivative of P_j vanishing at -1
    _Int[:, j] = Lg.legval(_tau, ci)
QCUM = _Int @ np.linalg.inv(_V)                   # int_{-1}^{tau_i} f = QCUM @ f(tau)
QEND = _wts                                       # int_{-1}^{1}


class Solver:
    def __init__(self, inputs, bdy_dir, ray=(1.0, 0.7), s0=0.1, N=60):
        self.Aex = load_atilde_exact(os.path.join(inputs, 'atilde.m'))
        A = load_atilde(os.path.join(inputs, 'atilde.m'))
        # BALANCING: atilde has entries up to ~1e5 that cancel; with J = D J', D = diag(sc),
        # the system dJ' = eps D^-1 A D J' has entries <~ 150 (scipy matrix_balance on sum_k |A_k|),
        # which removes most of the double-precision cancellation.
        _, (sc, _) = sla.matrix_balance(sum(np.abs(A[k]) for k in ACTIVE), permute=False, separate=True)
        self.sc = sc
        Ab = {k: A[k] * (1 / sc)[:, None] * sc[None, :] for k in ACTIVE}
        self.A = {k: sps.csr_matrix(Ab[k]) for k in ACTIVE}
        self.corner = Corner1(self.Aex, Ab, os.path.join(inputs, 'sol0.m'),
                              [os.path.join(bdy_dir, 'bdy_c1_w%d.m' % w) for w in range(1, WMAX + 1)])
        self.corner.C = [{ab: v / sc for ab, v in d.items()} for d in self.corner.C]
        X, Y = ray
        self.P0 = (s0 * X, s0 * Y)
        self.J0, _ = self.corner.frobenius(X, Y, s0, N)
        x0, y0 = self.P0
        self.r0 = np.sqrt(x0 * y0 * (1 - x0 - y0) + 0j)          # principal, > 0

    def start_values(self):
        """J at P0 (unbalanced)."""
        return [self.sc * v for v in self.J0]

    def evaluate(self, x, y, delta=(0.1, 0.2), panels=4, tol=1e-10, maxpanels=256, wmax=WMAX):
        prev = None
        M = panels
        while True:
            J = [self.sc * v for v in self._integrate(x, y, delta, M, wmax)]
            if prev is not None:
                err = max(np.max(np.abs(a - b) / np.maximum(1.0, np.abs(a))) for a, b in zip(J, prev))
                if err < tol or M >= maxpanels:
                    return J, dict(panels=M, err=err)
            prev = J
            M *= 2

    def _integrate(self, x, y, delta, M, wmax=WMAX):
        x0, y0 = self.P0
        Dx, Dy = x - x0, y - y0
        edges = np.linspace(0, 1, M + 1)
        J = [v.copy() for v in self.J0]
        r = self.r0
        for a, b in zip(edges[:-1], edges[1:]):
            h = (b - a) / 2
            t = (a + b) / 2 + h * _tau
            # paper's deformation: z_k(t) = x0_k + (t + 4 i delta_k t (1-t)) (x_k - x0_k)
            gx = t + 4j * delta[0] * t * (1 - t); gy = t + 4j * delta[1] * t * (1 - t)
            dgx = 1 + 4j * delta[0] * (1 - 2 * t); dgy = 1 + 4j * delta[1] * (1 - 2 * t)
            xs = x0 + gx * Dx; ys = y0 + gy * Dy
            dxs = dgx * Dx; dys = dgy * Dy
            dl, r = L.eval_nodes(xs, ys, dxs, dys, r)
            # weight recursion on this panel; Jn[w] = values at the nodes
            Jn = [np.repeat(J[0][:, None], NG, axis=1)]
            Jend = [J[0]]
            for w in range(1, wmax + 1):
                Fw = sum(self.A[k] @ (Jn[w - 1] * dl[k][None, :]) for k in ACTIVE)   # (371, NG)
                cum = (Fw @ QCUM.T) * h
                Jn.append(J[w][:, None] + cum)
                Jend.append(J[w] + (Fw @ QEND) * h)
            J = Jend
        return J
