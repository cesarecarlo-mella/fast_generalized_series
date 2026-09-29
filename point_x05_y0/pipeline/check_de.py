#!/usr/bin/env python3
"""
check_de.py -- integrability test of the letters at this expansion point
(independent of the boundary constants).

At every weight the recursion computes
    intt   = Int[ (atilde /. dlog_dt) . J , t ]
    tointy = -D[intt, y] + (atilde /. dlog_dy) . J
For a correct alphabet (d A = A ^ A holds only if every letter expansion is
right) tointy must NOT depend on t: all its t^k (k >= 1) coefficients vanish,
inside the region the truncation determines exactly (y^e with e <= ORDy - 1).
This is checked modulo one prime, weight by weight, with a zero boundary
(any boundary gives a solution of the same DE, so this tests the letters and
atilde, not the constants).

    python3 check_de.py --ordt 8 --ordy 8 --wmax 4 [--inputs inputs] [--corner 0]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import blowup_array as ba

ap = argparse.ArgumentParser()
ap.add_argument('--ordt', type=int, default=8)
ap.add_argument('--ordy', type=int, default=8)
ap.add_argument('--wmax', type=int, default=4)
ap.add_argument('--corner', type=int, default=0)
ap.add_argument('--inputs', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inputs'))
a = ap.parse_args()

C = ba.Corner(a.inputs, a.corner, a.ordt, a.ordy)
p = ba.gen_primes(1)[0]
E = ba.ModEngine(C, p)
J = E.vec_to_arrays(C.load_vector(os.path.join(a.inputs, 'sol0.m')))
ok = True
for w in range(1, a.wmax + 1):
    t0 = time.time()
    tags = sorted(t for t, x in J.items() if x.any())
    X = np.stack([J[t] for t in tags])
    Xb = np.where(X > p // 2, X - p, X).astype(np.float64)
    X2 = np.ascontiguousarray(Xb.transpose(1, 0, 2, 3).reshape(C.ncomp, len(tags) * C.S))
    multT = E.multiply(J, tags, X2, 't'); E.Tcache.clear()
    multY = E.multiply(J, tags, X2, 'y'); E.Tcache.clear()
    intt = E.integrate(multT, 1)
    tointy = E.dy(intt)
    for tg, x in multY.items():
        b = tointy.get(tg)
        if b is None:
            tointy[tg] = x
        else:
            b += x; np.remainder(b, p, out=b)
    bad = checked = 0
    logt = sum(1 for tg, x in {**intt, **tointy}.items() if tg[2] and x.any())
    for (pt, py, la, lb, ck), x in tointy.items():
        # grid index iy <-> exponent (iy - 1) + py/2 ; keep e <= ORDy - 1
        iymax = C.ORDy - 1 + 1 - (1 if py else 0)
        region = x[:, 2:, :iymax + 1]                 # t^k, k >= 1
        checked += region.size
        bad += int(np.count_nonzero(region))
    print("weight %d: t-dependent coefficients of tointy in the exact region: %d of %d  "
          "(Log[t] tags: %d)   %.1fs" % (w, bad, checked, logt, time.time() - t0), flush=True)
    ok &= (bad == 0 and logt == 0)
    J, _ = E.step(J, {})
print("INTEGRABLE (letters + atilde consistent at this point)" if ok else "!! NOT integrable -- a letter expansion is wrong")
sys.exit(0 if ok else 1)
