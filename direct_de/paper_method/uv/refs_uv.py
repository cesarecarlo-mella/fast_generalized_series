import sys, pickle, os
sys.path.insert(0, '../..')
import mpmath as mpm
from exact_w2 import Exact, G_quad
mpm.mp.dps = 30
E1, E2 = Exact(1, 30), Exact(2, 30)
out = {}
cv = lambda L: [(str(mpm.mpc(v).real), str(mpm.mpc(v).imag)) for v in L]
for line in open('lat.txt'):
    su, sv = line.split(); u, v = mpm.mpf(su), mpm.mpf(sv)
    x, y = -(1 - u - v) / v, -u / v
    out[(su, sv)] = (cv(E1(x, y, Gf=G_quad)), cv(E2(x, y, Gf=G_quad)))
pickle.dump(out, open('ref_uv_lat.pkl', 'wb')); print(len(out))
