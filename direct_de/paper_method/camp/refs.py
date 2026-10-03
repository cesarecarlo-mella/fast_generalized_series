"""exact weight-1,2 references (30 digits) for all campaign point sets.
Euclidean: closed-form GPLs.  Scattering: x+i0, y+i0, principal branches, weight-2 GPLs by quadrature."""
import sys, pickle, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import mpmath as mpm
from exact_w2 import Exact, G, G_quad
mpm.mp.dps = 30
E1, E2 = Exact(1, 30), Exact(2, 30)
name = sys.argv[1]
gf = G if name.startswith('eucl') else G_quad
out = {}
for line in open(os.path.join(os.path.dirname(__file__), name + '.txt')):
    sx, sy = line.split()
    x, y = mpm.mpf(sx), mpm.mpf(sy)
    cv = lambda L: [(str(mpm.mpc(v).real), str(mpm.mpc(v).imag)) for v in L]
    out[(sx, sy)] = (cv(E1(x, y, Gf=gf)), cv(E2(x, y, Gf=gf)))
pickle.dump(out, open(os.path.join(os.path.dirname(__file__), 'ref_%s.pkl' % name), 'wb'))
print(name, len(out))
