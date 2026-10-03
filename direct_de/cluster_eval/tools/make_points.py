"""lattice points (exact decimal strings, 30 digits) and round-robin chunks.
   eucl: x = i/N, y = j/N            (i, j >= 1, i + j <= N - 1)
   scat: u = i/N, v = j/N, w = 1-u-v (same lattice; x = -w/v, y = -u/v, z = 1/v)
A lattice point lying EXACTLY on the zero set of a letter (connection singular there, e.g. l13 = 0 at
(1/3, 4/9) for N = 9) is moved by 1e-6/N along the first coordinate (reported); the end point of the
path must not be a singular point of the DE.  N = 65 has no such points.
usage: make_points.py OUTDIR NLAT NCHUNK REGION..."""
import os, sys
from fractions import Fraction as F
import mpmath as mpm
mpm.mp.dps = 40
PLAIN = {1: {(1, 0): 1}, 2: {(0, 1): 1}, 3: {(0, 0): 1, (1, 0): -1, (0, 1): -1}, 4: {(0, 0): 1, (1, 0): -1},
         5: {(0, 0): 1, (0, 1): -1}, 6: {(1, 0): 1, (0, 1): 1}, 9: {(0, 0): 1, (0, 1): -1, (1, 1): -1},
         10: {(0, 0): 1, (1, 0): -1, (1, 1): -1}, 11: {(1, 0): 1, (0, 1): -1, (0, 2): 1},
         12: {(1, 0): -1, (2, 0): 1, (0, 1): 1}, 13: {(0, 0): 1, (1, 0): -2, (2, 0): 1, (0, 1): -1},
         16: {(1, 0): -1, (0, 0): 1, (0, 1): -2, (0, 2): 1}}


def singular(x, y):
    z = 1 - x - y; b = x * y * z
    ks = [k for k, p in PLAIN.items() if sum(c * x ** i * y ** j for (i, j), c in p.items()) == 0]
    ks += [k for k, a in ((7, x * z), (8, x * y)) if a * a + b == 0]
    return ks


out, N, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
os.makedirs(out, exist_ok=True)
lat = [(i, j) for i in range(1, N) for j in range(1, N) if i + j <= N - 1]
for reg in sys.argv[4:]:
    L = []
    for i, j in lat:
        a, b = F(i, N), F(j, N)
        x, y = (a, b) if reg == 'eucl' else (-(1 - a - b) / b, -a / b)
        ks = singular(x, y)
        sa = mpm.nstr(mpm.mpf(i) / N, 30)
        if ks:
            sa = mpm.nstr(mpm.mpf(i) / N + mpm.mpf('1e-6') / N, 30)
            print('  %s point (%d, %d)/%d on letter(s) %s -> first coordinate moved to %s' % (reg, i, j, N, ks, sa))
        L.append('%s %s' % (sa, mpm.nstr(mpm.mpf(j) / N, 30)))
    open(os.path.join(out, '%s.txt' % reg), 'w').write('\n'.join(L) + '\n')
    for k in range(K):
        open(os.path.join(out, '%s_%d.txt' % (reg, k)), 'w').write(''.join(l + '\n' for l in L[k::K]))
    print(reg, len(L), 'points,', K, 'chunks')
