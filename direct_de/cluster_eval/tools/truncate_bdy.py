"""boundary stream (paper layout) truncated to weights 0..W.   usage: truncate_bdy.py IN OUT W"""
import sys
L = open(sys.argv[1]).read().split('\n'); W = int(sys.argv[3])
o, n, ns = map(int, L[2].split())
assert W <= o
body = L[3:3 + (W + 1) * n]; sp = L[3 + (o + 1) * n:3 + (o + 1) * n + ns]
open(sys.argv[2], 'w').write('\n'.join(L[:2] + ['%d %d %d' % (W, n, ns)] + body + sp) + '\n')
