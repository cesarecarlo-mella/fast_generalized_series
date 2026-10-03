#!/usr/bin/env python3
"""
01_chain.py -- STEP 1: the recursion, weights 1..W, modulo ONE prime.

    python3 01_chain.py --workdir dist --corner 1 --k 7 --wmax 6

--k is the prime INDEX (prime = pipe_common.PRIME_TABLE[k]), so a cluster
array job can simply pass its task id.  Jobs are fully independent (no
shared state, no communication) and RESUMABLE: if residues for some weights
already exist for this prime, the chain restarts from the highest one.
Writes <workdir>/residues/J_c<c>_w<w>_k<kkk>.pkl  (one per weight).
"""
import os
# one BLAS thread per job: parallelism comes from running many primes
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import argparse
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blowup_array as ba
import pipe_common as pc

ap = argparse.ArgumentParser()
ap.add_argument('--workdir', required=True)
ap.add_argument('--corner', type=ba.parse_corner, required=True, help='1, 2, 3 or CC2, CC3')
ap.add_argument('--k', type=int, required=True)
ap.add_argument('--wmax', type=int, required=True)
a = ap.parse_args()

S = pc.load_pickle(pc.setup_path(a.workdir, a.corner))
if S['wmax'] < a.wmax:
    sys.exit("[01] setup only has boundaries up to weight %d -- rerun 00_prepare.py" % S['wmax'])
C = S['C']
p = pc.PRIME_TABLE[a.k]
tag = "[01 %s k%03d p=%d]" % (ba.ctag(a.corner), a.k, p)
t0 = time.time()
E = ba.ModEngine(C, p)

w0 = 0
for w in range(a.wmax, 0, -1):
    if os.path.isfile(pc.res_path(a.workdir, a.corner, w, a.k)):
        w0 = w; break
if w0 == a.wmax:
    print(tag, "all weights already present"); sys.exit(0)
J = E.vec_to_arrays(S['sol0']) if w0 == 0 else pc.load_residues(pc.res_path(a.workdir, a.corner, w0, a.k))
print(tag, "order %d/%d, start from weight %d  (setup %.1fs)" % (C.ORDt, C.ORDy, w0, time.time() - t0), flush=True)

for w in range(w0 + 1, a.wmax + 1):
    ts = time.time()
    J, tmul = E.step(J, E.vec_to_arrays(S['bdy'][w]))
    pc.save_residues(J, pc.res_path(a.workdir, a.corner, w, a.k))
    print(tag, "weight %d  %.1fs (multiply %.1fs)  %d tags" % (w, time.time() - ts, tmul, len(J)), flush=True)
import resource
print(tag, "DONE %.1fs   peak memory %.2f GB" % (time.time() - t0,
      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2))
