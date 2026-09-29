#!/usr/bin/env python3
"""
02_reconstruct.py -- STEP 2: exact rationals for one weight from all primes
available in <workdir>/residues, verified against 2 extra primes.

    python3 02_reconstruct.py --workdir dist --corner 0 --weight 4 --jobs 16

exit 0 : written  <workdir>/exact_c<c>_w<w>.m   (series in t = x - 1/2, y, Log[y])
         FORM pipeline's file -> 03_to_physical_vars.wl etc. work unchanged)
         and      <workdir>/rec_c<c>_w<w>.pkl   (exact slots)
exit 3 : not enough primes yet (some coefficients unverified) -- run more
         01_chain.py jobs (higher --k) and call this again.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blowup_array as ba
import pipe_common as pc

ap = argparse.ArgumentParser()
ap.add_argument('--workdir', required=True)
ap.add_argument('--corner', type=int, required=True)
ap.add_argument('--weight', type=int, required=True)
ap.add_argument('--jobs', type=int, default=os.cpu_count())
ap.add_argument('--no-exact', action='store_true', help='skip writing exact_c<c>_w<w>.m')
a = ap.parse_args()

c, w = a.corner, a.weight
tag = "[02 c%d w%d]" % (c, w)
ks = pc.primes_with_residues(a.workdir, c, w)
if len(ks) < 3:
    print(tag, "only %d primes available -- need at least 3" % len(ks)); sys.exit(3)
t0 = time.time()
results = [pc.load_residues(pc.res_path(a.workdir, c, w, k)) for k in ks]
primes = [pc.PRIME_TABLE[k] for k in ks]
out, bad, n2 = ba.reconstruct(results, primes, jobs=a.jobs)
del results
print(tag, "%d primes (+%d verification): %d coefficients, %d via common denominator, %d unverified  (%.1fs)"
      % (len(ks) - ba.NV, ba.NV, len(out), n2 - bad, bad, time.time() - t0), flush=True)
if bad:
    print(tag, "NEED_MORE_PRIMES"); sys.exit(3)

S = pc.load_pickle(pc.setup_path(a.workdir, c))
C = S['C']
pc.save_pickle(out, pc.rec_path(a.workdir, c, w))
if not a.no_exact:
    t1 = time.time()
    pc.write_exact(C, out, pc.exact_path(a.workdir, c, w))
    print(tag, "wrote exact_c%d_w%d.m  (%.1fs)" % (c, w, time.time() - t1), flush=True)
print(tag, "DONE %.1fs" % (time.time() - t0))
