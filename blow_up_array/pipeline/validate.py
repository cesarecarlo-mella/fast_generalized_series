#!/usr/bin/env python3
"""
validate.py -- exact, coefficient-by-coefficient comparison of this
pipeline's results (rec_c<c>_w<w>.pkl in --workdir) against
  --ref-exact DIR : exact_c<c>_w<w>.m from the FORM/Mathematica pipelines at
                    ANY order >= ours (truncated to ours on load), and/or
  --ref-phys DIR  : phys_c<c>_w<w>_ord<N>.m (e.g. cluster_results), compared
                    in the window SMALL <= ORDt, BIG <= ORDy - 2 ORDt (the
                    largest window our (ORDt,ORDy) rectangle fully determines)

    python3 validate.py --workdir dist -c 1 -w 1,2,3,4 --ref-exact ../../blow_up_form/pipeline/dist
    python3 validate.py --workdir dist -c 1 -w 5,6 --ref-phys ../../cluster_results/blow_up_form_2_order_20/out --phys-ord 20
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blowup_array as ba
import pipe_common as pc
import phys

ap = argparse.ArgumentParser()
ap.add_argument('--workdir', required=True)
ap.add_argument('-c', '--corners', default='1')
ap.add_argument('-w', '--weights', default='1,2')
ap.add_argument('--ref-exact', default=None)
ap.add_argument('--ref-phys', default=None)
ap.add_argument('--phys-ord', type=int, default=20)
a = ap.parse_args()

allok = True
for c in [int(x) for x in a.corners.split(',')]:
    C = pc.load_pickle(pc.setup_path(a.workdir, c))['C']
    X, Y = C.ORDt, C.ORDy - 2 * C.ORDt
    for w in [int(x) for x in a.weights.split(',')]:
        mine = pc.load_pickle(pc.rec_path(a.workdir, c, w))
        if a.ref_exact:
            f = os.path.join(a.ref_exact, 'exact_c%d_w%d.m' % (c, w))
            t0 = time.time()
            ref = ba.exact_to_slots(C, C.load_vector(f))
            onlyr, onlym, diff = phys.compare(mine, ref)
            bad = sorted(set(k[1] + 1 for k in onlyr + onlym + diff))
            allok &= not bad
            print("c%d w%d exact: ref %d, mine %d | missing %d extra %d different %d | bad comps %s  (%.0fs)"
                  % (c, w, len(ref), len(mine), len(onlyr), len(onlym), len(diff), bad[:15], time.time() - t0), flush=True)
        if a.ref_phys:
            f = os.path.join(a.ref_phys, 'phys_c%d_w%d_ord%d.m' % (c, w, a.phys_ord))
            t0 = time.time()
            ref, kept, total = phys.load_phys_ref(f, c, X, Y)
            mp_ = phys.slots_to_phys(C, mine, X, Y)
            onlyr, onlym, diff = phys.compare(mp_, ref)
            bad = sorted(set(k[0] + 1 for k in onlyr + onlym + diff))
            allok &= not bad
            print("c%d w%d phys (SMALL<=%d, BIG<=%d): ref %d, mine %d | missing %d extra %d different %d | bad comps %s  (%.0fs)"
                  % (c, w, X, Y, len(ref), len(mp_), len(onlyr), len(onlym), len(diff), bad[:15], time.time() - t0), flush=True)
print("ALL IDENTICAL" if allok else "MISMATCHES FOUND")
sys.exit(0 if allok else 1)
