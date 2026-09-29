#!/usr/bin/env python3
"""
00_prepare.py -- STEP 0: parse all inputs ONCE for one corner at one order
and store them (exact) in <workdir>/setup_c<c>.pkl for the prime jobs.

    python3 00_prepare.py --corner 1 --ordt 15 --ordy 45 --wmax 6 \
        --inputs ../../blow_up_form_2/pipeline/dist \
        --letters ../../blow_up_form_2/pipeline/master_letters --workdir dist

--inputs  : atilde.m, sol0.m, bdy_c<c>_w<w>.m   (a missing bdy file is an
            error unless --allow-missing-bdy, then it is a zero boundary)
--letters : dlog_dt_c<c>.m, dlog_dy_c<c>.m at ANY order >= (ordt, ordy):
            terms above the requested order are dropped on load (exact --
            same as select_letters.wl / chop_letters.wl), so the master
            letters can be used directly.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blowup_array as ba
import pipe_common as pc

ap = argparse.ArgumentParser()
ap.add_argument('--corner', type=int, required=True)
ap.add_argument('--ordt', type=int, required=True)
ap.add_argument('--ordy', type=int, required=True)
ap.add_argument('--wmax', type=int, required=True)
ap.add_argument('--inputs', required=True)
ap.add_argument('--letters', default=None)
ap.add_argument('--workdir', required=True)
ap.add_argument('--allow-missing-bdy', action='store_true',
                help='treat a missing bdy_c<c>_w<w>.m as a zero boundary instead of stopping')
a = ap.parse_args()
letters = a.letters or a.inputs
os.makedirs(a.workdir, exist_ok=True)

t0 = time.time()
# letters and atilde may live in different dirs: Corner reads both from one
# dir, so point it at the letters dir and feed atilde separately if needed
if os.path.abspath(letters) != os.path.abspath(a.inputs):
    import shutil, tempfile
    tmp = tempfile.mkdtemp(prefix='bua_')
    for f in ('atilde.m',):
        os.symlink(os.path.abspath(os.path.join(a.inputs, f)), os.path.join(tmp, f))
    for f in ('dlog_dt_c%d.m' % a.corner, 'dlog_dy_c%d.m' % a.corner):
        os.symlink(os.path.abspath(os.path.join(letters, f)), os.path.join(tmp, f))
    C = ba.Corner(tmp, a.corner, a.ordt, a.ordy)
    shutil.rmtree(tmp)
else:
    C = ba.Corner(a.inputs, a.corner, a.ordt, a.ordy)

sol0 = C.load_vector(os.path.join(a.inputs, 'sol0.m'))
bdy = {}
for w in range(1, a.wmax + 1):
    f = os.path.join(a.inputs, 'bdy_c%d_w%d.m' % (a.corner, w))
    if os.path.isfile(f):
        bdy[w] = C.load_vector(f)
    elif a.allow_missing_bdy:
        print("[00] WARNING: no %s -- using zero boundary at weight %d" % (f, w))
        bdy[w] = {}
    else:
        sys.exit("[00] ERROR: missing boundary file %s (add it to the inputs, or pass "
                 "--allow-missing-bdy to use a zero boundary)" % f)
pc.save_pickle({'C': C, 'sol0': sol0, 'bdy': bdy, 'ordt': a.ordt, 'ordy': a.ordy,
                'wmax': a.wmax}, pc.setup_path(a.workdir, a.corner))
print("[00 c%d] setup written (order %d/%d, weights <= %d) in %.1fs"
      % (a.corner, a.ordt, a.ordy, a.wmax, time.time() - t0))
