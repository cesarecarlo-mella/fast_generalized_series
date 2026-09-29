#!/usr/bin/env python3
# =============================================================================
#  run_all_families_bernoulli.py -- Bernoulli (conformal) resummation for every
#  family in all_families_array/, all corners, weights 1..W.
#
#  For each family:
#    1. if phys_c<c>_w<w>_ord<ORD>.m is missing but the recursion finished
#       (rec_c<c>_w<w>.pkl exists), write it first from the exact result
#       (only the variable change -- nothing is recomputed)
#    2. run code/run_bernoulli.py (Mathematica) on it
#  Existing ber_ files are skipped (use --redo to recompute them).
#
#  USAGE  (run with the venv python -- step 1 needs numpy/gmpy2)
#     venv/bin/python run_all_families_bernoulli.py -n 32 --ordt 15 --ordy 45 --ordp 15 --ordb 15
#     ... --only M3,N3      ... --weights 1,2,3      ... --kernel "MathKernel -script"
# =============================================================================
import argparse
import glob
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")
sys.path.insert(0, CODE)

ap = argparse.ArgumentParser()
ap.add_argument("-n", "--ncores", type=int, default=os.cpu_count())
ap.add_argument("--ordt", type=int, required=True, help="order of the recursion run (names the workdir dist_<ordt>_<ordy>)")
ap.add_argument("--ordy", type=int, required=True)
ap.add_argument("--ordp", type=int, required=True, help="physical order of the phys_ files (as --phys in the recursion)")
ap.add_argument("--ordb", type=int, required=True, help="Bernoulli truncation order")
ap.add_argument("--weights", default="1,2,3,4,5,6")
ap.add_argument("-c", "--corners", default="1,2,3")
ap.add_argument("--only", default=None)
ap.add_argument("--kernel", default="")
ap.add_argument("--redo", action="store_true", help="recompute ber_ files that already exist")
ap.add_argument("--python", default=sys.executable)
a = ap.parse_args()

fams = sorted((d for d in os.listdir(HERE) if os.path.isfile(os.path.join(HERE, d, "inputs", "atilde.m"))),
              key=lambda d: os.path.getsize(os.path.join(HERE, d, "inputs", "atilde.m")))
if a.only:
    fams = [f for f in fams if f in a.only.split(",")]
weights = [int(w) for w in a.weights.split(",")]
corners = [int(c) for c in a.corners.split(",")]
wdn = "dist_%d_%d" % (a.ordt, a.ordy)

report = []
T0 = time.time()
for F in fams:
    wd = os.path.join(HERE, F, wdn)
    if not os.path.isdir(wd):
        report.append((F, "no %s -- recursion not run" % wdn)); continue
    print("\n######## %s  (%s)" % (F, time.strftime("%H:%M:%S")), flush=True)
    status = []
    for w in weights:
        # 1) make sure the phys_ files exist for every corner
        for c in corners:
            ph = os.path.join(wd, "phys_c%d_w%d_ord%d.m" % (c, w, a.ordp))
            if os.path.isfile(ph):
                continue
            if not os.path.isfile(os.path.join(wd, "rec_c%d_w%d.pkl" % (c, w))):
                status.append("c%dw%d: no recursion result" % (c, w)); continue
            print("  writing %s from rec_c%d_w%d.pkl" % (os.path.basename(ph), c, w), flush=True)
            import pipe_common as pc, phys           # needs numpy/gmpy2 -> run with venv python
            C = pc.load_pickle(pc.setup_path(wd, c))['C']
            phys.write_phys(C, phys.slots_to_phys(C, pc.load_pickle(pc.rec_path(wd, c, w)), C.ORDt, a.ordp), ph)
        # 2) Bernoulli for the corners that have phys_ and (unless --redo) no ber_ yet
        cs = [c for c in corners
              if os.path.isfile(os.path.join(wd, "phys_c%d_w%d_ord%d.m" % (c, w, a.ordp)))
              and (a.redo or not os.path.isfile(os.path.join(wd, "ber_c%d_w%d_ordp%d_ordb%d.m" % (c, w, a.ordp, a.ordb))))]
        if not cs:
            continue
        cmd = [a.python, os.path.join(CODE, "run_bernoulli.py"), "-w", str(w), "-n", str(a.ncores),
               "-c", ",".join(map(str, cs)), "--ordp", str(a.ordp), "--ordb", str(a.ordb),
               "--workdir", wd, "--physdir", wd, "--outdir", wd]
        if a.kernel:
            cmd += ["--kernel", a.kernel]
        rc = subprocess.run(cmd).returncode
        if rc:
            status.append("w%d: run_bernoulli rc=%d" % (w, rc))
    have = len(glob.glob(os.path.join(wd, "ber_c*_ordp%d_ordb%d.m" % (a.ordp, a.ordb))))
    report.append((F, "%d/%d ber files%s" % (have, len(weights) * len(corners),
                                               ("  (" + "; ".join(status) + ")") if status else "")))

print("\n================ BERNOULLI SUMMARY (%.1f min) ================" % ((time.time() - T0) / 60))
for F, s in report:
    print("  %-4s %s" % (F, s))
