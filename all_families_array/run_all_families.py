#!/usr/bin/env python3
# =============================================================================
#  run_all_families.py -- run the array/modular blow-up recursion (weights
#  1..W, all corners) for EVERY family in all_families_array/, one after the
#  other, cheapest first (ordered by the size of <Family>/inputs/atilde.m).
#
#  Layout (fully self-contained -- copy the whole folder to the cluster):
#     code/                 the pipeline (run_array.py, 00/01/02, engine, ...)
#     letters/              master letters, shared by all families (15/45, 20 letters)
#     <Family>/inputs/      atilde.m, sol0.m, bdy_c<c>_w<w>.m
#     <Family>/dist_<T>_<Y>/  results: exact_c<c>_w<w>.m, phys_c<c>_w<w>_ord<ORD>.m
#
#  Families are independent: a failure is reported and the next family runs.
#  Everything is resumable (finished primes / weights / families are skipped).
#
#  USAGE
#     ./run_all_families.py -n 32 --ordt 15 --ordy 45 -w 6 --phys 15
#     ./run_all_families.py -n 32 --ordt 15 --ordy 45 -w 6 --only M3,N3,A
#     ./run_all_families.py -n 32 --ordt 15 --ordy 45 -w 6 --resume-from E1
#     ./run_all_families.py --dry-run --ordt 15 --ordy 45 -w 6
#  Any extra option after "--" is passed to run_array.py unchanged, e.g.
#     ./run_all_families.py -n 32 --ordt 15 --ordy 45 -w 6 -- --nprimes 40
# =============================================================================
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")
LETTERS = os.path.join(HERE, "letters")

argv = sys.argv[1:]
extra = []
if "--" in argv:
    i = argv.index("--"); extra = argv[i + 1:]; argv = argv[:i]
ap = argparse.ArgumentParser()
ap.add_argument("-n", "--ncores", type=int, default=os.cpu_count())
ap.add_argument("-w", "--weight", type=int, default=6)
ap.add_argument("-c", "--corners", default="1,2,3")
ap.add_argument("--ordt", type=int, required=True)
ap.add_argument("--ordy", type=int, required=True)
ap.add_argument("--phys", type=int, default=None, help="also write phys_..._ord<ORD>.m (e.g. = ordt)")
ap.add_argument("--only", default=None, help="comma list of families")
ap.add_argument("--resume-from", default=None, help="skip families before this one (in run order)")
ap.add_argument("--python", default=sys.executable)
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args(argv)

# fail fast (instead of once per family) if this interpreter lacks the packages
_chk = subprocess.run([a.python, "-c", "import numpy, gmpy2"], capture_output=True, text=True)
if _chk.returncode != 0:
    sys.exit("ERROR: %s cannot import numpy/gmpy2:\n  %s\n"
             "Run with the venv's interpreter, e.g.  venv/bin/python run_all_families.py ..."
             % (a.python, _chk.stderr.strip().splitlines()[-1] if _chk.stderr else "?"))

fams = sorted((d for d in os.listdir(HERE)
               if os.path.isfile(os.path.join(HERE, d, "inputs", "atilde.m"))),
              key=lambda d: os.path.getsize(os.path.join(HERE, d, "inputs", "atilde.m")))
if a.only:
    want = [f.strip() for f in a.only.split(",") if f.strip()]
    missing = [f for f in want if f not in fams]
    if missing:
        sys.exit("unknown families: %s  (have: %s)" % (missing, fams))
    fams = [f for f in fams if f in want]
if a.resume_from:
    if a.resume_from not in fams:
        sys.exit("--resume-from %s: not in %s" % (a.resume_from, fams))
    fams = fams[fams.index(a.resume_from):]

wd_name = "dist_%d_%d" % (a.ordt, a.ordy)
print("=== ALL FAMILIES  weights 1..%d | corners %s | order %d/%d | N = %d ===" %
      (a.weight, a.corners, a.ordt, a.ordy, a.ncores))
for f in fams:
    print("   %-4s  atilde %7.1f kB   -> %s/%s" %
          (f, os.path.getsize(os.path.join(HERE, f, "inputs", "atilde.m")) / 1e3, f, wd_name))
if a.dry_run:
    sys.exit(0)

ok, bad = [], []
T0 = time.time()
for f in fams:
    t0 = time.time()
    cmd = [a.python, os.path.join(CODE, "run_array.py"), "-w", str(a.weight), "-n", str(a.ncores),
           "-c", a.corners, "--ordt", str(a.ordt), "--ordy", str(a.ordy),
           "--workdir", os.path.join(HERE, f, wd_name),
           "--inputs", os.path.join(HERE, f, "inputs"), "--letters", LETTERS]
    if a.phys is not None:
        cmd += ["--phys", str(a.phys)]
    cmd += extra
    print("\n################  FAMILY %s  (%s)  ################" % (f, time.strftime("%H:%M:%S")), flush=True)
    rc = subprocess.run(cmd).returncode
    (ok if rc == 0 else bad).append(f)
    print("################  %s %s  (%.1f min)" % (f, "OK" if rc == 0 else "FAILED rc=%d" % rc,
                                                 (time.time() - t0) / 60), flush=True)

print("\n================================================================")
print("  done in %.1f min   OK: %s" % ((time.time() - T0) / 60, " ".join(ok) or "-"))
print("  FAILED: %s" % (" ".join(bad) or "-"))
print("================================================================")
sys.exit(0 if not bad else 1)
