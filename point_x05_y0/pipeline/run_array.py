#!/usr/bin/env python3
# =============================================================================
#  run_array.py -- driver for the array/modular recursion on ONE node, here for
#                  the expansion around (x, y) = (1/2, 0), "corner" 0:
#                      x = 1/2 + t ,  y = y      (no blow-up needed)
#                  (copy of blow_up_array/pipeline/run_array.py)
# -----------------------------------------------------------------------------
#    00_prepare.py      parse inputs once per corner         -> setup_c<c>.pkl
#    01_chain.py  x P   weights 1..W modulo prime k, -n in parallel
#                                                            -> residues/J_c<c>_w<w>_k<k>.pkl
#    02_reconstruct.py  exact rationals per weight            -> exact_c<c>_w<w>.m
#                       (series in t = x - 1/2 and y, with Log[y]; no change of
#                        variables is needed, so there is no phys_ step)
#  If a weight needs more primes than were run, more 01 jobs are launched
#  automatically (--addprimes at a time) and 02 is retried.
#
#  Unlike run_weight.py, -w W computes ALL weights up to W in one go (each prime
#  job runs the whole chain; nothing here needs the previous weight's exact
#  file). Everything is resumable: rerunning skips finished primes/weights, and
#  a later "-w 7" reuses the weight-6 residues.
#
#  USAGE
#     ./run_array.py -w 6 -n 32 --ordt 20 --ordy 20 --workdir dist_20_20
#         (inputs default to ./inputs: atilde, sol0, bdy_c0_w*, letters 20/20)
#     ./run_array.py -w 6 -n 32 --workdir dist_20_20             (order from .order)
#
#  REQUIRES: python3 >= 3.8, numpy, gmpy2   (no Mathematica, no FORM)
# =============================================================================
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

ap = argparse.ArgumentParser(description="array/modular recursion driver, point (1/2, 0)")
ap.add_argument("-w", "--weight", type=int, required=True, help="compute all weights 1..W")
ap.add_argument("-n", "--ncores", type=int, required=True, help="parallel prime jobs (and reconstruction workers)")
ap.add_argument("-c", "--corners", default="0")
ap.add_argument("--ordt", type=int)
ap.add_argument("--ordy", type=int)
ap.add_argument("--workdir", default=".")
ap.add_argument("--inputs", default=None, help="dir with atilde.m, sol0.m, bdy_c<c>_w<w>.m (default: ./inputs next to this script)")
ap.add_argument("--letters", default=None, help="dir with dlog_d{t,y}_c<c>.m at order >= requested (default: inputs)")
ap.add_argument("--nprimes", type=int, default=None, help="initial number of primes (default 4W+6, enough up to w6 at order 10)")
ap.add_argument("--addprimes", type=int, default=None, help="primes added per extra round (default max(4, ncores))")
ap.add_argument("--allow-missing-bdy", action="store_true",
                help="use a zero boundary for missing bdy_c<c>_w<w>.m (tests only)")
ap.add_argument("--python", default=sys.executable)
args = ap.parse_args()

WORK = os.path.abspath(args.workdir)
os.makedirs(WORK, exist_ok=True)
INPUTS = os.path.abspath(args.inputs or os.path.join(HERE, "inputs"))
LETTERS = os.path.abspath(args.letters or INPUTS)
W, NC = args.weight, args.ncores
corners = [int(s) for s in args.corners.split(",") if s.strip()]
nprimes = args.nprimes or (4 * W + 6)
addp = args.addprimes or max(4, NC)

orderFile = os.path.join(WORK, ".order")
if args.ordt is not None and args.ordy is not None:
    ORDt, ORDy = args.ordt, args.ordy
    if os.path.isfile(orderFile):
        old = tuple(int(x) for x in open(orderFile).read().split())
        if old != (ORDt, ORDy):
            sys.exit("ERROR: workdir %s holds a run at order %s -- use a fresh workdir" % (WORK, old))
    open(orderFile, "w").write("%d %d\n" % (ORDt, ORDy))
else:
    try:
        ORDt, ORDy = (int(x) for x in open(orderFile).read().split())
        print("    (order %d,%d from %s)" % (ORDt, ORDy, orderFile))
    except Exception:
        sys.exit("ERROR: no --ordt/--ordy given and no readable .order file in workdir")

import pipe_common as pc


def run(cmd, log):
    with open(log, "a") as fh:
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode


def logp(name):
    d = os.path.join(WORK, "logs")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


print("\n=== ARRAY RECURSION (1/2, 0)  weights 1..%d | corners %s | N = %d | ORDt %d ORDy %d ===" %
      (W, corners, NC, ORDt, ORDy))
print("    workdir %s\n    inputs  %s\n    letters %s\n" % (WORK, INPUTS, LETTERS))

overall_ok = True
for c in corners:
    tc = time.time()
    print("----------------------------------------------------------------")
    print("  CORNER %d   (%s)" % (c, time.strftime('%H:%M:%S')))
    print("----------------------------------------------------------------")
    if all(os.path.isfile(pc.exact_path(WORK, c, w)) and os.path.isfile(pc.rec_path(WORK, c, w))
           for w in range(1, W + 1)):
        print("  all weights 1..%d already reconstructed -- skipping" % W); continue

    # -- step 0 --------------------------------------------------------------
    sp = pc.setup_path(WORK, c)
    need = True
    if os.path.isfile(sp):
        S = pc.load_pickle(sp)
        need = (S['ordt'], S['ordy']) != (ORDt, ORDy) or S['wmax'] < W
        del S
    if need:
        print("  [0] 00_prepare")
        rc = run([args.python, os.path.join(HERE, "00_prepare.py"), "--corner", str(c), "--ordt", str(ORDt),
                  "--ordy", str(ORDy), "--wmax", str(W), "--inputs", INPUTS, "--letters", LETTERS,
                  "--workdir", WORK] + (["--allow-missing-bdy"] if args.allow_missing_bdy else []),
                 logp("log_00_c%d.txt" % c))
        if rc:
            print("  !! 00_prepare failed -- see logs/log_00_c%d.txt" % c); overall_ok = False; continue

    # -- steps 1+2 -----------------------------------------------------------
    def chain(k):
        for attempt in range(3):
            rc = run([args.python, os.path.join(HERE, "01_chain.py"), "--workdir", WORK, "--corner", str(c),
                      "--k", str(k), "--wmax", str(W)], logp("log_01_c%d_k%03d.txt" % (c, k)))
            if rc == 0 and os.path.isfile(pc.res_path(WORK, c, W, k)):
                return k, True
        return k, False

    todo_w = [w for w in range(1, W + 1)
              if not (os.path.isfile(pc.exact_path(WORK, c, w)) and os.path.isfile(pc.rec_path(WORK, c, w)))]
    target = nprimes
    ok_c = True
    while todo_w:
        have = set(pc.primes_with_residues(WORK, c, W))
        ks = [k for k in range(target) if k not in have]
        if ks:
            print("  [1] 01_chain  %d primes (k = %d..%d)  on %d cores" % (len(ks), ks[0], ks[-1], NC), flush=True)
            t1 = time.time(); done = 0
            with ThreadPoolExecutor(NC) as ex:
                futs = [ex.submit(chain, k) for k in ks]
                last = 0.0
                while True:
                    ndone = sum(f.done() for f in futs)
                    if time.time() - last > 60 or ndone == len(ks):
                        # progress = how many of these primes have reached each weight
                        per_w = [len(set(pc.primes_with_residues(WORK, c, w)) & set(ks)) for w in range(1, W + 1)]
                        print("     %6.0fs  primes finished %d/%d | reached weight: %s" %
                              (time.time() - t1, ndone, len(ks),
                               "  ".join("w%d:%d" % (w, n) for w, n in zip(range(1, W + 1), per_w))), flush=True)
                        last = time.time()
                    if ndone == len(ks):
                        break
                    time.sleep(5)
                for f in futs:
                    k, ok = f.result()
                    if not ok:
                        print("     !! prime k=%d failed 3x -- see logs/log_01_c%d_k%03d.txt" % (k, c, k))
        for w in list(todo_w):
            cmd = [args.python, os.path.join(HERE, "02_reconstruct.py"), "--workdir", WORK, "--corner", str(c),
                   "--weight", str(w), "--jobs", str(NC)]
            rc = run(cmd, logp("log_02_c%d_w%d.txt" % (c, w)))
            lines = [l for l in open(logp("log_02_c%d_w%d.txt" % (c, w))).read().splitlines() if "coefficients" in l]
            print("  [2] " + (lines[-1].strip() if lines else "02_reconstruct w%d rc=%d" % (w, rc)), flush=True)
            if rc == 0:
                todo_w.remove(w)
            elif rc == 3:
                break                  # higher weights need at least as many primes
            else:
                print("  !! 02_reconstruct failed (rc=%d) -- see logs/log_02_c%d_w%d.txt" % (rc, c, w))
                ok_c = False; todo_w = []; break
        if todo_w:
            if target >= len(pc.PRIME_TABLE) - addp:
                print("  !! prime table exhausted"); ok_c = False; break
            target += addp
            print("  -> weights %s need more primes: going to %d" % (todo_w, target), flush=True)
    overall_ok &= ok_c
    print("\n  corner %d weights 1..%d %s  (%.1f min)\n" % (c, W, "DONE" if ok_c else "INCOMPLETE", (time.time() - tc) / 60))

print("================================================================")
print("  weights 1..%d : %s" % (W, "ALL CORNERS DONE" if overall_ok else "INCOMPLETE -- see logs/"))
print("================================================================")
sys.exit(0 if overall_ok else 1)
