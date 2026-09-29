#!/usr/bin/env python3
# =============================================================================
#  run_bernoulli.py  --  BERNOULLI (CONFORMAL) RESUMMATION  (steps 5-6),
#                         independent of run_weight.py / run_to_physical.py --
#                         the blow_up_form(_2) equivalent of the old
#                         pipeline's run_emit_bernoulli.py.
# -----------------------------------------------------------------------------
#  Turns  phys_c<c>_w<w>_ord<ORDp>.m  (produced by run_to_physical.py, a
#  series in the physical variables x,y,z) into
#  ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>.m , the same series re-expanded in the
#  conformal variables:
#      x -> 1-Exp[-sx]   y -> 1-Exp[-sy]   z -> 1-Exp[-sz]
#  (only the corner's own two physical variables get transformed: c1: x,y;
#  c2: z,y; c3: x,z -- y is always y->sy, no corner-specific handling),
#  truncated to order ORDb in each new s-variable. Log[x],Log[y],Log[z] are
#  left as LogX,LogY,LogZ placeholders in the output (not resolved) -- see
#  05_bernoulli.wl header for the full rationale.
#  Completely separate from the recursion / physical-variable conversion --
#  run it whenever, as many times as you like, at whatever ORDb.
#
#       05_bernoulli.wl        c w p N ORDp ORDb   N procs  (slice by component)
#       06_merge_bernoulli.wl  c w ORDp ORDb N outdir   1  proc  stitch the slices
#
#  reads   :  <physdir>/phys_c<c>_w<w>_ord<ORDp>.m   (physdir default <workdir>/out)
#  writes  :  <outdir>/ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>.m  (outdir default physdir)
#
#  Works UNCHANGED against either blow_up_form/pipeline/ or
#  blow_up_form_2/pipeline/ -- copied as-is into both.
#
#  USAGE
#     ./run_bernoulli.py -w 2 -n 6 --ordp 15 --ordb 15
#     ./run_bernoulli.py -w 2 -n 6 --ordp 15 --ordb 15 -c 1,2 --workdir dist \
#                         --physdir dist/out --outdir dist/out \
#                         --kernel "MathKernel -script"
# =============================================================================

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time

SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))


def wl_path(name):
    return os.path.join(SCRIPTDIR, name)


def resolve_kernel(requested):
    for cand in ([requested] if requested else []) + \
                ["wolframscript -file", "MathKernel -script", "math -script",
                 "WolframKernel -script", "wolfram -script"]:
        if not cand:
            continue
        toks = shlex.split(cand)
        if shutil.which(toks[0]):
            return [shutil.which(toks[0])] + toks[1:]
    print('ERROR: no Mathematica kernel found. Pass --kernel "MathKernel -script"')
    sys.exit(1)


def spawn(kernel, script, cmdargs, workdir, logpath):
    fh = open(logpath, "w")
    cmd = kernel + [wl_path(script)] + [str(a) for a in cmdargs]
    return subprocess.Popen(cmd, cwd=workdir, stdout=fh, stderr=subprocess.STDOUT), fh


def run_streaming(kernel, script, cmdargs, workdir, logpath):
    cmd = kernel + [wl_path(script)] + [str(a) for a in cmdargs]
    print("    $ " + " ".join(cmd), flush=True)
    with open(logpath, "w") as fh:
        pr = subprocess.Popen(cmd, cwd=workdir, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in pr.stdout:
            print("    " + line.rstrip(), flush=True)
            fh.write(line)
        pr.wait()
    return pr.returncode


def nonempty(p, minb=2):
    return os.path.isfile(p) and os.path.getsize(p) >= minb


# ---------------------------------------------------------------- args
ap = argparse.ArgumentParser(description="Bernoulli/conformal resummation (steps 5-6), standalone")
ap.add_argument("-w", "--weight", type=int, required=True)
ap.add_argument("-n", "--ncores", type=int, required=True, help="parallel processes = slices")
ap.add_argument("-c", "--corners", default="1,2,3")
ap.add_argument("--ordp", type=int, required=True, help="order of the input phys_..._ord<ORDp>.m to read")
ap.add_argument("--ordb", type=int, required=True, help="truncation order in each new s-variable")
ap.add_argument("--workdir", default=".")
ap.add_argument("--physdir", default="", help="dir holding phys_ files (default <workdir>/out)")
ap.add_argument("--outdir", default="", help="dir for ber_ files (default same as --physdir)")
ap.add_argument("--kernel", default=os.environ.get("PIPELINE_KERNEL", ""))
ap.add_argument("--keep-scratch", action="store_true")
args = ap.parse_args()

KERNEL = resolve_kernel(args.kernel)
WORK = os.path.abspath(args.workdir)
PHYSDIR = os.path.abspath(args.physdir) if args.physdir else os.path.join(WORK, "out")
OUTDIR = os.path.abspath(args.outdir) if args.outdir else PHYSDIR
os.makedirs(OUTDIR, exist_ok=True)
weight, ncores, ORDp, ORDb = args.weight, args.ncores, args.ordp, args.ordb
corners = [int(s) for s in args.corners.split(",") if s.strip()]


def pp(name):
    return os.path.join(PHYSDIR, name)


print(f"\n=== BERNOULLI  weight {weight} | corners {corners} | N = {ncores} | ORDp {ORDp} ORDb {ORDb} ===")
print(f"    physdir {PHYSDIR}")
print(f"    outdir  {OUTDIR}")
print(f"    kernel  {' '.join(KERNEL)}\n")

missing = [f"phys_c{c}_w{weight}_ord{ORDp}.m" for c in corners
           if not nonempty(pp(f"phys_c{c}_w{weight}_ord{ORDp}.m"))]
if missing:
    print("REFUSING TO START -- missing in physdir:")
    for m in missing:
        print("   ", m)
    sys.exit(1)

overall_ok = True
for c in corners:
    print(f"\n----------------------------------------------------------------")
    print(f"  CORNER {c}   weight {weight}   ORDp {ORDp}  ORDb {ORDb}   ({time.strftime('%H:%M:%S')})")
    print(f"----------------------------------------------------------------")
    t_c = time.time()
    ber_f = os.path.join(OUTDIR, f"ber_c{c}_w{weight}_ordp{ORDp}_ordb{ORDb}.m")
    if nonempty(ber_f):
        print(f"  ber_c{c}_w{weight}_ordp{ORDp}_ordb{ORDb}.m exists -- skipping"); continue

    # -- step 5 : N parallel bernoulli processes, live progress, retry x3 ---
    print(f"  [1/2] 05_bernoulli  x{ncores}")

    def prog(p):
        try:
            return int(open(pp(f"prog05_c{c}_w{weight}_p{p}.txt")).read().replace('"', " ").split()[0])
        except Exception:
            return 0

    def slice_ok(p):
        return nonempty(pp(f"ber_c{c}_w{weight}_ordp{ORDp}_ordb{ORDb}_p{p}_{ncores}.m"))

    pending = list(range(ncores))
    for attempt in range(1, 4):
        if not pending:
            break
        if attempt > 1:
            print(f"     retry {attempt} for slices {pending}")
        procs = {}
        for p in pending:
            for f in (pp(f"ber_c{c}_w{weight}_ordp{ORDp}_ordb{ORDb}_p{p}_{ncores}.m"),
                      pp(f"prog05_c{c}_w{weight}_p{p}.txt")):
                if os.path.isfile(f):
                    os.remove(f)
            procs[p] = spawn(KERNEL, "05_bernoulli.wl", [c, weight, p, ncores, ORDp, ORDb],
                             PHYSDIR, pp(f"log_05_c{c}_w{weight}_p{p}.txt"))
        t0 = time.time()
        while any(pr.poll() is None for pr, _ in procs.values()):
            time.sleep(5)
            done = sum(prog(p) for p in procs)
            alive = sum(1 for pr, _ in procs.values() if pr.poll() is None)
            print(f"    {done} components done  {alive}/{ncores} procs  "
                  f"{time.time()-t0:5.0f}s", flush=True)
        for pr, fh in procs.values():
            pr.wait(); fh.close()
        pending = [p for p in pending if not slice_ok(p)]
    if pending:
        print(f"  !! bernoulli slices {pending} failed after 3 attempts")
        overall_ok = False; continue

    # -- step 6 : merge ---------------------------------------------------
    print(f"\n  [2/2] 06_merge_bernoulli")
    rc = run_streaming(KERNEL, "06_merge_bernoulli.wl", [c, weight, ORDp, ORDb, ncores, OUTDIR], PHYSDIR,
                       pp(f"log_06_c{c}_w{weight}_ordp{ORDp}_ordb{ORDb}.txt"))
    if rc != 0 or not nonempty(ber_f):
        print(f"  !! merge failed (rc={rc})"); overall_ok = False; continue

    if not args.keep_scratch:
        for p in range(ncores):
            for nm in (f"ber_c{c}_w{weight}_ordp{ORDp}_ordb{ORDb}_p{p}_{ncores}.m",
                       f"prog05_c{c}_w{weight}_p{p}.txt"):
                if os.path.isfile(pp(nm)):
                    os.remove(pp(nm))

    print(f"\n  corner {c} weight {weight} ordp {ORDp} ordb {ORDb} DONE  ({(time.time()-t_c)/60:.1f} min)"
          f"  ->  {os.path.relpath(ber_f)}")

print("\n================================================================")
print(f"  bernoulli weight {weight} ordp {ORDp} ordb {ORDb} : "
      f"{'ALL CORNERS DONE' if overall_ok else 'INCOMPLETE -- see logs'}")
print("================================================================")
sys.exit(0 if overall_ok else 1)
