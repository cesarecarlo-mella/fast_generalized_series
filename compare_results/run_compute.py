#!/usr/bin/env python3
# =============================================================================
#  run_compute.py  --  STAGE 1 (expensive): raw per-component differences,
#                       parallel driver
# -----------------------------------------------------------------------------
#  Evaluates the blow-up pipeline's output (phys_c<c>_w<w>_ord<ORDp>.m native +
#  ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>.m Bernoulli/conformal) against the
#  reference (exact solution for weight 1/2, higher-order truncation for
#  weight 6) over a shared grid, and saves the RAW signed differences + the
#  reference vector -- no digit-count, MIN/MEDIAN, or chopping decision is
#  made here. That's run_analyze.py's job, deliberately kept separate and
#  cheap so tuning the analysis never requires re-running this step.
#
#       gen_grid.wl        N gridfile                1  proc  (once, cached)
#       01_eval_grid.wl    mode w p N gridfile datadir ORDp ORDb   N procs
#       02_merge_raw.wl    mode w ORDp ORDb N outdir              1  proc
#
#  reads   :  <datadir>/phys_c<c>_w<w>_ord<ORDp>.m ,
#             <datadir>/ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>.m
#  writes  :  <outdir>/rawdata_<mode>_w<w>_ordp<ORDp>_ordb<ORDb>.m
#
#  USAGE
#     ./run_compute.py --datadir /path/to/blow_up_form_2_order_20/out \
#                       --weights 1,2 --modes c1only,coverage \
#                       --npoints 1500 -n 8
#     ./run_compute.py --datadir ... --weights 6 --ordp 18 --ordb 19 -n 1  (self-conv, see MEMORY WARNING below)
#
#  MEMORY WARNING (weight 6 specifically): every one of the -n parallel
#  processes independently loads its OWN full copy of the needed phys_/ber_
#  data before touching a single point -- fine for small weight-2 series,
#  dangerous for weight-6 series (100s of MB on disk, much bigger once
#  parsed). Keep -n small for weight 6; see 01_eval_grid.wl / the incident
#  notes in this project's history before raising it.
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


CORNERS_NEEDED = {"c1only": [1], "c2only": [2], "c3only": [3], "coverage": [1, 2, 3]}

# ---------------------------------------------------------------- args
ap = argparse.ArgumentParser(description="STAGE 1: raw per-component differences vs reference, saved to disk (no digit-count decisions made here)")
ap.add_argument("--weights", default="1,2", help="comma list, only 1 and/or 2 have an exact reference")
ap.add_argument("--modes", default="c1only,c2only,c3only,coverage",
                 help="comma list: c1only, c2only, c3only, coverage")
ap.add_argument("-n", "--ncores", type=int, required=True, help="parallel processes = slices")
ap.add_argument("--npoints", type=int, default=1500, help="grid points over the triangle (default 1500)")
ap.add_argument("--datadir", required=True, help="dir with phys_c<c>_w<w>_ord<ORDp>.m / ber_..._ordp<ORDp>_ordb<ORDb>.m")
ap.add_argument("--ordp", type=int, default=15, help="which phys_ truncation order to read (default 15); weight 6: LOWER order")
ap.add_argument("--ordb", type=int, default=15, help="which ber_ truncation order to read (default 15); weight 6: HIGHER order (reference)")
ap.add_argument("--griddir", default="", help="where the shared grid file lives (default --outdir)")
ap.add_argument("--outdir", default="", help="dir for rawdata_ results (default: this script's own dir, i.e. compare_results/results)")
ap.add_argument("--kernel", default=os.environ.get("PIPELINE_KERNEL", ""))
ap.add_argument("--keep-scratch", action="store_true")
args = ap.parse_args()

KERNEL = resolve_kernel(args.kernel)
DATADIR = os.path.abspath(args.datadir)
OUTDIR = os.path.abspath(args.outdir) if args.outdir else os.path.join(SCRIPTDIR, "results")
GRIDDIR = os.path.abspath(args.griddir) if args.griddir else OUTDIR
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(GRIDDIR, exist_ok=True)
ncores, npoints, ORDp, ORDb = args.ncores, args.npoints, args.ordp, args.ordb
weights = [int(s) for s in args.weights.split(",") if s.strip()]
modes = [s.strip() for s in args.modes.split(",") if s.strip()]
for m in modes:
    if m not in CORNERS_NEEDED:
        print(f"ERROR: unknown mode {m!r} (must be c1only, c2only, c3only, or coverage)"); sys.exit(1)
for w in weights:
    if w not in (1, 2, 6):
        print(f"ERROR: weight {w} not supported -- 1,2 (exact reference) or 6 (self-convergence)"); sys.exit(1)
if any(w == 6 for w in weights) and len(weights) > 1:
    print("ERROR: weight 6 (self-convergence) can't be mixed with 1/2 (exact) in one call -- "
          "ORDp/ORDb mean different things (low/high order vs native/Bernoulli order). Run separately.")
    sys.exit(1)


def dp(name):
    return os.path.join(DATADIR, name)


def op(name):
    return os.path.join(OUTDIR, name)


print(f"\n=== COMPUTE (stage 1/2)  weights {weights} | modes {modes} | N = {ncores} | npoints {npoints} "
      f"| ORDp {ORDp} ORDb {ORDb} ===")
print(f"    datadir {DATADIR}")
print(f"    outdir  {OUTDIR}")
print(f"    kernel  {' '.join(KERNEL)}\n")

# -- step 0 : shared grid (once, cached by point count) ---------------------
gridfile = os.path.join(GRIDDIR, f"grid_N{npoints}.m")
if not nonempty(gridfile):
    print(f"[0/2] generating shared grid ({npoints} points) -> {gridfile}")
    rc = run_streaming(KERNEL, "gen_grid.wl", [npoints, gridfile], GRIDDIR,
                       os.path.join(GRIDDIR, f"log_grid_N{npoints}.txt"))
    if rc != 0 or not nonempty(gridfile):
        print("ERROR: gen_grid.wl failed"); sys.exit(1)
else:
    print(f"[0/2] shared grid already exists: {gridfile} -- reusing")

overall_ok = True
for mode in modes:
    corners = CORNERS_NEEDED[mode]
    for w in weights:
        print(f"\n----------------------------------------------------------------")
        print(f"  MODE {mode}   weight {w}   ORDp {ORDp}  ORDb {ORDb}   ({time.strftime('%H:%M:%S')})")
        print(f"----------------------------------------------------------------")
        t_c = time.time()
        out_f = op(f"rawdata_{mode}_w{w}_ordp{ORDp}_ordb{ORDb}.m")
        if nonempty(out_f):
            print(f"  rawdata_{mode}_w{w}_ordp{ORDp}_ordb{ORDb}.m exists -- skipping"); continue

        missing = []
        if w == 6:
            needed = [f"phys_c{c}_w{w}_ord{o}.m" for c in corners for o in (ORDp, ORDb)] + \
                     [f"ber_c{c}_w{w}_ordp{o}_ordb{o}.m" for c in corners for o in (ORDp, ORDb)]
        else:
            needed = [f"phys_c{c}_w{w}_ord{ORDp}.m" for c in corners] + \
                     [f"ber_c{c}_w{w}_ordp{ORDp}_ordb{ORDb}.m" for c in corners]
        for nm in needed:
            if not nonempty(dp(nm)):
                missing.append(nm)
        if missing:
            print(f"  REFUSING -- missing in datadir:")
            for m in missing:
                print("   ", m)
            overall_ok = False; continue

        # -- step 1 : N parallel eval processes, live progress, retry x3
        print(f"  [1/2] 01_eval_grid  x{ncores}")

        def prog(p):
            try:
                return int(open(op(f"progEVAL_{mode}_w{w}_p{p}.txt")).read().replace('"', " ").split()[0])
            except Exception:
                return 0

        def slice_ok(p):
            return nonempty(op(f"raw_{mode}_w{w}_p{p}_{ncores}.m"))

        pending = list(range(ncores))
        for attempt in range(1, 4):
            if not pending:
                break
            if attempt > 1:
                print(f"     retry {attempt} for slices {pending}")
            procs = {}
            for p in pending:
                for f in (op(f"raw_{mode}_w{w}_p{p}_{ncores}.m"), op(f"progEVAL_{mode}_w{w}_p{p}.txt")):
                    if os.path.isfile(f):
                        os.remove(f)
                procs[p] = spawn(KERNEL, "01_eval_grid.wl",
                                 [mode, w, p, ncores, gridfile, DATADIR, ORDp, ORDb],
                                 OUTDIR, op(f"log_01_{mode}_w{w}_p{p}.txt"))
            t0 = time.time()
            while any(pr.poll() is None for pr, _ in procs.values()):
                time.sleep(5)
                done = sum(prog(p) for p in procs)
                alive = sum(1 for pr, _ in procs.values() if pr.poll() is None)
                print(f"    {done} points done  {alive}/{ncores} procs  "
                      f"{time.time()-t0:5.0f}s", flush=True)
            for pr, fh in procs.values():
                pr.wait(); fh.close()
            pending = [p for p in pending if not slice_ok(p)]
        if pending:
            print(f"  !! eval slices {pending} failed after 3 attempts")
            overall_ok = False; continue

        # -- step 2 : merge ----------------------------------------------
        print(f"\n  [2/2] 02_merge_raw")
        rc = run_streaming(KERNEL, "02_merge_raw.wl", [mode, w, ORDp, ORDb, ncores, OUTDIR], OUTDIR,
                           op(f"log_02_{mode}_w{w}.txt"))
        if rc != 0 or not nonempty(out_f):
            print(f"  !! merge failed (rc={rc})"); overall_ok = False; continue

        if not args.keep_scratch:
            for p in range(ncores):
                for nm in (f"raw_{mode}_w{w}_p{p}_{ncores}.m", f"progEVAL_{mode}_w{w}_p{p}.txt"):
                    if os.path.isfile(op(nm)):
                        os.remove(op(nm))

        print(f"\n  {mode} weight {w} DONE  ({(time.time()-t_c)/60:.1f} min)  ->  {os.path.relpath(out_f)}")

print("\n================================================================")
print(f"  compute : {'ALL DONE' if overall_ok else 'INCOMPLETE -- see logs'}")
print("================================================================")
sys.exit(0 if overall_ok else 1)
