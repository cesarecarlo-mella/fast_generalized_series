#!/usr/bin/env python3
# =============================================================================
#  run_multi_order.py  --  ONE call, the full two-stage pipeline (run_compute.py
#                           + run_analyze.py) across MULTIPLE orders, auto-
#                           chopping any order below --ceiling that isn't
#                           already sitting in --datadir.
# -----------------------------------------------------------------------------
#  For each order in --orders:
#    1. chop_phys_ber.wl  (skipped if the order == --ceiling, or the files
#       already exist) -- per corner, per weight, cheap truncation only.
#    2. run_compute.py    (stage 1, expensive: raw differences vs reference)
#    3. run_analyze.py    (stage 2, cheap: MIN/MEDIAN with --chop threshold)
#  All orders share ONE grid file (same --outdir, same --npoints), so every
#  order is directly comparable point-for-point.
#
#  USAGE
#     ./run_multi_order.py --datadir <dir with phys_/ber_ at ceiling order> \
#         --ceiling 20 --orders 12,16,20 --weights 1,2 \
#         --modes c1only,c2only,c3only,coverage --npoints 3000 -n 8 \
#         --chop 1e-4 --outdir ../cluster_results/comparison_v2
#
#  Weight 6 (self-convergence): --orders here means the SAME list run_compute.py
#  takes as its ORDp=ORDb pair for weight 1/2, but weight 6 needs a LOW/HIGH
#  PAIR instead (see run_compute.py). Pass --w6-pairs "18,19;19,20" to run
#  those pairs instead of --orders/--weights for weight 6 -- run it separately
#  from weight 1/2 (same restriction as run_compute.py itself).
#
#  Weight 6 ALSO auto-builds (once per corner/pair, skipped if already
#  present) a coefficient cache via 00_build_coeffcache.wl before each
#  run_compute.py call -- this is the dominant cost of weight-6 evaluation
#  (~15-35 min/corner, expression-dependent, NOT grid-dependent), and every
#  parallel point-slicing worker 01_eval_grid.wl spawns would otherwise redo
#  it from scratch. Written straight into --datadir (required location --
#  01_eval_grid.wl looks for it next to phys_/ber_, not in --outdir).
#
#  RESULTS: land in --outdir as compare_<mode>_w<w>_ordp<ORDp>_ordb<ORDb>.m
#  (the final MIN/MEDIAN digit-agreement tables) -- copy these back to a
#  local compare_results/results/ (or wherever plot_compare.wl's resultsDir
#  points) to visualize. rawdata_*.m in the same --outdir are the
#  intermediate per-point differences (needed to rerun 03_analyze.wl with a
#  different --chop without recomputing), not the final numbers themselves.
# =============================================================================

import argparse
import os
import shlex
import shutil
import subprocess
import sys

SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
CORNERS_FOR_MODE = {"c1only": [1], "c2only": [2], "c3only": [3], "coverage": [1, 2, 3]}


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


def run(cmd, label):
    print(f"\n==== {label} ====", flush=True)
    print(f"$ {' '.join(str(c) for c in cmd)}", flush=True)
    rc = subprocess.call([str(c) for c in cmd])
    if rc != 0:
        print(f"ERROR: {label} failed (rc={rc})"); sys.exit(1)


def nonempty(p, minb=2):
    return os.path.isfile(p) and os.path.getsize(p) >= minb


ap = argparse.ArgumentParser(description="multi-order driver for the two-stage compute+analyze pipeline, with auto-chopping of sub-ceiling orders")
ap.add_argument("--datadir", required=True, help="dir with phys_c<c>_w<w>_ord<ceiling>.m / ber_..._ordp<ceiling>_ordb<ceiling>.m")
ap.add_argument("--ceiling", type=int, required=True, help="the order phys_/ber_ actually exist at (real data, no chopping needed)")
ap.add_argument("--orders", default="12,16,20", help="comma list of orders to run (weight 1/2 only -- any < ceiling get chopped first)")
ap.add_argument("--w6-pairs", default="", help="semicolon-separated low,high pairs for weight 6, e.g. '18,19;19,20' -- run instead of --orders/--weights, separately")
ap.add_argument("--weights", default="3")
ap.add_argument("--modes", default="c1only,c2only,c3only,coverage")
ap.add_argument("--npoints", type=int, required=True)
ap.add_argument("-n", "--ncores", type=int, required=True)
ap.add_argument("--chop", type=float, default=1e-4, help="analyze-stage chop threshold (default 1e-4)")
ap.add_argument("--outdir", required=True)
ap.add_argument("--kernel", default=os.environ.get("PIPELINE_KERNEL", ""))
args = ap.parse_args()

KERNEL = resolve_kernel(args.kernel)
DATADIR = os.path.abspath(args.datadir)
OUTDIR = os.path.abspath(args.outdir)
os.makedirs(OUTDIR, exist_ok=True)
modes = [s.strip() for s in args.modes.split(",") if s.strip()]
corners_needed = sorted({c for m in modes for c in CORNERS_FOR_MODE[m]})
kernel_args = ["--kernel", args.kernel] if args.kernel else []
COMPUTE = os.path.join(SCRIPTDIR, "run_compute.py")
ANALYZE = os.path.join(SCRIPTDIR, "run_analyze.py")
CHOP = os.path.join(SCRIPTDIR, "chop_phys_ber.wl")
CACHE = os.path.join(SCRIPTDIR, "00_build_coeffcache.wl")


def chop_if_needed(weight, order):
    if order == args.ceiling:
        return
    for c in corners_needed:
        physF = os.path.join(DATADIR, f"phys_c{c}_w{weight}_ord{order}.m")
        berF = os.path.join(DATADIR, f"ber_c{c}_w{weight}_ordp{order}_ordb{order}.m")
        if nonempty(physF) and nonempty(berF):
            continue
        run(KERNEL + [CHOP, c, weight, args.ceiling, order, DATADIR, DATADIR],
            f"chop weight {weight}, corner {c}, order {order}")


def cache_if_needed_w6(lo, hi):
    # 01_eval_grid.wl's weight-6 fast path looks for coeffcache_c<c>_w6_ordp<lo>_ordb<hi>.m
    # in --datadir specifically (NOT --outdir) -- it must live next to phys_/ber_, matching
    # what every 01_eval_grid.wl worker will pass as its own datadir. Building it here, once,
    # before run_compute.py spawns N parallel workers, means every worker loads the same
    # cache instead of each independently redoing the ~15-35 min/corner Expand+
    # PowerExpand+CoefficientRules step (see 00_build_coeffcache.wl's own header).
    for c in corners_needed:
        cacheF = os.path.join(DATADIR, f"coeffcache_c{c}_w6_ordp{lo}_ordb{hi}.m")
        if nonempty(cacheF):
            print(f"  coeffcache_c{c}_w6_ordp{lo}_ordb{hi}.m exists -- skipping")
            continue
        run(KERNEL + [CACHE, c, lo, hi, DATADIR, DATADIR],
            f"build coeff cache corner {c}, weight 6, N={lo},{hi}")


if args.w6_pairs:
    pairs = [tuple(int(x) for x in p.split(",")) for p in args.w6_pairs.split(";") if p.strip()]
    print(f"=== MULTI-ORDER (weight 6 self-convergence) ===  pairs {pairs}  modes {modes}  ncores {args.ncores}")
    for lo, hi in pairs:
        for o in (lo, hi):
            chop_if_needed(6, o)
    for lo, hi in pairs:
        print(f"\n---------------- weight 6: N={lo} vs N={hi}  (building coeff cache first) ----------------")
        cache_if_needed_w6(lo, hi)
        print(f"\n---------------- weight 6: N={lo} vs N={hi} ----------------")
        run([sys.executable, COMPUTE, "--datadir", DATADIR, "--weights", "6", "--modes", ",".join(modes),
             "--npoints", args.npoints, "--ordp", lo, "--ordb", hi, "-n", args.ncores,
             "--outdir", OUTDIR] + kernel_args, f"compute weight 6 N={lo},{hi}")
        run([sys.executable, ANALYZE, "--datadir", OUTDIR, "--weights", "6", "--modes", ",".join(modes),
             "--ordp", lo, "--ordb", hi, "--chop", args.chop, "--outdir", OUTDIR] + kernel_args,
            f"analyze weight 6 N={lo},{hi}")
else:
    orders = sorted({int(s) for s in args.orders.split(",") if s.strip()})
    weights = [int(s) for s in args.weights.split(",") if s.strip()]
    print(f"=== MULTI-ORDER ===  orders {orders}  weights {weights}  modes {modes}  ncores {args.ncores}  chop {args.chop}")
    for O in orders:
        if O > args.ceiling:
            print(f"ERROR: order {O} > ceiling {args.ceiling}"); sys.exit(1)
        for w in weights:
            chop_if_needed(w, O)
    for O in orders:
        print(f"\n---------------- order {O} ----------------")
        run([sys.executable, COMPUTE, "--datadir", DATADIR, "--weights", args.weights, "--modes", ",".join(modes),
             "--npoints", args.npoints, "--ordp", O, "--ordb", O, "-n", args.ncores,
             "--outdir", OUTDIR] + kernel_args, f"compute order {O}")
        run([sys.executable, ANALYZE, "--datadir", OUTDIR, "--weights", args.weights, "--modes", ",".join(modes),
             "--ordp", O, "--ordb", O, "--chop", args.chop, "--outdir", OUTDIR] + kernel_args,
            f"analyze order {O}")

print("\n=== MULTI-ORDER DONE ===")
print(f"results in {OUTDIR}")
