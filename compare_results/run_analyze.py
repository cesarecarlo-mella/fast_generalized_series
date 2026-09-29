#!/usr/bin/env python3
# =============================================================================
#  run_analyze.py  --  STAGE 2 (cheap): turn rawdata_*.m (from run_compute.py)
#                       into compare_*.m (digit-count deliverable, MIN + MEDIAN)
# -----------------------------------------------------------------------------
#  Single pass over already-computed raw differences -- no symbolic evaluation,
#  just numeric aggregation -- so this is fast and safe to re-run as many
#  times as you like with a different --chop threshold, entirely independent
#  of run_compute.py. That's the whole point of the two-stage split: tune the
#  analysis without ever re-running the expensive grid evaluation again.
#
#  For each point, per variable (native/Bernoulli): a component is EXCLUDED
#  ("chopped") from that point's MIN/MEDIAN if |ref[component]| < --chop --
#  i.e. that component's reference value is close enough to zero there that
#  a relative-error digit count is not a meaningful measure of convergence
#  (see 03_analyze.wl header for the concrete confirmed example).
#
#  USAGE
#     ./run_analyze.py --datadir results --weights 1,2 --modes c1only,coverage \
#                       --ordp 20 --ordb 20 --chop 1e-4 --outdir results
# =============================================================================

import argparse
import os
import shlex
import shutil
import subprocess
import sys

SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
CORNERS_NEEDED = {"c1only": [1], "c2only": [2], "c3only": [3], "coverage": [1, 2, 3]}


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


def nonempty(p, minb=2):
    return os.path.isfile(p) and os.path.getsize(p) >= minb


def wl_num(x):
    # Python's float repr ("1e-06") is not valid Mathematica syntax -- ToExpression
    # silently mis-parses it as the symbol `e` (confirmed: ToExpression["1e-06"] -> -6 + e,
    # not 10^-6), which made the chop threshold below a no-op. Mathematica's own
    # scientific-notation marker is "*^", not "e".
    return repr(float(x)).replace("e", "*^")


ap = argparse.ArgumentParser(description="STAGE 2: rawdata_*.m -> compare_*.m (MIN + MEDIAN digit counts), with a tunable chop threshold")
ap.add_argument("--datadir", required=True, help="dir with rawdata_<mode>_w<w>_ordp<>_ordb<>.m (from run_compute.py)")
ap.add_argument("--weights", default="1,2")
ap.add_argument("--modes", default="c1only,c2only,c3only,coverage")
ap.add_argument("--ordp", type=int, required=True)
ap.add_argument("--ordb", type=int, required=True)
ap.add_argument("--chop", type=float, default=1e-4, help="exclude a component from MIN/MEDIAN if |ref| < this (default 1e-4)")
ap.add_argument("--outdir", default="", help="dir for compare_ results (default: --datadir)")
ap.add_argument("--kernel", default=os.environ.get("PIPELINE_KERNEL", ""))
args = ap.parse_args()

KERNEL = resolve_kernel(args.kernel)
DATADIR = os.path.abspath(args.datadir)
OUTDIR = os.path.abspath(args.outdir) if args.outdir else DATADIR
os.makedirs(OUTDIR, exist_ok=True)
weights = [int(s) for s in args.weights.split(",") if s.strip()]
modes = [s.strip() for s in args.modes.split(",") if s.strip()]

print(f"=== ANALYZE (stage 2/2) ===  chop={args.chop}")
print(f"  datadir {DATADIR}\n  outdir  {OUTDIR}")

ok_all = True
for mode in modes:
    for w in weights:
        rawf = os.path.join(DATADIR, f"rawdata_{mode}_w{w}_ordp{args.ordp}_ordb{args.ordb}.m")
        if not nonempty(rawf):
            print(f"  SKIP {mode} w{w}: missing {os.path.basename(rawf)} (run run_compute.py first)")
            ok_all = False; continue
        cmd = KERNEL + [os.path.join(SCRIPTDIR, "03_analyze.wl"),
                         mode, w, args.ordp, args.ordb, DATADIR, OUTDIR, wl_num(args.chop)]
        print(f"\n$ {' '.join(str(c) for c in cmd)}")
        rc = subprocess.call([str(c) for c in cmd])
        if rc != 0:
            print(f"  !! analyze failed for {mode} w{w} (rc={rc})"); ok_all = False

print("\n================================================================")
print(f"  analyze : {'ALL DONE' if ok_all else 'INCOMPLETE -- see above'}")
print("================================================================")
sys.exit(0 if ok_all else 1)
