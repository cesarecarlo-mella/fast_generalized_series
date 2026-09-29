#!/bin/bash
# ---------------------------------------------------------------------------
#  One node, all corners, weights 1..6  --  the simplest way to run.
#  run_array.py runs -n prime jobs side by side and adds primes by itself if
#  a weight needs more. Resubmitting the same script resumes where it stopped.
#
#  Memory: about 2 GB per prime job at order 15/45 and weight 6 (measured),
#  about 3-4 GB at 20/60 (estimated), so --mem >= ncores * that.
# ---------------------------------------------------------------------------
#SBATCH --job-name=bua_w6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=bua_%j.log

# module load python/3.11          # <-- adapt to your cluster
# pip install --user numpy gmpy2   # once

# sbatch runs a spooled copy of this script, so locate the pipeline via the
# submit directory:  cd blow_up_array/pipeline && sbatch slurm_single_node.sh
HERE=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}
SRC=${SRC:-$HERE/inputs}   # atilde.m, sol0.m, bdy_*, dlog letters (shipped in ./inputs)

python3 "$HERE/run_array.py" -w 6 -n "$SLURM_CPUS_PER_TASK" -c 1,2,3 \
    --ordt 15 --ordy 45 --workdir "$HERE/dist_15_45" \
    --inputs "$SRC" --letters "$SRC" --phys 15
