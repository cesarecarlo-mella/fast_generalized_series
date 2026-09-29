#!/bin/bash
# ---------------------------------------------------------------------------
#  One node, the point (x, y) = (1/2, 0), weights 1..6, order 20/20.
#  run_array.py runs -n prime jobs side by side and adds primes by itself if
#  a weight needs more. Resubmitting the same script resumes where it stopped.
#  Needs inputs/bdy_c0_w1.m .. bdy_c0_w6.m (boundary constants at the point).
# ---------------------------------------------------------------------------
#SBATCH --job-name=pt_w6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=pt_%j.log

HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-$HERE/venv/bin/python}       # a venv with numpy + gmpy2 (see README)

$PY "$HERE/run_array.py" -w 6 -n "$SLURM_CPUS_PER_TASK" --ordt 20 --ordy 20 --workdir "$HERE/dist_20_20"
