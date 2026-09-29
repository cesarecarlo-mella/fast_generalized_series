#!/bin/bash
# ---------------------------------------------------------------------------
#  Many nodes: one SLURM ARRAY TASK PER PRIME (the prime jobs are fully
#  independent), then one reconstruction job per weight that waits for them.
#
#     ./slurm_multinode.sh <corner> <W> <ORDt> <ORDy> [NPRIMES] [FIRST_K]
#  (corner is always 0 here: the point (x, y) = (1/2, 0))
#
#  e.g.  ./slurm_multinode.sh 0 6 20 20            # primes k = 0..39
#  If a reconstruction log says NEED_MORE_PRIMES, add primes and
#  reconstruct again:
#        ./slurm_multinode.sh 0 6 20 20 10 40      # primes k = 40..49
#  (reconstruction always uses every prime found in workdir/residues/)
# ---------------------------------------------------------------------------
set -e
C=$1; W=$2; OT=$3; OY=$4; NP=${5:-40}; K0=${6:-0}
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${SRC:-$HERE/inputs}
WD=${WD:-$HERE/dist_${OT}_${OY}}
PY=${PY:-python3}
mkdir -p "$WD/logs"

# step 0: parse inputs once (fast, seconds)
if [ ! -f "$WD/setup_c$C.pkl" ]; then
  $PY "$HERE/00_prepare.py" --corner $C --ordt $OT --ordy $OY --wmax $W \
      --inputs "$SRC" --letters "$SRC" --workdir "$WD"
fi

# step 1: one array task per prime
K1=$((K0 + NP - 1))
JID=$(sbatch --parsable --job-name=pt1_c$C --array=$K0-$K1 \
      --cpus-per-task=1 --mem=4G --time=04:00:00 \
      --output="$WD/logs/slurm_01_c${C}_k%a.txt" \
      --wrap="$PY $HERE/01_chain.py --workdir $WD --corner $C --k \$SLURM_ARRAY_TASK_ID --wmax $W")
echo "submitted prime array $JID (k=$K0..$K1)"

# step 2: reconstruct each weight once all primes are done
for w in $(seq 1 $W); do
  sbatch --job-name=pt2_c${C}_w$w --dependency=afterok:$JID \
      --cpus-per-task=16 --mem=64G --time=04:00:00 \
      --output="$WD/logs/slurm_02_c${C}_w$w.txt" \
      --wrap="$PY $HERE/02_reconstruct.py --workdir $WD --corner $C --weight $w --jobs 16"
done
echo "submitted reconstruction jobs for weights 1..$W"
