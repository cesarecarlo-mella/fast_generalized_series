#!/bin/bash
# Submit everything:  one SLURM array per (region, precision), NCHUNK tasks each (1 core per task),
# then the analysis job (multi-core) once all arrays are done.
#     ./00_setup.sh            # once
#     ./slurm_submit.sh
# Resubmitting skips completed chunks.  Time per point at weight 6 (one core):
#     double ~1-3 s,   dd ~20-50 s   ->  with NCHUNK=64: ~32 points/task, dd tasks ~15-30 min.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/config.sh"
mkdir -p $HERE/logs
DEPS=""
for R in $REGIONS; do
  for P in double dd; do
    T=${TIME_DOUBLE:-01:00:00}; [ $P = dd ] && T=${TIME_DD:-04:00:00}
    JID=$(sbatch --parsable --job-name=ce_${R}_${P} --array=0-$((NCHUNK-1)) \
          --cpus-per-task=1 --mem=1G --time=$T ${SBATCH_EXTRA} \
          --output=$HERE/logs/${R}_${P}_%a.log \
          --wrap="$HERE/01_run_chunk.sh $R $P \$SLURM_ARRAY_TASK_ID")
    echo "submitted $R $P: array $JID"
    DEPS="$DEPS:$JID"
  done
done
JA=$(sbatch --parsable --job-name=ce_analyze --dependency=afterany$DEPS \
     --cpus-per-task=${NCPU_ANALYZE:-16} --mem=32G --time=02:00:00 ${SBATCH_EXTRA} \
     --output=$HERE/logs/analyze.log \
     --wrap="$PY $HERE/02_analyze.py --jobs \${SLURM_CPUS_PER_TASK}")
echo "submitted analysis: $JA (after all arrays)"
