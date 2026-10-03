#!/bin/bash
# same as slurm_submit.sh on one machine:   NPAR=32 ./run_local.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/config.sh"
NPAR=${NPAR:-$(nproc)}
for R in $REGIONS; do for P in double dd; do for K in $(seq 0 $((NCHUNK-1))); do echo "$R $P $K"; done; done; done \
  | xargs -P $NPAR -L 1 $HERE/01_run_chunk.sh
$PY $HERE/02_analyze.py --jobs $NPAR
