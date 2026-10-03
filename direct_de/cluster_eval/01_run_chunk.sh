#!/bin/bash
# one chunk of points:   ./01_run_chunk.sh REGION PREC CHUNK      (REGION eucl|scat, PREC double|dd)
# writes out/REGION_PREC_CHUNK.out ; skipped if already complete (resubmission resumes).
set -e
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/config.sh"
R=$1; P=$2; K=$3
chart=xy; [ $R = scat ] && chart=uv
if [ $P = dd ]; then ERR=$ERR_DD; REL=$REL_DD; else ERR=$ERR_DOUBLE; REL=$REL_DOUBLE; fi
IN=$HERE/points/${R}_${K}.txt; OUT=$HERE/out/${R}_${P}_${K}.out
mkdir -p $HERE/out
n=$(grep -c . $IN || true)
if [ -f $OUT ] && [ "$(grep -c '^#' $OUT)" = "$n" ]; then echo "complete: $OUT"; exit 0; fi
export DIFFEQS_MAX_TIME=$MAX_TIME
[ "$REL" != "0" ] && export DIFFEQS_EPS_REL=$REL
$HERE/build/bin/h_${chart}_${P} x $HERE/build/bdy_${R}.txt $ERR $DELTA < $IN > $OUT.tmp
mv $OUT.tmp $OUT
echo "done $OUT ($n points)"
