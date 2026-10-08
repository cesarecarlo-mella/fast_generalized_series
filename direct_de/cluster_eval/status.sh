#!/bin/bash
# progress of cluster_eval:  ./status.sh   (or: watch -n 30 ./status.sh)
cd "$(dirname "$0")"
source ./config.sh 2>/dev/null
echo "queue:"; squeue -u $USER -h -o "%T %j" | sort | uniq -c
for R in $REGIONS; do for P in double dd; do
  f=${R}_${P}
  n=$(cat out/${f}_*.out out/${f}_*.out.tmp 2>/dev/null | grep -c '^#')
  c=$(ls out/${f}_*.out 2>/dev/null | wc -l)
  tot=$(grep -c . points/${R}.txt 2>/dev/null)
  printf "%-12s %5d / %s points   %4d / %d chunks\n" $f $n "$tot" $c $NCHUNK
done; done
echo "failed points: $(cat out/*.out* 2>/dev/null | grep -c FAILED)"
[ -f results/summary.md ] && echo "DONE -> results/summary.md"
