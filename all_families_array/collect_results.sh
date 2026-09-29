#!/bin/bash
# ---------------------------------------------------------------------------
#  collect_results.sh -- pack ONLY the result files of every family into one
#  archive (no residues/, no setup pickles, no logs, no rec_*.pkl), ready to
#  download.
#
#     ./collect_results.sh dist_15_45                 -> results_dist_15_45.tar.gz
#     ./collect_results.sh dist_15_45 phys ber        -> only phys_* and ber_*
#  (default kinds: exact phys ber)
# ---------------------------------------------------------------------------
set -e
WD=${1:?usage: collect_results.sh <workdir name, e.g. dist_15_45> [kinds...]}
shift || true
KINDS=${@:-exact phys ber}
cd "$(dirname "$0")"
LIST=$(mktemp)
for F in */"$WD"; do
  for k in $KINDS; do
    ls "$F"/${k}_c*_w*.m 2>/dev/null || true
  done
done > "$LIST"
N=$(wc -l < "$LIST")
[ "$N" -gt 0 ] || { echo "no result files found under */$WD"; exit 1; }
OUT=results_${WD}.tar.gz
tar -czf "$OUT" -T "$LIST"
rm -f "$LIST"
echo "packed $N files into $(pwd)/$OUT  ($(du -h "$OUT" | cut -f1))"
