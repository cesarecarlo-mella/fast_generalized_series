#!/bin/bash
cd "$(dirname "$0")/.."
export DIFFEQS_MAX_TIME=300
B=data/bdy_w6_dd_CC1ray1_uv.txt
for h in A B; do (
  ./h_family_c_uv_double x $B 1e-10 0.1 0.2 < uv/w6pts_$h.txt > uv/s6_dbl10_$h.out
  ./h_family_c_uv_double x $B 5e-11 0.1 0.2 < uv/w6pts_$h.txt > uv/s6_dbl511_$h.out
  ./h_family_c_uv_dd     x $B 1e-14 0.1 0.2 < uv/w6pts_$h.txt > uv/s6_dd14_$h.out
  echo done > uv/done_w6_$h ) & done
wait
