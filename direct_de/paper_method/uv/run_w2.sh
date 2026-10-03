#!/bin/bash
cd "$(dirname "$0")/.."
export DIFFEQS_MAX_TIME=300
( ./h_family_c_uv_double x data/bdy_w2_dd_CC1ray1_uv.txt 1e-12 0.1 0.2 < uv/lat.txt > uv/s2_dbl12.out; echo done > uv/done_s2 ) &
( ./h_family_c_double x data/bdy_w2_dd_ray1.txt 1e-12 0.1 0.2 < camp/eucl2000.txt > uv/e2_dbl12.out; echo done > uv/done_e2 ) &
wait
