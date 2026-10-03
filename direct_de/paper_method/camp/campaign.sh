#!/bin/bash
# $1 = A|B (half of the points; two workers run in parallel, one per core)
cd "$(dirname "$0")/.."
export DIFFEQS_MAX_TIME=900
h=$1
E6=data/bdy_w6_dd_ray1.txt; S6=data/bdy_w6_dd_CC1ray1.txt
E2=data/bdy_w2_dd_ray1.txt; S2=data/bdy_w2_dd_CC1ray1.txt
run() { echo "$(date +%T) start $1" >> camp/log_$h; }
# --- Euclidean, weight 6, 500 points
run e6_dbl10; ./h_family_c_double x $E6 1e-10 0.1 0.2 < camp/eucl500_$h.txt > camp/e6_dbl10_$h.out
run e6_dbl12; ./h_family_c_double x $E6 1e-12 0.1 0.2 < camp/eucl500_$h.txt > camp/e6_dbl12_$h.out
run e6_dd14;  ./h_family_c_dd     x $E6 1e-14 0.1 0.2 < camp/eucl500_$h.txt > camp/e6_dd14_$h.out
# --- Euclidean, weight 2, 2000 points
run e2_dbl12; ./h_family_c_double x $E2 1e-12 0.1 0.2 < camp/eucl2000_$h.txt > camp/e2_dbl12_$h.out
# --- scattering, weight 6, 500 points
run s6_dbl10; ./h_family_c_double x $S6 1e-10 0.1 0.2 < camp/scat500_$h.txt > camp/s6_dbl10_$h.out
run s6_dbl12; ./h_family_c_double x $S6 1e-12 0.1 0.2 < camp/scat500_$h.txt > camp/s6_dbl12_$h.out
run s6_dd14;  ./h_family_c_dd     x $S6 1e-14 0.1 0.2 < camp/scat500_$h.txt > camp/s6_dd14_$h.out
# --- scattering, weight 2, 2000 points
run s2_dbl12; ./h_family_c_double x $S2 1e-12 0.1 0.2 < camp/scat2000_$h.txt > camp/s2_dbl12_$h.out
# --- weight 2 in dd (2000 points each region)
run e2_dd16;  ./h_family_c_dd     x $E2 1e-16 0.1 0.2 < camp/eucl2000_$h.txt > camp/e2_dd16_$h.out
run s2_dd16;  ./h_family_c_dd     x $S2 1e-16 0.1 0.2 < camp/scat2000_$h.txt > camp/s2_dd16_$h.out
run end; echo done > camp/done_$h
