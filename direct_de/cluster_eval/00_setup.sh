#!/bin/bash
# build the four solvers, write the lattice points + chunks, truncate the boundaries to WEIGHT.
# Run once on a login node (or interactively):   ./00_setup.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/config.sh"
B=$HERE/build; mkdir -p $B/bin $B/obj
$PY $HERE/src/patch_boost.py $B "$BOOST_INC"
INC="-I$B/boostpatch ${BOOST_INC:+-I$BOOST_INC} -I$QD_DIR/include"
FLAGS="-std=c++17 -O3 -ffp-contract=off -DNDEBUG -DBOOST_UBLAS_NDEBUG ${ARCHFLAGS:--march=native}"
# -ffp-contract=off is REQUIRED: FMA contraction breaks QD's double-double arithmetic.
# ARCHFLAGS: set e.g. ARCHFLAGS="-march=x86-64-v3" if login and compute nodes differ.
echo "compiling the connection data (static, ~1 min) ..."
$CXX -std=c++17 -O1 -c -o $B/obj/hgen_data.o $HERE/src/gen/hgen_data.cpp
for chart in xy uv; do
  D=""; [ $chart = uv ] && D="-DUV_CHART"
  $CXX $FLAGS $INC $D -DREAL=double -o $B/bin/h_${chart}_double $HERE/src/h_family_c.cpp $B/obj/hgen_data.o &
  $CXX $FLAGS $INC $D -DDIFFEQS_QD -DREAL=dd_real -o $B/bin/h_${chart}_dd $HERE/src/h_family_c.cpp $B/obj/hgen_data.o $QD_DIR/lib/libqd.a &
done
wait
ls -la $B/bin
$PY $HERE/tools/make_points.py $HERE/points $NLAT $NCHUNK $REGIONS
$PY $HERE/tools/truncate_bdy.py $HERE/data/bdy_w6_dd_ray1.txt        $B/bdy_eucl.txt $WEIGHT
$PY $HERE/tools/truncate_bdy.py $HERE/data/bdy_w6_dd_CC1ray1_uv.txt  $B/bdy_scat.txt $WEIGHT
# smoke test: one point per region and precision
echo "0.3 0.2" | $B/bin/h_xy_double x $B/bdy_eucl.txt 1e-10 0.1 0.2 --quiet
echo "0.3 0.2" | $B/bin/h_uv_double x $B/bdy_scat.txt 1e-10 0.1 0.2 --quiet
echo "setup done"
