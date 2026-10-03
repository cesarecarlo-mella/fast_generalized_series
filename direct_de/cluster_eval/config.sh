# ---------------------------------------------------------------------------
# cluster_eval configuration (sourced by every script)
# ---------------------------------------------------------------------------
WEIGHT=${WEIGHT:-6}            # integrate weights 1..WEIGHT (<= 6)
NLAT=${NLAT:-65}               # lattice spacing 1/NLAT; (NLAT-1)(NLAT-2)/2 points per region (65 -> 2016)
NCHUNK=${NCHUNK:-64}           # array tasks per (region, precision)
REGIONS=${REGIONS:-"eucl scat"}   # eucl: x, y > 0 (variables x, y);  scat: x, y < 0, z > 1 (variables u, v)
# requested local error of Bulirsch-Stoer:  |err_i| <= ERR + REL * |f_i|   (REL = 0: paper default, absolute)
ERR_DOUBLE=${ERR_DOUBLE:-1e-12}; REL_DOUBLE=${REL_DOUBLE:-1e-12}
ERR_DD=${ERR_DD:-1e-16};         REL_DD=${REL_DD:-0}
DELTA=${DELTA:-"0.1 0.2"}      # complex deformation of the path (paper: 0.1 0.2)
MAX_TIME=${MAX_TIME:-3600}     # per point [s]; a point that exceeds it is written as FAILED
# toolchain
CXX=${CXX:-g++}
BOOST_INC=${BOOST_INC:-}       # directory that contains boost/ (empty: compiler default search path)
QD_DIR=${QD_DIR:-$HERE/third_party/qd}   # include/qd/*.h and lib/libqd.a (shipped: x86_64 Linux build)
PY=${PY:-python3}
