# cluster_eval -- H family, 2016 lattice points per region, double + dd reference (paper method)

Evaluates all 371 canonical masters, weights 1..6, of the H family with the method of
Czakon & Tancredi (arXiv:2606.30354): one Bulirsch-Stoer system for all eps orders **and** the
square root r = Sqrt[xyz] (dr = r db/2b), connection matrix compiled in. Each point is run twice:

* in `double` (production precision);
* in `dd_real` (double-double, QD) as the reference.

`02_analyze.py` compares double against dd for every weight and component, compares both against the
exact weight-1,2 GPLs, and reports the evaluation times.

| region | variables (integration and points) | boundary | lattice |
|---|---|---|---|
| `eucl` | x, y (x, y > 0, x + y < 1) | corner-1 constants, P0 = (0.1, 0.07), r = +Sqrt | x = i/N, y = j/N |
| `scat` | u, v (x = -w/v, y = -u/v, z = 1/v, w = 1-u-v) | CC1 constants (`transport.boundaries/H/dist/bdy_CC1_w*.m`), (x, y) = (-0.1, -0.07), r = -Sqrt\|xyz\| | u = i/N, v = j/N |

Here i, j >= 1 and i + j <= N - 1. With N = 65 that is 2016 points per region. Scattering-region
plots use the axes (u, w).

## Quick start

    ./00_setup.sh          # build (~1 min) + points + boundaries + smoke test; on the login node
    ./slurm_submit.sh      # 2 regions x 2 precisions x NCHUNK array tasks, then the analysis job
    # ... or on a single machine:   NPAR=32 ./run_local.sh

All results are written to `results/` (see below). Resubmitting skips the chunks that are already
complete, so a resubmission resumes where it stopped.

## Configuration (`config.sh`, every value can also be set as an environment variable)

| variable | default | meaning |
|---|---|---|
| `WEIGHT` | 6 | integrate weights 1..WEIGHT |
| `NLAT` | 65 | lattice spacing 1/NLAT (65 -> 2016 points) |
| `NCHUNK` | 64 | array tasks per (region, precision), round-robin split |
| `ERR_DOUBLE`, `REL_DOUBLE` | 1e-12, 1e-12 | double: \|err_i\| <= ERR + REL\|f_i\| |
| `ERR_DD`, `REL_DD` | 1e-16, 0 | dd reference (absolute, as in the paper) |
| `DELTA` | "0.1 0.2" | complex deformation of the path |
| `MAX_TIME` | 3600 | per-point limit [s]; points over the limit are written as FAILED |
| `BOOST_INC` | (compiler default) | directory containing `boost/` (odeint, header only; tested with 1.83) |
| `QD_DIR` | `third_party/qd` | QD headers + `libqd.a` (shipped: x86_64 Linux, Ubuntu 24.04 build) |
| `ARCHFLAGS` | `-march=native` | set e.g. `-march=x86-64-v3` if the login node differs from the compute nodes |
| `TIME_DOUBLE`, `TIME_DD`, `SBATCH_EXTRA`, `NCPU_ANALYZE` | 1h, 4h, -, 16 | SLURM settings, e.g. `SBATCH_EXTRA="--partition=xyz --account=abc"` |

**Double tolerance.** In double the absolute criterion of the paper cannot go below about 1e-10 at
weight 6: the w6 values are O(10-1000), and Bulirsch-Stoer then rejects steps forever. That is why the
default is absolute + relative 1e-12, which is converged and 5x faster than absolute 1e-10.

## Cost (one core per point, weight 6)

| | double | dd_real (1e-16) |
|---|---|---|
| time per point | ~1-3 s (Euclidean), ~2 s (scattering) | ~15-45 s |
| 2016 points | ~1 core-hour | ~15-25 core-hours |

With `NCHUNK=64`, a dd task holds about 32 points and takes 10-25 minutes. `TIME_DD=04:00:00` leaves
plenty of margin.

## Output (`results/`)

* `summary.md`: per region:
  * timing (median, mean, 10%, 90%, max, total core-hours, evaluations, ms per evaluation);
  * double vs dd for each weight: worst component per point (worst point, 1%, 10%, median), the
    median component, the fraction of points with >= 12, 10 and 8 digits, and the maximum absolute
    error on the vanishing components;
  * double and dd vs exact at weights 1 and 2. The dd line certifies the reference.
* `maps_<region>.pdf/png`: double vs dd, one column per weight, rows worst and median component.
* `exact_<region>.pdf/png`: weight 2 vs exact, for double and dd.
* `per_function_<region>.pdf/png`: per component, the median and worst point (w2 vs exact; w6 double vs dd).
* `timing.pdf/png`: histograms of the time per point.
* `digits.npz`: all per-point / per-weight / per-component digit arrays, coordinates and times.
* `ref_<region>.npz`: cache of the exact references.

Digits of a component are -log10(|a-b|/|b|), using only components with |b| >= 1e-4. The 136
functions that vanish identically at weights 1-2, and any other zeros, are covered by the
absolute-error column instead. A NaN or inf in a result counts as 0 digits.

## Layout

    config.sh  00_setup.sh  01_run_chunk.sh  slurm_submit.sh  run_local.sh  02_analyze.py
    src/          diffeqs.hpp (solver, interface of the paper), h_family_c.cpp (driver; -DUV_CHART: u, v),
                  gen/hgen_data.{cpp,hpp} (connection matrices as exact hi/lo doubles), gen/hgen_letters.hpp,
                  patch_boost.py (explicit double casts in boost's bulirsch_stoer.hpp for dd_real; algorithm unchanged)
    data/         dd boundaries, weights 0-6 (from frobenius.cpp: Euclidean corner 1; scattering CC1, written in (u, v))
    exact/        exact_w2.py + eps_1.m, eps_2.m (exact weight-1,2 GPL solution)
    tools/        make_points.py, truncate_bdy.py, plot_w6_paper.py (plot style)
    third_party/  QD 2.3.x (BSD-LBNL licence, see COPYRIGHT): headers + static libqd.a

Requirements:
* C++17 compiler;
* boost headers (odeint);
* python3 with numpy, mpmath and matplotlib.

If the shipped `libqd.a` does not link on your cluster, because the C++ ABI or glibc are too old,
build QD 2.3.24 from https://www.davidhbailey.com/dhbsoftware/ and set `QD_DIR` to its install
prefix.

**Important.** Keep `-ffp-contract=off`. FMA contraction breaks the double-double arithmetic.

**End points on letter zeros.** A lattice point lying exactly on a letter's zero set (an apparent
singularity of the DE, where the connection is infinite) cannot be an end point. `make_points.py`
moves such points by 1e-6/N and prints them. N = 65 has none; for N = 9, 8 points per region are
moved.

## Regenerating the inputs

These are only needed if the DE or the boundaries change. Work from `series.all/direct_de/paper_method`:
* `gen_code.py` writes `gen/`;
* `frobenius.cpp` and `export_data.py` write the boundaries.
