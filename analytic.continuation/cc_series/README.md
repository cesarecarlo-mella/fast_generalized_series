# cc_series — series expansions covering the continued triangle (H family)

Continued region: x, y < 0, z > 1, parametrised by the triangle
x = -w/v, y = -u/v, z = 1/v, u + w + v = 1 (u, w, v > 0).
In the original variables X = |x| = w/v, Y = |y| = u/v.
Continuation convention: Log x -> log|x| + iπ, Log y -> log|y| + iπ, √(xyz) -> -√|xyz|.

## Why seven expansions

The corners CC2 (w → 1) and CC3 (u → 1) are LINES in the original variables
(x → -∞ at any y, and y → -∞ at any x). From weight 2 the masters depend on the
ratio y = -u/v (resp. x = -w/v) there. The corner expansion splits this into
log-powers times towers in the ratio, and the towers converge only on one side
of |ratio| = 1. So each corner needs two charts, one for each side, and each
chart's boundary constants must be defined with the same order of limits as
the chart.

| chart | blow-up | variables | covers | boundary constants |
|---|---|---|---|---|
| CC1  | none (old corner-1 series, continued) | X, Y | X + Y < 1 (u + w < v) | old c1 + continuation (`tools/cc1_from_c1.py`) |
| CC2  | u = t v², BIG = v (u ≪ v) | u, v | Y < 1, near w = 1 | `transport.boundaries` (u → 0 first), all weights |
| CC3  | w = t v², BIG = v (w ≪ v) | w, v | X < 1, near u = 1 | `transport.boundaries` (w → 0 first), all weights |
| CC2e | v = t u², BIG = u (v ≪ u) | u, v | Y > 1, near w = 1, reaches the edge v = 0 | weights 1-2 only (see below) |
| CC3e | v = t w², BIG = w (v ≪ w) | w, v | X > 1, near u = 1, reaches the edge v = 0 | weights 1-2 only (see below) |
| CC2r | none: x = -1/p, y = -1 - p + t | t, p | band around Y = 1 near w = 1 | weights 1-2 only (see below) |
| CC3r | none: y = -1/p, x = -1 - p + t | t, p | band around X = 1 near u = 1 | weights 1-2 only (see below) |

The towers in CC2/CC2e converge only for |ratio| ≠ 1, because there is a singularity at
ratio −1 (y = +1). That leaves a band around Y = 1 (and around X = 1). The ratio-line charts
are centred on that band, at x = −∞, y = −1. The letter 1 − x(1+y) passes through that point;
it is a coordinate line there (= t/p) and is spurious: the pipeline produces no Log[t] at all
(checked at weights 1-2). Their output is `exact_CC2r_w<w>.m` in (t, p). There is no phys
conversion: evaluate in t, p, with p = 1/X and t = 1 − Y + p. CC2r is good for |t| < 0.5,
p < 0.35. Its weight-2 constants contain I Pi Log[2] (y = −1 gives 1 − y = 2).

Selection rule (`check/coverage2.py`):

- If X + Y < 1, use CC1. Its log(1 − x − y) limits it to |x| + |y| < 1.
- Otherwise, on the CC2 side (w ≥ u): use CC2r if |t| < 0.5 and p < 0.35; else CC2 if Y < 1; else CC2e.
- On the CC3 side, the mirror of this.

## Test (weights 1-2, whole triangle, 406-point grid; `results_w2/`)

| order | points with min digits ≥ 2 | ≥ 4 | ≥ 6 |
|---|---|---|---|
| N = 20, seven charts, rule (`cov7_N20`) | w1: 406/406, w2: 391/406 | w1: 400, w2: 348 | w1: 359, w2: 282 |
| N = 20, seven charts, best chart per point | w2: 403/406 | w2: 368 | w2: 293 |
| N = 12, five charts, no ratio charts (`cov5_N12`) | w1: 404/406, w2: 356/406 | w1: 340, w2: 235 | w1: 206, w2: 98 |

Remaining weak spots:

- The middle of the Y = 1 and X = 1 bands, where p > 0.35 and the ratio charts no longer converge.
- The centre (u ≈ w ≈ 0.3).

Both are the analogue of the old centre, which Bernoulli variables fixed. Adding more ratio points
along the band, or Bernoulli-type variables, would help.

## Layout

- `letters/build_cc_letters.py` writes the dlog series of the ORIGINAL letters composed
  with each chart, exact, into `inputs/`. Run it as `--corners CC2,CC3,CC2e,CC3e`.
- `letters/check_map.py` checks CC3 exactly against the old c3 letters through the alphabet map.
- `pipeline/` is blow_up_array/pipeline generalised to named regions (`-c CC2,CC3,CC2e,CC3e`).
  It uses the original H `atilde.m` unchanged.
- `tools/cc1_from_c1.py` writes `phys_CC1_*` from `phys_c1_*` (exact continuation).
- `check/exact_continued.py` gives the EXACT weight ≤ 2 solution continued to any point of the region. It
  starts from the GPL solution `eps_{1,2}.m` at (x, y) = (0.1, 0.1), rotates both phases by +π near the
  origin, then integrates the canonical DE exactly along a complex path in the triangle coordinates.
  A ± imaginary bump gives the same result (to 1e-10), so there is no singularity of the masters inside the region.
- `check/fit_r.py` does the same for CC2r/CC3r (basis Pi^2, Log[2]^2, I Pi Log[2]).
- `check/rcharts.py`, `check/coverage2.py`: evaluation of the ratio charts and the seven-chart coverage test.
- NOTE: the check scripts contain the cloud paths used for the test runs; adapt them before running them elsewhere.
- `check/fit_consts.py` fits the CC2e/CC3e weight-1,2 constants from the continued exact solution, at two points each,
  and recognises them exactly (q·iπ, q·π², denominators dividing 4320).
- `check/coverage.py` evaluates every chart against the continued exact solution on a grid of the whole triangle.

## Run

    cd pipeline
    python3 run_array.py -w 6 -n 32 -c CC2,CC3 --ordt 20 --ordy 60 --workdir dist_20_60 --phys 20
    python3 run_array.py -w 2 -n 32 -c CC2e,CC3e --ordt 20 --ordy 60 --workdir dist_20_60 --phys 20   # weights 1-2 only for now
    python3 run_array.py -w 2 -n 32 -c CC2r,CC3r --ordt 20 --ordy 20 --workdir dist_r_20              # no --phys (no blow-up)

## Open

The CC2e/CC3e/CC2r/CC3r boundary constants for weights ≥ 3 need a transport with the right order of limits
(v ≪ u at CC2, v ≪ w at CC3). The weight ≤ 2 ones were fitted from the exact solution, which is fine
for a test but not for production. Options:

- Transport along the exceptional divisor of the CC2/CC3 corner (ratio s = u/(u+v) from 0 to 1). The
  odd letters give √(s(1-s)) there, and rationalising gives letters {0, ±i}.
- Transport along the edge v = 0 from a point where the constants are known.
