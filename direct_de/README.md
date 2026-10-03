# direct_de — direct numerical solution of the canonical DE (H, Euclidean region)

Strategy of Czakon–Tancredi (arXiv:2606.30354): no series, integrate dJ = eps A J
numerically along complex-deformed paths.  The boundary is the shuffle-regularised
corner-1 boundary (bdy_c1_w<w>.m, Log x, Log y -> 0).

## Files
- `common.py`   : load atilde / boundary vectors (numeric and exact rational).
- `letters.py`  : the 14 H letters; ray Taylor series of the dlogs (FFT); dlogs on a complex path
                  with Sqrt[xyz] tracked continuously.
- `corner.py`   : corner-1 start.
    * shuffle regularisation: near (0,0)  J = sum log^a x log^b y C_ab + O(x, y),
      C_ab = Ax^a Ay^b c_{w-a-b}/(a! b!), Ax = a_1 (letter x), Ay = a_2 (letter y), built EXACTLY.
      Checked against the x^0 y^0 log part of the old corner-1 series (w1-4).
      x+y (a_6) and x-y (a_11 + a_12) annihilate all C_ab: spurious at the corner.
    * along a ray x = s X, y = s Y: regularised constants c_ray = sum C_ab log^a X log^b Y (log s -> 0),
      then a log-power (Frobenius) series in s to s0 (N = 60 terms).
- `solver.py`   : `Solver(inputs, bdy_dir, ray=(1, 0.7), s0=0.1).evaluate(x, y)` -> weights 0..6.
    * path z_k(t) = P0_k + (t + 4 i delta_k t(1-t)) (x_k - P0_k), delta = (0.1, 0.2)
    * weight recursion with 20-point Gauss-Legendre panels and a spectral cumulative-integration matrix,
      panels doubled until two refinements agree.
    * matrix balancing (atilde entries up to 1.6e5 -> ~150).
- `test_c1.py`, `test_twostart.py`, `frob_test.py`: tests.

## Tests (double precision)
| test | median digits | worst component |
|---|---|---|
| vs old corner-1 series, w1-4, near the corner | 13.3-14.9 | 7.3-10.9 |
| vs exact GPL solution, w1-2, across the triangle | 12.7-14.7 | 7.3-10.4 |
| two different starts (ray, s0), w1-6, across the triangle | 12.2-14.7 | 6.1-11.6 |

The worst components are limited by double precision: rows of atilde have entries up to ~1e5 that cancel
in sums over letters (balancing does not help for that), so small components lose ~5-8 digits. This is
the same floor as in the paper (1e-6 - 1e-9 in double, quad precision for more).

Time: ~15-25 s per point (pure numpy, 2 cores), dominated by the 256-panel refinement
the error criterion asks for; start-up (exact C_ab + Frobenius) ~30 s once.

## Precision map, weights 1-2 (`prec_grid.py`, `plot_prec.py` -> `direct_de_precision_w2.{pdf,png}`)
Paper grid grid_N2000.m (2016 points), direct solver (start: corner 1, x = start point P0) vs the exact GPL solution:

| | worst | 10% | median | best |
|---|---|---|---|---|
| w1, min over components | 7.9 | 9.6 | 10.3 | 11.1 |
| w1, median over components | 14.1 | 14.4 | 14.6 | 15.0 |
| w2, min over components | 5.1 | 7.3 | 8.2 | 9.5 |
| w2, median over components | 12.6 | 13.0 | 13.4 | 14.1 |

Uniform over the whole triangle (no region problem, unlike the series); the floor is double-precision
cancellation. About 0.4 s per point at weight <= 2 (2016 points: 13 min on 2 cores).
`evaluate(x, y, wmax=2)` integrates only up to the requested weight.
