# compare_results_w3: series against the exact solution at weight 3

This folder is a copy of `compare_results/`, the pipeline behind the paper figures. The only
change is that the comparison against the exact result is done at **weight 3 (ε³)** instead of
weight 2. The weight-2 folder is unchanged.

## Exact reference

`exact_solutions/H_solution_w3.m` is a copy of `../H_exact_sol/H_solution_w3.m`:
- `HW3[i] = {ε⁰, ε¹, ε², ε³}` for the 371 canonical masters `BASISLC[H,"",i]`, Euclidean region.
- Only Log, PolyLog[2|3] of rational arguments, π², ζ₃.
- Master 157 at ε³ uses Ti₂, Ti₃ (Ti_n(u) = Im Li_n(i u), defined in the file).

Checked here: its ε² part agrees with the weight-2 reference used so far (`eps_2.m`, as evaluated
by Mathematica in `results_w2/rawdata_*`) to 4e-14 absolute at 150 random grid points.

## Two ways to run it (same output files)

### A. Python, fast (what produced `cluster_results/results_w3_consistent`)

```bash
venv/bin/python fast_exact.py \
    --datadir ../cluster_results/blow_up_array_results/out \
    --weight 3 --ceiling 20 --orders 12 16 20 \
    --modes c1only,c2only,c3only,coverage \
    --grid grids/grid_N2000.m --exact exact_solutions/H_solution_w3.m \
    --outdir ../cluster_results/results_w3_consistent -n 8
```

The definitions are identical to `results_w2_consistent`, the data of the paper figures:
- **Order N**: the order-20 phys_/ber_ files truncated to both corner variables ≤ N
  (as in `chop_phys_ber.wl`).
- **Difference**: series_N − exact, signed.
- **Digits**: min(16, −log10 |diff|/|ref|).
- **Chop**: components with |ref| < 1e-4 are dropped.
- **Statistics**: MIN and MEDIAN over components.
- **Modes**: c1only, c2only, c3only, and coverage (nearest corner, the same `pickC` rule, ties
  included).
- **Grid**: grid_N2000, 2016 points.

New files:
- `exact_eval.py`: evaluates `HW3[i][[w+1]]` with mpmath at 30 digits, at the **exact rational**
  grid points.
- `hp_eval.py`: evaluates the series with exact (double-double) coefficients and longdouble sums.
  Without this, double precision caps the result at about 11–14 digits where the true agreement is
  16.
  - This needs x86-64 (cluster or cloud). On arm64 it stops with an error.
- `fast_exact.py`: the driver. It writes `rawdata_*` and `compare_*` in the usual format.

**Validation.** Run at `--weight 2`, it reproduces `results_w2_consistent`, the paper data, for
c1only and coverage at N = 12, 16, 20. The digits agree to within 0.005 at all 2016 points, and
the corner assignment is identical.

Two details matter:
- **The grid point must be exact.** Some components are very steep: shifting x by 2e-17 moves
  one of them by 6e-15. So x, y, z must come from the same rationals for the series and for the
  exact solution.
- **Ties.** `pickC` takes the last maximum when coordinates tie, e.g. on x = y.

### B. Mathematica, the original path

`run_multi_order.py`, `run_compute.py`, `01_eval_grid.wl`, `02_merge_raw.wl` and `03_analyze.wl`
now accept weight 3, which is the default in this folder. `01_eval_grid.wl` loads
`H_solution_w3.m` and uses `HW3[i][[4]]`. Weights 1, 2 and 6 behave as before. Example:

```bash
./run_multi_order.py --datadir ../cluster_results/blow_up_array_results/out \
    --ceiling 20 --orders 12,16,20 --weights 3 \
    --modes c1only,c2only,c3only,coverage --npoints 2000 -n 8 \
    --chop 1e-4 --outdir ../cluster_results/results_w3
```

## Figures (same scripts, `--weight 3`)

```bash
python3 plot_w6_paper.py --indir ../cluster_results/results_w3_consistent --outdir figures/paper \
    --weight 3 --pair 20,20 --compare-pairs 12,12 16,16 20,20 --modes c1only
python3 plot_grid.py --indir ../cluster_results/results_w3_consistent --outdir figures/paper \
    --weight 3 --mode coverage --rep bernoulli --pairs 12,12 16,16 20,20
```

The outputs are `figures/paper/w3_exact_orders_c1only_min.pdf` and
`figures/paper/w3_exact_coverage_bernoulli_grid.pdf`, the weight-3 versions of paper figures 1
and 2.

## Note on what is stored

`cluster_results/results_w3_consistent/` contains:
- the 12 `compare_*_w3_*.m` files (c1only, c2only, c3only and coverage × N = 12, 16, 20);
- `cache/exact_w3_grid_N2000.npz`, the exact ε³ values on the grid.

The `rawdata_*_w3_*.m` files, about 42 MB each, were not copied back. Re-running
`fast_exact.py` with the same arguments regenerates them, together with the parsed series tables,
in about 11 minutes on 2 cores.
