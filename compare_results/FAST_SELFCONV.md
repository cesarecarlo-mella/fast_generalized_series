# Fast weight-6 self-convergence check: `fast_selfconv.py`

This replaces `chop_phys_ber.wl` + `00_build_coeffcache.wl` + `01_eval_grid.wl` +
`02_merge_raw.wl` + `03_analyze.wl` for the self-convergence case (weight 6, low vs high
order). It writes the **same** `rawdata_*.m` and `compare_*.m` files, so `plot_compare.wl`
works unchanged. No Mathematica is needed, only python3 with numpy and gmpy2 (the
`blow_up_array` venv works).

```bash
venv/bin/python fast_selfconv.py --datadir <dir with phys_c*_w6_ord20.m, ber_c*_w6_ordp20_ordb20.m> \
    --weight 6 --ceiling 20 --pairs 18,19 19,20 --modes c1only,c2only,c3only,coverage \
    --grid results/grid_N1000.m --outdir results_w6_fast -n 6
```

## How it works

1. **Parse once.** Each **ceiling** file (6 files: 3 corners, phys and ber) is parsed once
   into a numeric monomial table, cached as `.npz` in `<outdir>/cache`. The files are
   parsed in parallel. Lower orders are exponent masks, not separate files.
2. **One evaluation pass.** Every monomial is evaluated on the whole grid in a single
   vectorised pass and binned by its level, max(deg_SMALL, deg_BIG). Every order N is a
   cumulative sum over the bins, so **all order pairs and all modes come from one pass**.
3. **Direct differences.** low − high is summed directly from the dropped levels. This
   avoids the ~1e-11 round-off floor the old path had, because it no longer subtracts two
   large totals.
4. **Same analysis as `03_analyze.wl`**: the digit formula, `--chop` (default 1e-4), and
   MIN/MEDIAN.

## Measured (H weight 6, corner 1, grid_N1000 = 990 points, cloud VM)

| | old Mathematica path | fast_selfconv.py |
|---|---|---|
| one mode, one pair | 50 workers × ~1900 s ≈ **27 CPU-hours** | — |
| parse (once, cached) | — | 14 min (1 core per file) |
| evaluation, **all 21 orders, all pairs** | — | ~3.3 min per corner (1 core) |

## Checked against the old `results_w6` rawdata (c1only, 18/19)

- **Bernoulli:** agrees to ~1e-12 (median relative difference).
- **Native (phys):** the old files turned out to use the **old y-only (BIG) truncation**,
  with x left at order 20. That is exactly the flaw `chop_phys_ber.wl`'s header says makes
  the native check blind to non-convergence in x.
  - With `--phys-trunc big` the tool reproduces the old numbers: max relative difference
    2e-5, and only where the old values sit at their round-off floor. The median digits
    agree to 0.01.
  - The default, `--phys-trunc both`, is the corrected definition, and its native results
    differ from the old ones. **The old weight-6 native self-convergence plots should be
    regenerated.**
