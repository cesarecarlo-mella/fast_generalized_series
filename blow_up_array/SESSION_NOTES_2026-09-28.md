# Session notes 2026-09-28: array/modular blow-up pipeline

(The same notes are saved in the claude.ai project "Fast_series" as `claude/blow_up_array_pipeline.md`.)

## Folders (all in series.all/)

- **`blow_up_array/pipeline/`**: H family, production, self-contained.
  `venv/bin/python run_array.py -w 6 -n 60 -c 1,2,3 --ordt 20 --ordy 60 --workdir dist_20_60 --phys 20`
- **`all_families_array/`**: 21 families.
  `venv/bin/python run_all_families.py -n 64 --ordt 15 --ordy 45 -w 6 --phys 15`.
  The letters are 20 letters at 15/45, so 15/45 is the maximum order.
- **`blow_up_array_fast/pipeline/`**: experimental, on hold. The integrability shortcut makes
  the chain 2.15x faster, bit-identical.

## Validated (exact)

- **H at 10/30:** weights 1-4 match `blow_up_form`. Weights 1-5 match the cluster phys
  files (x, y <= 10).
- **H at 20/60:** weights 1-2 match the cluster phys files (x, y <= 20).
- **Families:** M3 weights 1-4 and N3 weights 1-2 match the FORM results. I3 runs through
  weights 1-6.

## Fixed today

- **Cache-key bug** in `ModEngine.T`: junk in the truncation-edge slots from weight 6 on.
  The fix is in all local copies. The cluster still needs the new `blowup_array.py`, and
  H weight 6 should be rerun to clean the extra terms out of `exact_*_w6.m`.
- **Missing bdy file** is now an error (it used to be silently zero).
- **Parser** accepts `cnst[6, k][5]`.
- **I3 weight-6 boundaries** rebuilt from the corrected `all.ints/boundaries_all/I3_bc.m`.
  Only the 103 cnst components changed. The old files are kept as `*.old_with_cnst`.
  The other families' `_bc.m` files are not newer than their `dist/` builds.

## Cluster (newton1, /scratch/ge45tuy)

- **Python:** always use `venv/bin/python`. `all_families_array/venv` is a symlink to
  `blow_up_array/pipeline/venv`.
- **Files to copy:** `run_bernoulli.py`, `05_bernoulli.wl`, `06_merge_bernoulli.wl`, and the
  fixed `blowup_array.py` / `run_array.py` / `00_prepare.py`.
- **Bernoulli, per weight:**
  `./run_bernoulli.py -w $w -n 32 -c 1,2,3 --ordp 20 --ordb 20 --workdir dist_20_60 --physdir dist_20_60 --outdir dist_20_60`.
  It needs Mathematica.

## Performance (one prime chain, weights 1-6)

| order | time | peak memory |
|---|---|---|
| 10/30 | 54 s | 1.0 GB |
| 15/45 | 130 s | 1.9 GB |
| 20/60 | 356 s | 3.1 GB |

H at 20/60 took about 35 min per corner on the cluster.

Ideas not yet done: per-slot common denominators (about half the primes), the P/Q recurrence
for rational letters, triangular truncation, and 60-bit primes via FLINT.
