# all_families_array -- array/modular recursion for every integral family

This does the same job as `all_families/` (`run_all_families_recursion.py` plus each family's
`run_exact_all_weights.py` and `run_to_physical.py`), using the `blow_up_array` pipeline:
python3 with numpy and gmpy2, no Mathematica and no FORM for the recursion.
The folder is self-contained, so copy it to the cluster as it is.

    code/                  pipeline (same as blow_up_array/pipeline)
    letters/               master letters shared by all families: 20 letters, built to 15/45
    <Family>/inputs/       atilde.m, sol0.m, bdy_c{1,2,3}_w{1..6}.m   (from all_families/<F>/pipeline/dist)
    <Family>/dist_T_Y/     results: exact_c<c>_w<w>.m, phys_c<c>_w<w>_ord<ORD>.m, logs/

There are 21 families: A B1 B2 C D E1 E2 F1 F2 G1 G2 G3 G4 I1 I2 I3 L3 M3 N3 P3 R2.

## Run

    python3 -m venv venv && venv/bin/pip install numpy gmpy2      # once
    venv/bin/python run_all_families.py -n 64 --ordt 15 --ordy 45 -w 6 --phys 15
    venv/bin/python run_all_families.py --dry-run --ordt 15 --ordy 45    # show the order
    venv/bin/python run_all_families.py -n 64 --ordt 15 --ordy 45 --only M3,N3
    venv/bin/python run_all_families.py -n 64 --ordt 15 --ordy 45 -- --nprimes 64   # extra run_array.py options

Families run cheapest first. A failure is reported and the next family still runs. Everything
can be resumed. To run a single family by hand:

    venv/bin/python code/run_array.py -w 6 -n 64 --ordt 15 --ordy 45 --workdir M3/dist_15_45 \
        --inputs M3/inputs --letters letters --phys 15

## Differences from the H pipeline (blow_up_array)

- **Letters.** There are 20 letters instead of 14, identical for all families, built at
  **15/45**, so any order up to 15/45 works. For a higher order, rebuild them with
  `all_families/pipeline_master/00_build_letters.wl` and put the new files in `letters/`.
- **Component count** varies by family, from about 80 for M3 to about 330 for I3.
- **Boundaries** are taken from each family's `all_families/<F>/pipeline/dist`. I3's weight-6
  files were stale, containing 103 `cnst[6,k][5]` placeholders per corner, and were rebuilt on
  28 Sep from the corrected `all.ints/boundaries_all/I3_bc.m`; the old files are kept as
  `*.old_with_cnst`. If a `*_bc.m` changes again, rebuild without Mathematica:
  `python3 code/build_bdy_from_bc.py --bc ../../all.ints/boundaries_all/<F>_bc.m --atilde <F>/inputs/atilde.m --check <F>/inputs --out <F>/inputs`
  (`--check` shows which files change before they are overwritten).
- **Bernoulli.** Run `code/run_bernoulli.py` per family with `--physdir <F>/dist_T_Y`. It needs
  Mathematica; see the blow_up_array notes.

## Validated

Every comparison was exact, coefficient by coefficient, against the FORM results already in
`all_families/`:

- M3, corners 1-3, weights 1-4 (order 4/12): all identical.
- N3, corners 1-3, weights 1-2 (order 4/12): all identical.
- I3 (the largest family), corner 1, weights 1-6 at 6/18: runs through and verifies (with the old
  boundaries). There is no reference result to compare it against.
