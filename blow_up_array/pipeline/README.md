# blow_up_array/pipeline -- cluster pipeline for the array/modular blow-up recursion

This replaces `blow_up_form_2/pipeline/{01_recur_step, 02_reassemble, 03_to_physical_vars,
04_merge_physical}` and produces the same `exact_c<c>_w<w>.m` and `phys_c<c>_w<w>_ord<ORD>.m`
files. It needs only **python3, numpy and gmpy2**: no Mathematica, no FORM.
The inputs are the existing `atilde.m`, `sol0.m`, `bdy_c<c>_w<w>.m` and `dlog_d{t,y}_c<c>.m`.
The master letters (25/75) can be used directly, because terms above the requested order are
dropped on load.

    pip install --user numpy gmpy2

**Self-contained:** `inputs/` holds everything needed: `atilde.m`, `sol0.m`, `bdy_c{1,2,3}_w{1..6}.m`,
and the master letters `dlog_d{t,y}_c{1,2,3}.m` at 25/75, so any order up to 25/75 works.
Copy this `pipeline/` folder to the cluster and run it. For weight 7 or higher, add
`bdy_c<c>_w7.m` etc. to `inputs/`.

## Steps

| script | what | parallelism |
|---|---|---|
| `00_prepare.py` | parse the inputs once per corner -> `setup_c<c>.pkl` | 1 process, seconds |
| `01_chain.py --k K` | the whole chain for weights 1..W, **modulo prime number K** -> `residues/J_c<c>_w<w>_k<K>.pkl` | one independent job per prime |
| `02_reconstruct.py --weight w` | exact rationals from all primes found, checked against 2 extra primes -> `exact_c<c>_w<w>.m`, `rec_c<c>_w<w>.pkl`, and with `--phys ORD` also `phys_c<c>_w<w>_ord<ORD>.m`. Exits with code 3 (`NEED_MORE_PRIMES`) if more primes are needed. | `--jobs N` |
| `run_array.py` | single-node driver (like `run_weight.py`): runs 00, then 01 on `-n` cores, then 02, and adds primes automatically until every weight verifies | `-n` |
| `validate.py` | exact comparison with reference `exact_*.m` (any higher order) and/or `phys_*.m` | |

Parallelism is **over primes**. Each prime job is a complete, independent run of the
recursion, with nothing shared and no communication. The number of prime jobs you need grows
with weight, not with the number of components:

| weight | primes needed (order 10/30, measured) |
|---|---|
| 1-5 | 20 |
| 6 | 26 |

The default starting count is 4W+6, and more are added automatically if needed.
Primes are indexed by `--k`, so on a cluster the array task id is the prime.

## Run

Single node (`slurm_single_node.sh`):

    ./run_array.py -w 6 -n 32 -c 1,2,3 --ordt 15 --ordy 45 --workdir dist_15_45 --phys 15

Many nodes (`slurm_multinode.sh`, one array task per prime, then one reconstruction job per
weight):

    ./slurm_multinode.sh 1 6 15 45          # corner 1, weights 1..6, primes k=0..39

Everything can be resumed. Finished primes and weights are skipped, and a later `-w 7`
restarts every prime from its stored weight-6 residues. Use a new workdir for a new order.

## Measured cost (corner 1, one prime job, one core, 2-core cloud VM)

| weight | order 10/30 | order 15/45 |
|---|---|---|
| 1 | 0.5 s | 1.5 s |
| 2 | 1.4 s | 3.1 s |
| 3 | 3.1 s | 6.9 s |
| 4 | 7.0 s | 17.2 s |
| 5 | 14.4 s | 37.3 s |
| 6 | 27.0 s | 63.7 s |
| **chain 1..6** | **54 s, 1.0 GB peak** | **130 s, 1.9 GB peak** |

Total CPU time is roughly (number of primes) x (chain time). At 15/45 and weight 6 that is
about 30 x 130 s, or roughly 1 CPU-hour per corner, plus reconstruction (minutes, parallel).
For comparison, the FORM pipeline at 15/45 used about 72 CPU-min for weight 3 alone and about
150 cores x 60 min for weight 6 (both for corner 1).

## Validated

Every comparison was exact, coefficient by coefficient, and all were identical:

- corner 1, order 10/30, weights 1-4, against `blow_up_form/pipeline/dist/exact_c1_w*.m`
- corner 1, order 10/30, weights 1-5, against `cluster_results/blow_up_form_2_order_20/out/phys_c1_w*_ord20.m`
  in the window x, y <= 10
- corners 2 and 3, weights 1-2, against the `blow_up_form` exact files, using the 25/75 master letters
- the files the pipeline writes (`exact_*.m`, `phys_*.m`) read back to exactly the same content

## How it works

See `blowup_array.py`, in short:

- **Series storage.** Each series is stored as dense (t, BIG) grids, one per
  (parity, log powers, constant monomial).
- **Multiply.** `rowVec.J = sum_l dlog_l (sum_j atilde_ij,l J_j)`, computed as two dense
  float64 BLAS matrix products.
- **Integrations.** These act directly on the grids.
- **Arithmetic.** Everything is done modulo primes below 2^22. Every rational is a single
  machine integer, and the BLAS partial sums stay below 2^53, so the result is exact.
- **Recovering rationals.** CRT plus rational reconstruction, with a second pass that uses
  the common denominator of each component.
