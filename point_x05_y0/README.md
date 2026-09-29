# point_x05_y0 — H family expanded around (x, y) = (1/2, 0)

This is the array/modular pipeline (a copy of `blow_up_array/pipeline`) applied to a new
expansion point:

    x = 1/2 + t ,   y = y            (no blow-up)

The letters diverge only in y (`Log[y]`, and `Sqrt[y]` in l7, l8). They are regular in x,
so the result is a double series in `t = x - 1/2` and `y`, with `Log[y]` and half-integer
powers of `y`. There are no `Log[t]` terms.

The files use "corner 0" as the label of this point (`dlog_dt_c0.m`, `bdy_c0_w<w>.m`,
`exact_c0_w<w>.m`, ...), so all the scripts work unchanged. `-c 0` is the default everywhere.

## Task 1: is a blow-up needed? No.

`python3 letters/check_point.py` lists every curve on which a dlog of the 14 active letters
(1–13, 16) is singular and checks which of them pass through the point. That means the zeros
of each letter's polynomial, plus, for l7 and l8, the zeros of `b` and of `a^2 + b`.

- Only `y = 0` passes through the point: `l2 = y`, and the square-root argument
  `b = x(1-x-y) y` of l7 and l8. For l8 it also enters through `a^2 + b = y·x(xy+z)`.
- In every case the letter is `y^1` times a factor that is nonzero at the point.
- With a single singular curve there is nothing to resolve, so no blow-up is needed. Every
  dlog is `y^-1` or `y^-1/2` times a double power series in `(t, y)`.
- The script also checks `alphabet.py` against `letters.dict.m`; they agree to 1e-15.

Nearest singularities, which bound where the series converge:

| direction | radius | set by |
|---|---|---|
| t at y = 0 | 1/2 | x = 0 and x = 1 (l1, l3, l4, l6, l10–l13, l16, l7) |
| y at x = 1/2 | 1/4 | l12 = x² − x + y and l13 = (x−1)² − y, which both vanish at (1/2, 1/4) |

So along `x = 1/2` the series converge for `y < 1/4`. Away from `x = 1/2` the region shrinks.
In the physical triangle it is roughly bounded by the l12 curve `y = x(1−x)` and the l13 curve
`y = (1−x)²`.

## Task 2: letters to order 20

`python3 letters/build_letters.py 20 20` writes `pipeline/inputs/dlog_dt_c0.m`,
`dlog_dy_c0.m` and `info_c0.m`. It runs in under a second.

- **Arithmetic.** Exact rationals (gmpy2) with a truncation of `t^i`, `i <= 20` and `y^e`,
  `e <= 20`. For l7 and l8, `e` is a half-integer.
- **Plain letters.** For `P = y^m Q`: `d_t log P = Q_x/Q` and `d_y log P = m/y + Q_y/Q`.
- **l7, l8.** `d_v log((a − I√b)/(a + I√b)) = I (2 a_v b − a b_v) / (√b (a² + b))`, written as
  `I y^(n−d−1/2) · (power series)`. The square root is the principal one, as in
  `letters.dict.m`, so no sign fixing is needed.
- **Numeric check.** Against numerical derivatives of the letters themselves (mpmath,
  40 digits) at (t, y) = (1/50, 3/100). The maximum relative difference is 3e-19, which is
  the truncation error.
- **Integrability check.** `python3 pipeline/check_de.py --ordt 8 --ordy 8 --wmax 4` checks,
  weight by weight, that the y-integration constant `−∂_y ∫A_t J dt + A_y J` of the recursion
  does not depend on t. This holds only if every letter expansion and atilde are consistent,
  and it does not need the boundary constants.
  - Result: 0 violations out of about 5×10^5 checked coefficients, at weights 1–4.
  - A deliberately perturbed l7 is caught at weight 3.

## Task 3: the pipeline

The scripts are the same as in `blow_up_array/pipeline` (python3, numpy, gmpy2; no
Mathematica). The recursion is unchanged:

    J_w(t, y) = ∫_0^t A_t J_{w-1} dt  +  ∫_0^y [ −∂_y(first term) + A_y J_{w-1} ]_{t-independent} dy  +  bdy_w

In words:
- Integrate in t from the point at fixed y.
- Then integrate in y along x = 1/2, starting from y = 0. There, `y^-1 Log[y]^k` integrates
  to `Log[y]^(k+1)/(k+1)` with no constant.
- Add the boundary constants `bdy_w`.

    cd pipeline
    ./run_array.py -w 6 -n 32 --ordt 20 --ordy 20 --workdir dist_20_20

or `sbatch slurm_single_node.sh`, or `./slurm_multinode.sh 0 6 20 20`.

Output: `dist_20_20/exact_c0_w<w>.m`, a list of 371 series in `t`, `y`, `Log[y]` with the
constants kept symbolic. Here `t` means `x − 1/2`. No change of variables is needed, so there
is no `phys_` step.

### Boundary conditions

The boundary constants come from the solutions on the line y = 0:
`.../H/solutions/eps_sep/w<w>/sol<w>_1.m`. These are HPLs `G[a..., t]` with `t = x` and
`a ∈ {0, 1}`.
- This line is the right one: at weight 1 the coefficients of `G[0,t]` and `G[1,t]` in
  `sol1_1.m` equal the `Log[x]` and `Log[1−x]` coefficients predicted by atilde·sol0, for
  all 371 components.
- `bdy_c0_w<w>.m` is `sol<w>_1.m` at `t = 1/2`, written by:

      python3 pipeline/make_bdy.py --sols <.../H/solutions/eps_sep>

- The values `G[a1, ..., an, 1/2]` are kept exact, as symbolic atoms next to `Pi` and
  `Zeta[n]`. The recursion carries every constant monomial as a label, so nothing needs to be
  evaluated. There are only 98 constant monomials even at weight 6.
- `pipeline/gpl_num.py` evaluates the atoms numerically (40 digits, series in x). A
  PolyLogTools reduction to `Log[2]`, `PolyLog[n, 1/2]`, … is optional.
- In the full expansion around (1/2, 0), `bdy_w` is the `y^0 Log[y]^0` coefficient of
  `J_w(1/2, y)`.

### Checks against known solutions

- **`check_x.py` (line y = 0, weights 1–4, order 10/10).** Compares the `y^0 Log[y]^0` part of
  our series with `sol<w>_1.m` at `x = 1/2 + t`, for t = ±0.05 and 0.1. They agree to at most
  2e-7, which is the order-10 truncation error `(2t)^11` at t = 0.1.
- **`check_exact.py` (full 2D, weights 1–2).** Compares with the exact solution
  `full_sol/eps_sep/eps_<w>.m`, as in `solution.nb`. The GPLs in x and y are evaluated
  numerically at (x, y) = (0.53, 0.02), (0.47, 0.03), (0.55, 0.05) and (0.45, 0.01), for all
  371 components. This includes every `Log[y]` and `Sqrt[y]` term.

  | order | weight 1 | weight 2 |
  |---|---|---|
  | 10/10 | 8e-12 | 2e-10 |
  | 20/20 | 1e-21 | 3e-20 |

### Cost

One prime chain, weights 1–6 at order 20/20, with zero boundary (so fewer constant tags than
the real run): 13 s and 0.34 GB on one core. For comparison, H at 20/60 took about 6 min per
prime. The grid is 22×22 instead of 22×62, so the real run should be much faster than the H
corners: about 30 primes × a few minutes.

## Files

    letters/alphabet.py        active letters as exact polynomials; the point
    letters/series2.py         exact truncated double series (mul, inv, sqrt)
    letters/build_letters.py   TASK 2: dlog series -> pipeline/inputs/
    letters/check_point.py     TASK 1: singular curves through the point, radii
    letters/letters.dict.m     copy of a_tildes/letters.dict.m (for the check)
    pipeline/                  TASK 3: 00_prepare, 01_chain, 02_reconstruct, run_array,
                               slurm scripts, modcmp
    pipeline/make_bdy.py       bdy_c0_w<w>.m from eps_sep/w<w>/sol<w>_1.m at t = 1/2
    pipeline/gpl_num.py        numerical HPLs G[a..., x], a in {0,1} (40 digits)
    pipeline/check_de.py       integrability test (no boundary needed)
    pipeline/check_x.py        check on y = 0 against sol<w>_1.m
    pipeline/check_exact.py    full check against full_sol/eps_sep/eps_<w>.m (w <= 2)
    pipeline/inputs/           atilde.m, sol0.m (same as H), dlog_d{t,y}_c0.m, info_c0.m,
                               bdy_c0_w1..6.m
