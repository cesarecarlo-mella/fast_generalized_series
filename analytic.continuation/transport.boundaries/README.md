# transport.boundaries

Computes boundary constants at three new kinematic regions — **CC1**, **CC2**,
**CC3** — by analytically continuing OUT from the existing corner-1 boundary
data of a family's differential equation.

## What this is

CC1, CC2, CC3 are all reached starting from the *same* known point (corner
1's own boundary, already available from `all.ints/boundaries_all` / the
`all_families/` recursion output), walked along three different straight-line
paths in a reparametrized kinematic chart:

- **CC1** — direct continuation from corner 1 itself. No differential
  equation is solved: the existing `phys_c1_w<w>_ord<N>.m` series is
  re-expressed through the substitution `Log[x] -> Log[w]-Log[1-w-u]+I*Pi`
  (same idea for `Log[y]`), Taylor-expanded in the new `(u,w)` around 0, and
  the surviving logs are dropped. This picks up the `I*Pi` monodromy phase
  from continuing `x`/`y` across their branch cut — nothing else changes.
- **CC2**, **CC3** — reached via the "high energy" reparametrization
  `x -> -w/(1-u-w)`, `y -> -u/(1-u-w)`: an effective 1D differential equation
  is derived along `w` (`u=0` fixed, → CC2) and along `u` (`w=0` fixed, →
  CC3), integrated weight-by-weight with `GIntegrate`, seeded at every weight
  by CC1's own value (since `u=w=0` in this parametrization IS corner 1),
  then evaluated at the path's far endpoint via the project's standard
  change-of-fibration-basis (`change_of_fibrationBasis/fibrtto1.m`,
  `1-t` → `t->0`, drop the resulting divergent all-zero-string `G[0,...,0,t]`).

**CC2 and CC3 are *not* the analytic continuation of the existing blow-up
corners C2/C3 used by `../../all_families/`.** They share the "2"/"3" labels
by coincidence of numbering, not by relation — hence the double-`C` naming
(`CC1`/`CC2`/`CC3`) throughout this pipeline's output, to keep them visually
and unambiguously distinct from `all_families/.../dist/bdy_c<c>_w<w>.m`.

## What this is explicitly NOT

**This performs analytic continuation only — it does NOT perform the
translation from the Euclidean region to the physical decay region.**

Everything here is a purely mathematical continuation of the differential
equation's solution within its own variable space (picking up monodromy
phases, re-expanding around new points). It says nothing about, and does not
implement, the physics-level map from the Euclidean kinematic point where the
DEQ/recursion is naturally solved to the physical (decay) kinematic region —
that is a separate, distinct step, not performed anywhere in this pipeline.
Do not treat CC1/CC2/CC3 output as physical/decay-region results without that
additional, separate translation.

## Layout

```
transport.boundaries/
  README.md               (this file)
  <Family>/                one folder per family (currently: H only)
    00_transport_boundaries.wl
    dist/
      bdy_CC1_w0.m .. bdy_CC1_w6.m
      bdy_CC2_w0.m .. bdy_CC2_w6.m
      bdy_CC3_w0.m .. bdy_CC3_w6.m
```

Each `bdy_CC<n>_w<w>.m` is a length-`Nf` list of pure constants (`Pi`,
`Zeta[n]`), matching the weight-`w` boundary at that region — same
`dist/`/naming convention as `all_families/.../pipeline/dist/bdy_c<c>_w<w>.m`.

## Usage

```
cd <Family>/
wolframscript -file 00_transport_boundaries.wl        # weights 1..6 (default)
wolframscript -file 00_transport_boundaries.wl 1       # weights 1..1 only (cheap sanity run)
```

The weight-cutoff argument is **positional**, not a `--flag` — `wolframscript
-file` silently drops `--`-prefixed arguments instead of forwarding them to
the script, so a flag-based interface doesn't work here.

## Prerequisites

- `../../all.ints/a_tildes/<Family>_atilde.m`, `letters.dict.m`
- `../../all.ints/boundaries_all/<Family>_bc.m` (for `sol0`, the corner-independent weight-0 boundary)
- `../../change_of_fibrationBasis/fibrtto1.m`
- `../../cluster_results/blow_up_form_2_order_20/out/phys_c1_w<w>_ord20.m`,
  weights 1-6 — the CC1 step needs this family's *own* corner-1 physical
  series already computed at real order. Currently only exists for H
  (from the cluster). Generalizing to another family needs that family's own
  `phys_c1_w<w>_ord<N>.m` produced first (via `all_families/`), before this
  script can run for it.

## Correctness checks already performed (H)

Before writing `00_transport_boundaries.wl`, both non-trivial routines were
independently verified, not just eyeballed:

- **`fibrtto1.m` regularization completeness**: all 126 rules (complete
  `{0,1}`-word alphabet, weight 1 through 6) reduce to a pure constant under
  the regularized `t->0` limit — checked via `GetGs` on every rule's RHS
  after substitution, zero leftover `G[...]` or `t`-dependence at any weight.
- **`pathW`/`pathU` (the effective transport connections)**: cross-checked
  against an independent, naively-constructed derivative for all 20 letters.
  18/20 match exactly; the 2 sqrt-sector exceptions (l7, l8) were verified via
  an actual numerical limit (not the symbolic `Series` shortcut) — both real
  and imaginary parts of `d(Log l7,l8)/dw` genuinely shrink to 0 as `u->0`.

## Performance

H (`Nf=371`), all 6 weights, all three regions, **no parallelism** (single
Mathematica kernel): ~9.7 minutes total on a laptop. CC1 (the direct
continuation) is the expensive part (~513s, dominated by `Series`/`Expand` on
the large order-20 physical series); CC2/CC3 (`GIntegrate` along the
transported DEQ) are cheap (~36s combined). This step is expected to stay
laptop-feasible even for the largest families — the real bottleneck is
producing each family's own `phys_c1_w<w>_ord<N>.m` in the first place, which
needs the full (cluster-run, parallelized) `all_families/` recursion.

## Planned next steps (not yet built)

- Generalize `00_transport_boundaries.wl` to every family (needs each
  family's own `phys_c1_w<w>_ord<N>.m` first).
- A new recursion pipeline that actually integrates the series *in* the CC1,
  CC2, CC3 regions (this folder currently produces boundary constants only,
  not full series).
