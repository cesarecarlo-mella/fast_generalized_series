# Boundary constants for the continued-region patches: how to get them analytically

Region: x, y < 0, z > 1. The chart is x = −w/v, y = −u/v, z = 1/v, with u + w + v = 1.

Continuation: Log x → log|x| + iπ, Log y → log|y| + iπ, and √(xyz) → −√|xyz|.

The corners CC2 (w → 1) and CC3 (u → 1) are **lines** in the original variables: x → −∞ at any y, and y → −∞ at any x.
From weight 2 the masters depend on the ratio there (y = −u/v at CC2, x = −w/v at CC3). Because of that,
**a corner constant (the log⁰·log⁰ coefficient) depends on the order of limits.** Each patch's constants must be
computed in the same ordering as the patch's chart.

Below, everything is written for CC2. CC3 is the same with u ↔ w (x ↔ y).

## 1. Status

| patch | chart | ordering | constants | status |
|---|---|---|---|---|
| CC1 | old corner-1 series, continued | — | exact continuation of `phys_c1` | all weights |
| CC2 | u = t v², BIG v | u ≪ v | `transport.boundaries` (along u = 0) | all weights |
| CC3 | w = t v², BIG v | w ≪ v | `transport.boundaries` (along w = 0) | all weights |
| CC2e | v = t u², BIG u | v ≪ u | fitted from the exact solution | **w ≤ 2 only** |
| CC3e | v = t w², BIG w | v ≪ w | fitted | **w ≤ 2 only** |
| CC2r | x = −1/p, y = −1 − p + t | p → 0 at y = −1 | fitted | **w ≤ 2 only** |
| CC3r | y = −1/p, x = −1 − p + t | p → 0 at x = −1 | fitted | **w ≤ 2 only** |

## 2. The divisor transport (gives CC2e and CC2r at every weight)

All three CC2-type patches sit on the **exceptional divisor** of the CC2 corner. Blow it up once:

    u = s ρ ,   v = (1 − s) ρ ,   w = 1 − ρ ,      s = u/(u+v) ∈ [0, 1]   (s = y/(y − 1))

- ρ → 0 is the corner.
- s labels the ratio: s → 0 is u ≪ v (CC2), s = 1/2 is y = −1 (CC2r), and s → 1 is v ≪ u (CC2e).

**Restricted DE.** Take every original letter composed with this chart, expand at ρ = 0, and keep only the
s-dependence. The term m_k log ρ is dropped (regularised), and the constant parts of the letters drop out of the dlog:

    d/ds F(s) = ε Σ_k a_k  ∂_s log l_k(s, ρ = 0)  F(s)

**Even letters.** On the divisor the even letters factor into s, 1 − s and constants. Check this with the
letter-builder: run `build_cc_letters.py` with a W = 1 chart u = s ρ and read off the ρ⁰ part. If the check
holds, F is a combination of HPLs/GPLs with letters {0, 1}, built weight by weight with
G(a, w⃗; s) = ∫₀ˢ G(w⃗; t)/(t − a) dt. This is the same GIntegrate recursion as `00_transport_boundaries.wl`.

**Odd letters (l7, l8).** On the divisor l8 → (q − i)/(q + i) with q = √(s/(1 − s)), so its dlog is
−i ds/√(s(1 − s)) and is not rational in s. Rationalise:

    s = q²/(1 + q²)   (q ∈ [0, ∞)),   or   s = sin²θ

In q the alphabet becomes {0, i, −i}. The odd masters 156 and 157 start at 0 (zero constants everywhere), but they
pick up arcsin-type terms along the divisor. Those feed the even masters at weight ≥ 2 through products
(odd letter) × (odd master).

**Start (s → 0).** The CC2 constants from the transport are the regularised s → 0 values. This holds if the towers
vanish at s → 0, which they do: they are power series in u/v with no constant term.

> Caveat: CC2 uses the W = 2 chart (u = t v², i.e. u ≪ v²), while the W = 1 divisor reaches s → 0 with
> v² ≪ u ≪ v. They agree only if the tangent curve u ~ c v² (the reason for W = 2) is spurious. The
> corner-1-continued and CC3 agreement supports this, but it must be checked: transport in both orderings at
> weight ≤ 2 and compare with the fitted constants.

**End points:**

- **CC2e (s → 1).** Substitute s = 1 − τ (`fibrtto1.m`-style rules, or the q-variable version) and drop
  G(0,…,0; τ), i.e. log(1 − s) → 0. The constants come out in MZVs: π², ζ3, ζ5, iπ powers, and possibly level-4
  constants (Catalan) from the {±i} letters. In q, s → 1 is q → ∞: use the q → 1/q map.
- **CC2r (s = 1/2, q = 1).** The values are G(…; 1/2) and G(±i, …; 1): log 2, π, Li_n(1/2), and at weight 2
  iπ·log 2 (seen in the fit). Then convert from the divisor normalisation (log ρ → 0) to the chart normalisation
  (log p → 0). On the divisor p = −1/x = v/w = (1 − s)ρ/(1 − ρ), so log p = log ρ + log(1 − s) + O(ρ). At s = 1/2:

      Σ_k log^k ρ F_k   with   log ρ = log p + log 2   ⇒   constant_p = Σ_k (log 2)^k F_k

  The coordinate t: the letter 1 − x(1+y) = t/p is spurious (no Log[t], checked at w ≤ 2), so there is no t-regularisation.

**Cross-checks available:**

- At w ≤ 2 the result must reproduce `inputs/bdy_CC2e_w{1,2}.m` and `bdy_CC2r_w{1,2}.m`. Those were fitted from the
  exact solution and recognised exactly: CC2e is q·iπ, q·π²; CC2r adds iπ·log 2; all denominators divide 4320.
- CC2e and CC3e (and CC2r/CC3r) are related by u ↔ w, which permutes the masters. Same check for CC2 vs CC3.

## 3. Numerical route (any weight, no new analytics)

`direct_de/` now integrates the full canonical DE, weights 0–6, from the corner-1 boundaries. With the path
continued into x, y < 0 exactly as in `check/exact_continued.py` (rotate both phases by +π near the origin),
it gives J at any continued point and any weight. Then, for each patch:

1. Run the pipeline with zero boundary at weight w (lower weights already exact) and evaluate the patch series
   at two points well inside its convergence domain.
2. c_w = J_direct(P) − series(P; c_w = 0). The two points must agree to the solver's accuracy.
3. Recognise c_w exactly with PSLQ in the weight-w basis. For CC2e/CC3e that is MZVs with iπ; for CC2r/CC3r add
   log 2 and Li_n(1/2) (as in `check/fit_r.py`, which does weight 2).

This needs the direct solver at higher precision than double for weights ≥ 4: the worst components are only good
to about 1e-8 in double.

## 4. Order of limits: summary rule

The constant of a patch is the log⁰·log⁰ coefficient in **its own** chart's ordering. Translating between orderings
means transporting along the divisor in s, because the difference is the asymptotics of the ratio functions
(for example Li₂(−v/u): 0 vs −π²/6 − ½ log²(v/u)). Changing the log variable at the same ordering
(log v vs log p vs log ρ) only re-expands the logs (log p = log ρ + const) and is purely algebraic.
