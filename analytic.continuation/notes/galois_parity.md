# The sign of √(xyz) as a Galois (ℤ₂) symmetry of the canonical DE

H family (3L1M), analytic continuation to x, y < 0. Written 3 October 2026.
The LaTeX version of the same note is `galois_parity.tex` / `.pdf` in this folder.

## Setup

The 371 H masters satisfy dg = ε dA g, with A = Σ_k a_k log l_k(x, y). The alphabet has one
square root, r = √(xyz) with z = 1 − x − y, so it lives in the quadratic extension K(r) of
K = ℚ(x, y).

- **Even letters** l1–l6, l9–l20 are rational (in K).
- **Odd letters** are l7 = (xz − i r)/(xz + i r) and l8 = (xy − i r)/(xy + i r).

## The Galois action

Gal(K(r)/K) = {1, σ} with σ: r → −r.

- Even letters are fixed by σ.
- Odd letters go to their inverses, so dlog l7 → −dlog l7 and dlog l8 → −dlog l8.

The masters split the same way. Masters 156 and 157 (the ones with a square-root prefactor) are
odd; the other 369 are even. Let P be diagonal with −1 at 156 and 157 and +1 elsewhere.

## The matrix is graded

Write A in even/odd blocks (A_ee, A_eo, A_oe, A_oo). I checked this entry by entry on
`H_atilde`:

- The diagonal blocks (even–even, odd–odd) contain only even letters.
- All 134 nonzero odd↔even entries contain only l7 and l8.
- l7 and l8 appear nowhere else.

So **σ(A) = P·A·P**.

The odd–odd entries, for example A_156,156, do **not** change sign. They multiply an odd master
and feed the derivative of an odd master, so the two signs cancel.

## Consequence

If g solves dg = ε dA g with boundary value g0, then P·g solves the conjugate equation
d g′ = ε dσ(A) g′ with boundary value P·g0. The proof is one line:
d(Pg) = ε P dA g = ε (P dA P)(Pg) = ε dσ(A)(Pg), using P² = 1.

The odd boundary constants vanish at every weight, including weight 0 (`sol0`). This holds at
the old corners c1, c2, c3 and at the continued points CC1, CC2, CC3. So P·g0 = g0, and flipping
the sign of r changes only the overall sign of g156 and g157. The other 369 masters are
identical.

Block by block:

- **Even rows:** dg_e = ε (dA_ee g_e + dA_eo g_o). A flipped odd letter multiplies a flipped odd
  master, so nothing changes.
- **Odd rows:** d(−g_o) = ε (dA_oo (−g_o) + (−dA_oe) g_e), which is the original equation times −1.

**Symbol level.** g^(w) = ∫ dA g^(w−1) + c^(w), and the odd constants are zero. By induction on
the weight, every symbol term of an even master has an even number of odd letters, and every
term of an odd master an odd number. The masters are eigenvectors of σ with eigenvalue ±1.

## The continued region: a minus sign in front of the root, not under it

The chart is x = −w/v, y = −u/v, z = 1/v with u + w + v = 1. It covers x, y < 0, z > 1.

- **The radicand is positive at both ends:** xyz > 0 for x, y > 0, and xyz = uw/v³ > 0 in the
  continued region. It is negative in between (for example x < 0 < y), where r is imaginary.
- **Each crossing contributes a phase e^(±iπ/2).** With the continuation used for the CC
  boundaries, log x → log|x| + iπ and log y → log|y| + iπ:
  √x → i√|x| and √y → i√|y|, so **r → −√|xyz|**.
- Opposite phases for x and y would give +√|xyz|. So the sign of r is tied to the iπ choice,
  which comes from the i0 prescription.

**Choice adopted:** r = −√|xyz| in the continued region.

- `cc_series/letters/build_cc_letters.py` now defaults to `--odd-sign -1` (`--odd-sign +1`
  gives the principal root).
- Compared with the principal root, only l7 and l8 flip in the CC2/CC3 letter files.
- `check_map.py` (exact comparison with the old corner letters, with the sign applied to l7 and
  l8) still finds every coefficient identical.
- This matters only when g156 and g157 are combined with quantities defined in the original
  region, such as the prefactor √(xyz) itself.
