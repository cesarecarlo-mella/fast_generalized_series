# blow_up_array -- slot/array prototype of the blow-up recursion

Same maths and truncation as `blow_up_form*/pipeline/01_recur_step*.wl` + `02_reassemble.wl`,
but without symbolic expressions, Mathematica or FORM:

- Every series is a set of dense (t, BIG) coefficient grids. There is one grid per tag
  (parity_t, parity_BIG, Log[t] power, Log[BIG] power, constant monomial in Pi/Zeta/I).
  Exponents are counted in half-units, and the grid starts at -1 so the dlog poles fit.
- The product is computed through the letters:
  `rowVec.J = sum_l dlog_l * (sum_j atilde_ij,l J_j)`.
  That is one dense matrix product over components, followed by one dense
  shift-matrix product per letter (both float64 BLAS).
- The integrations and `-D[., BIG]` are closed-form maps acting directly on the grids.
- All arithmetic is done modulo primes below 2^22. Every rational becomes one machine
  integer, sums need no gcds, and BLAS partial sums stay below 2^53, so the arithmetic is
  exact. Exact rationals are recovered once at the end by CRT plus rational
  reconstruction, with a second pass that uses the component's common denominator.
  Each coefficient is checked against 2 extra primes.

Run (corner 1, order 10/30, weights 1..6, compared with the blow_up_form exact files
and the cluster phys files):

    python3 run.py --data ../blow_up_form/pipeline/dist --ref ../blow_up_form/pipeline/dist \
       --phys-ref ../cluster_results/blow_up_form_2_order_20/out \
       --corner 1 --ordt 10 --ordy 30 --wmax 6 --nprimes 20 --jobs 4

Requires python3, numpy and gmpy2 (`pip install gmpy2`).
See `run_c1_ord10_w1-6.log` for the first test run.
