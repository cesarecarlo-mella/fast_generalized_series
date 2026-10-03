# paper_method -- the H family with the method of Czakon & Tancredi (arXiv:2606.30354)

One first-order system for all eps orders **and** the square root r = Sqrt[x y z]
(dr = r db / 2b), integrated with `boost::numeric::odeint::bulirsch_stoer` along the
complexified straight path of the paper (eq. 4.1, delta = (0.1, 0.2)), local error =
max absolute error.  Precision: `double` or QD `dd_real`.  Full description: `doc/method.tex`.

## Files

| file | what |
|---|---|
| `diffeqs.hpp` | generic solver, interface of the paper (`diffeqs<Real>`, `connection_t`, `vector_field_t`, `path_t`, boundary stream, `evaluate`, exceptions); optional compiled mat-vec hook |
| `h_family_c.cpp` | **production driver**: connection matrix compiled in (`gen/hgen_data.cpp`), letters generated (`gen/hgen_letters.hpp`) |
| `h_family.cpp`, `h_data.hpp` | same, matrix read at run time from `data/atilde.txt` |
| `h_family_gen.cpp` | fully unrolled generated code (one statement per entry) -- correct but 1.4-2.7x SLOWER (8-13 MB of instructions per call); kept for reference |
| `gen_code.py` | writes `gen/` (all variants) from `atilde.m`, coefficients as exact (hi, lo) double pairs |
| `frobenius.cpp` | boundary at a regular point from the corner-1 (or CC1) regularised constants, log-power series along a ray, letter series by exact power-series algebra; any precision |
| `export_data.py` | `data/atilde.txt`, `data/letters.txt`, ray constants (`export_ray_constants`, also for CC1), exact w<=2 boundaries |
| `qdpatch/` | boost 1.83 `bulirsch_stoer.hpp` with explicit `double` casts (ambiguous size_t -> dd_real otherwise); algorithm unchanged |
| `camp/` | campaign scripts: `campaign.sh` (runs), `refs.py` (exact w<=2 references, also continued), `analyze.py` (tables + plots) |
| `run_cpp.py` | python driver |

## Build (cloud / Linux; needs boost >= 1.83 headers, libqd)

    python3 gen_code.py                  # gen/ (from ../../ccpipe/inputs/atilde.m)
    make -f build.mk -j2 h_family_c_double h_family_c_dd frobenius_dd

Compile times (2 cores): `h_family_c_*` 8 s (matrix as static data); unrolled variant 150-160 s per precision.
**dd_real needs `-ffp-contract=off`** (FMA contraction breaks the double-double error-free transformations;
the extrapolation then never converges).

## Run

    # boundary (dd) at P0 = s0 (X, Y) from corner-1 constants; RSIGN -1 for the continued region (r = -Sqrt|xyz|)
    python3 -c "import export_data as E; E.export_ray_constants('1','0.7')"
    ./frobenius_dd data data/raycst_1_0.7.txt 1 0.7 0.1 60 > data/bdy_w6_dd_ray1.txt
    # continued region: CC1 constants (transport.boundaries/H/dist/bdy_CC1_w*.m)
    python3 -c "import export_data as E; E.export_ray_constants('-1','-0.7', bdy_dir='.../transport.boundaries/H/dist', prefix='bdy_CC1', tag='CC1')"
    ./frobenius_dd data data/raycstCC1_-1_-0.7.txt -1 -0.7 0.1 60 -1 > data/bdy_w6_dd_CC1ray1.txt
    # evaluate (points "x y" on stdin)
    echo "0.3 0.2" | ./h_family_c_dd x data/bdy_w6_dd_ray1.txt 1e-14 0.1 0.2

%RESULTS%
