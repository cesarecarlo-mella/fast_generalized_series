// h_family_gen.cpp -- COMPILED connection (gen/, from gen_code.py).  H family (371 canonical masters, 14 letters, one square root r = Sqrt[x y z],
// z = 1 - x - y) solved with diffeqs<Real> (method of arXiv:2606.30354).
//
//   ./h_family_<prec> DATA BOUNDARY ERROR DELTAX DELTAY [--quiet] < points > results
// points : one "x y" per line (decimal strings, read at full working precision)
// results: per point "# x y steps evals rejected time_s" then order*371 + 1 lines "re im"
//          (f_{.,1} ... f_{.,order}, r);  with --quiet only the "#" line.
// env DIFFEQS_MAX_TIME=<s> : per-point time limit (diffeqs_max_time_exception -> "FAILED").
//
// The square root is NOT evaluated: r is a component of the system, dr/dt = r (db/dt) / (2 b),
// started from the boundary value in BOUNDARY (its branch there is the only input; e.g. +Sqrt in
// the Euclidean region, -Sqrt|xyz| in the continued region).
#include "diffeqs.hpp"
#include "gen/hgen.hpp"
#include "gen/hgen_letters.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>

#ifndef REAL
#define REAL double
#endif
using Real = REAL;
using D = diffeqs<Real>;
using cplx = D::cplx;
int main(int argc, char** argv) {
#ifdef DIFFEQS_QD
    unsigned int old_cw; fpu_fix_start(&old_cw);
#endif
    if (argc < 6) { std::cerr << "usage: DATA BOUNDARY ERROR DELTAX DELTAY [--quiet] < points\n"; return 1; }
    const bool quiet = (argc > 6 && !std::strcmp(argv[6], "--quiet"));
    // argv[1] (DATA) is ignored: everything is compiled in
    Real error = from_string<Real>(argv[3]);
    D::real_vec delta = {from_string<Real>(argv[4]), from_string<Real>(argv[5])};
    D::matrix pattern;
    pattern.n = hgen::N; pattern.val.assign(hgen::NNZ, cplx(0));
    D::connection_t connection = [&](const D::cplx_vec& z, const D::cplx_vec& dz, const D::cplx_vec_range& f, D::matrix& A) {
        cplx c[hgen::NL];
        hgen_dlogs<Real>(z[0], z[1], dz[0], dz[1], f[0], c);
        hgen_fill<Real>(c, A.val.data());
    };
    // special function: dr = r db / (2 b)   (eq. 2.6)
    D::vector_field_t vector_field = [&](const D::cplx_vec& z, const D::cplx_vec& dz, const D::cplx_vec_range& f, D::cplx_vec_range& df) {
        const cplx x = z[0], y = z[1], dx = dz[0], dy = dz[1];
        const cplx zz = cplx(1) - x - y, dzz = -dx - dy;
        const cplx b = x * y * zz, db = dx * y * zz + x * dy * zz + x * y * dzz;
        df[0] = f[0] * db / (Real(2) * b);
    };

    std::ifstream bf(argv[2]);
    D solver(connection, vector_field, bf, pattern);
    solver.set_matvec([](const D::matrix& A, const cplx* in, cplx* out) { hgen_mv<Real>(A.val.data(), in, out); });
    solver.set_max_evals(5000000);
    if (const char* mt = std::getenv("DIFFEQS_MAX_TIME")) solver.set_max_time(std::atoi(mt));

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line); std::string sx, sy;
        if (!(is >> sx >> sy)) continue;
        D::real_vec x1 = {from_string<Real>(sx), from_string<Real>(sy)};
        D::cplx_vec f;
        try {
            solver.evaluate(x1, delta, error, f);
        } catch (diffeqs_exception& e) {
            std::cout << "# " << sx << " " << sy << " FAILED " << e.what() << std::endl; continue;
        }
        std::cout << "# " << sx << " " << sy << " " << solver.steps() << " " << solver.evals() << " "
                  << solver.rejected() << " " << solver.time() << "\n";
        if (!quiet)
            for (auto& v : f) std::cout << to_string<Real>(v.real()) << " " << to_string<Real>(v.imag()) << "\n";
        std::cout.flush();
    }
#ifdef DIFFEQS_QD
    fpu_fix_end(&old_cw);
#endif
    return 0;
}
