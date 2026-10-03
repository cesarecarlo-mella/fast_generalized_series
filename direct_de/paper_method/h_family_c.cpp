// h_family_c.cpp -- connection matrix COMPILED IN as static data (gen/hgen_data.cpp, exact hi/lo doubles) + letters as generated code.  H family (371 canonical masters, 14 letters, one square root r = Sqrt[x y z],
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
#include "gen/hgen_data.hpp"
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

// integration variables.  Default: (x, y).  -DUV_CHART: scattering-region chart (u, v),
//   x = -w/v, y = -u/v, z = 1/v, w = 1 - u - v  (the path is straight in (u, v)); the DE is the
//   same canonical DE, the letters are composed with the chart (chain rule on the tangent).
static inline void to_xy(const D::cplx_vec& z, const D::cplx_vec& dz, cplx& x, cplx& y, cplx& dx, cplx& dy) {
#ifdef UV_CHART
    const cplx u = z[0], v = z[1], du = dz[0], dv = dz[1];
    const cplx w = cplx(1) - u - v, dw = -du - dv;
    x = -w / v;  dx = -(dw * v - w * dv) / (v * v);
    y = -u / v;  dy = -(du * v - u * dv) / (v * v);
#else
    x = z[0]; y = z[1]; dx = dz[0]; dy = dz[1];
#endif
}

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
        cplx x, y, dx, dy; to_xy(z, dz, x, y, dx, dy);
        hgen_dlogs<Real>(x, y, dx, dy, f[0], c);
        hgen::fill<Real>(c, A.val.data());
    };
    // special function: dr = r db / (2 b)   (eq. 2.6)
    D::vector_field_t vector_field = [&](const D::cplx_vec& z, const D::cplx_vec& dz, const D::cplx_vec_range& f, D::cplx_vec_range& df) {
        cplx x, y, dx, dy; to_xy(z, dz, x, y, dx, dy);
        const cplx zz = cplx(1) - x - y, dzz = -dx - dy;
        const cplx b = x * y * zz, db = dx * y * zz + x * dy * zz + x * y * dzz;
        df[0] = f[0] * db / (Real(2) * b);
    };

    std::ifstream bf(argv[2]);
    D solver(connection, vector_field, bf, pattern);
    solver.set_matvec([](const D::matrix& A, const cplx* in, cplx* out) { hgen::mv<Real>(A.val.data(), in, out); });
    solver.set_max_evals(5000000);
    if (const char* mt = std::getenv("DIFFEQS_MAX_TIME")) solver.set_max_time(std::atoi(mt));
    if (const char* er = std::getenv("DIFFEQS_EPS_REL")) solver.set_relative_error(from_string<Real>(er));

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
