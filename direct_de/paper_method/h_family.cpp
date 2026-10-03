// h_family.cpp -- H family (371 canonical masters, 14 letters, one square root r = Sqrt[x y z],
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
#include "h_data.hpp"
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
using H = HData<Real>;

static inline cplx ipow(const cplx& x, int k) { cplx r(1); for (int q = 0; q < k; ++q) r *= x; return r; }

static cplx peval(const H::Poly& p, const cplx& x, const cplx& y) {
    cplx s(0);
    for (size_t n = 0; n < p.c.size(); ++n) s += p.c[n] * ipow(x, p.i[n]) * ipow(y, p.j[n]);
    return s;
}
static cplx pdt(const H::Poly& p, const cplx& x, const cplx& y, const cplx& dx, const cplx& dy) {
    cplx s(0);
    for (size_t n = 0; n < p.c.size(); ++n) {
        int a = p.i[n], b = p.j[n];
        if (a) s += p.c[n] * Real(a) * ipow(x, a - 1) * ipow(y, b) * dx;
        if (b) s += p.c[n] * Real(b) * ipow(x, a) * ipow(y, b - 1) * dy;
    }
    return s;
}

int main(int argc, char** argv) {
#ifdef DIFFEQS_QD
    unsigned int old_cw; fpu_fix_start(&old_cw);
#endif
    if (argc < 6) { std::cerr << "usage: DATA BOUNDARY ERROR DELTAX DELTAY [--quiet] < points\n"; return 1; }
    const bool quiet = (argc > 6 && !std::strcmp(argv[6], "--quiet"));
    H h(argv[1]);
    Real error = from_string<Real>(argv[3]);
    D::real_vec delta = {from_string<Real>(argv[4]), from_string<Real>(argv[5])};
    const cplx I(0, 1);
    const int nl = (int)h.letter.size();

    D::matrix pattern;
    pattern.n = h.n; pattern.rowptr = h.rowptr; pattern.col = h.col; pattern.val.assign(h.col.size(), cplx(0));

    // connection: A_t = sum_k a_k d/dt log l_k(z(t)); odd letters use the special function r = f[0]
    D::connection_t connection = [&](const D::cplx_vec& z, const D::cplx_vec& dz, const D::cplx_vec_range& f, D::matrix& A) {
        const cplx x = z[0], y = z[1], dx = dz[0], dy = dz[1];
        const cplx zz = cplx(1) - x - y, dzz = -dx - dy;
        const cplx r = f[0];
        const cplx b = x * y * zz, db = dx * y * zz + x * dy * zz + x * y * dzz;
        std::fill(A.val.begin(), A.val.end(), cplx(0));
        for (int q = 0; q < nl; ++q) {
            const int k = h.letter[q];
            cplx c;
            if (H::is_odd(k)) {                       // d log (a - i r)/(a + i r) with r^2 = b
                cplx a, da;
                if (k == 7) { a = x * zz; da = dx * zz + x * dzz; }
                else        { a = x * y;  da = dx * y + x * dy; }
                c = I * (Real(2) * b * da - a * db) / (r * (a * a + b));
            } else {
                const H::Poly& p = h.plain.at(k);
                c = pdt(p, x, y, dx, dy) / peval(p, x, y);
            }
            for (auto& [pos, v] : h.lpos_re[q]) A.val[pos] += c * v;      // complex * real
            for (auto& [pos, v] : h.lpos_cx[q]) A.val[pos] += c * v;
        }
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
