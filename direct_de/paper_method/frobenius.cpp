// frobenius.cpp -- boundary values at a non-singular point P0 = s0 (X, Y) from the corner-1
// (x, y -> 0) shuffle-regularised constants, at working precision (double / dd_real / qd_real).
//
// Along the ray x = s X, y = s Y the DE is dJ/ds = eps (R/s + sum_i B_i s^i) J with
// R = sum_k m_k a_k (m_k = order of vanishing of l_k at s = 0) and B_i = sum_k B_k[i] a_k.
// J_w(s) = sum_k log^k(s) sum_i P[w][k][i] s^i, built weight by weight; the s-regularised
// constants (log s -> 0) are the ray constants c^ray_w (data/raycst_X_Y.txt, export_data.py).
// The letter Taylor coefficients B_k[i] are computed EXACTLY by power-series algebra
// (no FFT): plain letters q'/q, odd letters i (2 b a' - a b') / (r (a^2 + b)) with
// r = s Sqrt[X Y] Sqrt[1 - s (X + Y)].
//
//   ./frobenius_<prec> DATA RAYCST X Y S0 N [RSIGN] > boundary   (paper's boundary-stream layout, order 6)
// RSIGN = +1 (default): r = +Sqrt[xyz] on the ray;  -1: r = -Sqrt|xyz| (continued region x, y < 0).
#include "h_data.hpp"
#include <iostream>

#ifndef REAL
#define REAL double
#endif
using Real = REAL;
using cplx = std::complex<Real>;
using H = HData<Real>;
using Ser = std::vector<cplx>;

static Ser smul(const Ser& a, const Ser& b, size_t L) {
    Ser c(L, cplx(0));
    for (size_t i = 0; i < std::min(L, a.size()); ++i)
        for (size_t j = 0; j < b.size() && i + j < L; ++j) c[i + j] += a[i] * b[j];
    return c;
}
static Ser sdiv(const Ser& a, const Ser& b, size_t L) {     // b[0] != 0
    Ser c(L, cplx(0));
    for (size_t i = 0; i < L; ++i) {
        cplx s = i < a.size() ? a[i] : cplx(0);
        for (size_t j = 1; j <= i && j < b.size(); ++j) s -= b[j] * c[i - j];
        c[i] = s / b[0];
    }
    return c;
}
static Ser sder(const Ser& a) { Ser c(a.size(), cplx(0)); for (size_t i = 1; i < a.size(); ++i) c[i - 1] = Real(double(i)) * a[i]; return c; }
static Ser sadd(Ser a, const Ser& b, cplx f = cplx(1)) { if (b.size() > a.size()) a.resize(b.size(), cplx(0)); for (size_t i = 0; i < b.size(); ++i) a[i] += f * b[i]; return a; }
static size_t val(const Ser& a, const Real& tiny) { size_t v = 0; while (v < a.size() && std::abs(a[v]) <= tiny) ++v; return v; }
static Ser shift(const Ser& a, size_t v) { return Ser(a.begin() + std::min(v, a.size()), a.end()); }

int main(int argc, char** argv) {
#ifdef DIFFEQS_QD
    unsigned int old_cw; fpu_fix_start(&old_cw);
#endif
    if (argc < 7) { std::cerr << "usage: DATA RAYCST X Y S0 N\n"; return 1; }
    H h(argv[1]);
    const Real X = from_string<Real>(argv[3]), Y = from_string<Real>(argv[4]), s0 = from_string<Real>(argv[5]);
    const int N = std::atoi(argv[6]);
    const Real rsign = (argc > 7) ? Real(std::atoi(argv[7])) : Real(1);
    const int n = h.n, nl = (int)h.letter.size();
    const Real tiny = Real(1e-28);
    const size_t L = N + 12;
    const cplx I(0, 1);

    // ---- letter series along the ray
    std::vector<int> m(nl, 0);
    std::vector<Ser> B(nl);
    auto ray = [&](const H::Poly& p) {
        Ser c(8, cplx(0));
        for (size_t t = 0; t < p.c.size(); ++t) {
            Real v = p.c[t]; for (int q = 0; q < p.i[t]; ++q) v *= X; for (int q = 0; q < p.j[t]; ++q) v *= Y;
            c[p.i[t] + p.j[t]] += v;
        }
        return c;
    };
    const Real S = X + Y;
    Ser x_ = {cplx(0), cplx(X)}, y_ = {cplx(0), cplx(Y)}, z_ = {cplx(1), cplx(-S)};
    Ser bpoly = smul(smul(x_, y_, L), z_, L);
    // Sqrt[1 - S s] = sum binom(1/2, k) (-S s)^k
    Ser g(L, cplx(0)); { Real cb = 1, ps = 1; for (size_t k = 0; k < L; ++k) { g[k] = cplx(cb * ps); cb = cb * (Real(0.5) - Real(double(k))) / Real(double(k + 1)); ps *= -S; } }
    for (int q = 0; q < nl; ++q) {
        int k = h.letter[q];
        if (!H::is_odd(k)) {
            Ser c = ray(h.plain.at(k));
            size_t v = val(c, tiny); m[q] = (int)v;
            Ser qq = shift(c, v); qq.resize(L, cplx(0));
            B[q] = sdiv(sder(qq), qq, L);
        } else {
            Ser a = (k == 7) ? smul(x_, z_, L) : smul(x_, y_, L);
            Ser num = sadd(smul(smul(bpoly, sder(a), L), Ser{cplx(2)}, L), smul(a, sder(bpoly), L), cplx(-1));
            Ser den = smul(smul(Ser{cplx(0), cplx(rsign * sqrt(X * Y))}, g, L), sadd(smul(a, a, L), bpoly), L);
            size_t vd = val(den, tiny), vn = val(num, tiny);
            if (vn < vd) { std::cerr << "odd letter has a pole on the ray\n"; return 2; }
            Ser nn = shift(num, vd), dd = shift(den, vd);
            B[q] = sdiv(nn, dd, L - vd - 1);
            for (auto& v : B[q]) v *= I;
            m[q] = 0;
        }
        B[q].resize(N + 1, cplx(0));
    }

    // ---- ray constants
    std::ifstream cf(argv[2]); int wmax; cf >> wmax;
    std::vector<std::vector<cplx>> cr(wmax + 1, std::vector<cplx>(n));
    for (auto& v : cr) for (auto& e : v) { std::string a, b; cf >> a >> b; e = cplx(from_string<Real>(a), from_string<Real>(b)); }
    if (!cf) { std::cerr << "bad raycst\n"; return 1; }

    // ---- recursion: P[w][k][i][n]
    using Mat = std::vector<std::vector<cplx>>;               // [i][component]
    auto zeros = [&]() { return Mat(N + 1, std::vector<cplx>(n, cplx(0))); };
    std::vector<std::map<int, Mat>> P(wmax + 1);
    P[0][0] = zeros(); P[0][0][0] = cr[0];
    std::vector<cplx> conv(n), tmp(n);
    for (int w = 1; w <= wmax; ++w) {
        std::map<int, Mat> G;
        for (auto& [k, Pk] : P[w - 1]) {
            Mat g = zeros();
            for (int i = 0; i <= N; ++i) {
                for (int q = 0; q < nl; ++q) {
                    // coefficient vector multiplying a_q at s^(i-1):  m_q P_k[i] + sum_j B_q[j] P_k[i-1-j]
                    std::fill(conv.begin(), conv.end(), cplx(0));
                    bool any = false;
                    if (m[q]) { for (int c = 0; c < n; ++c) conv[c] += Real(double(m[q])) * Pk[i][c]; any = true; }
                    for (int j = 0; j < i; ++j) {
                        const cplx bq = B[q][j];
                        if (bq == cplx(0)) continue;
                        const auto& Pv = Pk[i - 1 - j];
                        for (int c = 0; c < n; ++c) conv[c] += bq * Pv[c];
                        any = true;
                    }
                    if (any) h.apply(q, conv.data(), g[i].data());
                }
            }
            G[k] = std::move(g);
        }
        std::map<int, Mat> nw;
        auto slot = [&](int k) -> Mat& { auto it = nw.find(k); if (it == nw.end()) it = nw.emplace(k, zeros()).first; return it->second; };
        for (auto& [k, g] : G) {
            for (int i = 0; i <= N; ++i) {
                if (i == 0) { Mat& t = slot(k + 1); for (int c = 0; c < n; ++c) t[0][c] += g[0][c] / Real(double(k + 1)); continue; }
                Real fk = 1;                                       // k!/(k-j)!
                for (int j = 0; j <= k; ++j) {
                    if (j) fk *= Real(double(k - j + 1));
                    Real coef = fk; for (int e = 0; e <= j; ++e) coef /= Real(double(i));
                    if (j % 2) coef = -coef;
                    Mat& t = slot(k - j);
                    for (int c = 0; c < n; ++c) t[i][c] += coef * g[i][c];
                }
            }
        }
        Mat& t0 = slot(0); for (int c = 0; c < n; ++c) t0[0][c] += cr[w][c];
        P[w] = std::move(nw);
        std::cerr << "weight " << w << " done\n";
    }

    // ---- evaluate at s0 and write the boundary stream
    const Real ls = log(s0);
    const Real x0 = s0 * X, y0 = s0 * Y;
    std::cout << "2\n" << to_string<Real>(x0) << " " << to_string<Real>(y0) << "\n" << wmax << " " << n << " 1\n";
    for (int w = 0; w <= wmax; ++w) {
        std::vector<cplx> Jc(n, cplx(0));
        for (auto& [k, Pk] : P[w]) {
            Real lk = 1; for (int e = 0; e < k; ++e) lk *= ls;
            std::vector<cplx> acc(n, cplx(0));
            for (int i = N; i >= 0; --i) for (int c = 0; c < n; ++c) acc[c] = acc[c] * s0 + Pk[i][c];
            for (int c = 0; c < n; ++c) Jc[c] += lk * acc[c];
        }
        for (auto& v : Jc) std::cout << to_string<Real>(v.real()) << " " << to_string<Real>(v.imag()) << "\n";
    }
    const Real r0 = rsign * sqrt(x0 * y0 * (Real(1) - x0 - y0));
    std::cout << to_string<Real>(r0) << " " << to_string<Real>(Real(0)) << "\n";
#ifdef DIFFEQS_QD
    fpu_fix_end(&old_cw);
#endif
    return 0;
}
