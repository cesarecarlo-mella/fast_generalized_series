// h_data.hpp -- H family data: letters and connection matrices, read from data/ (export_data.py).
//   plain letters  : polynomials in x, y                       (letters.txt)
//   odd letters    : l7, l8 = (a - i r)/(a + i r), r^2 = b = x y z, a = x z (l7), a = x y (l8)
//   connection     : A = sum_k a_k dlog l_k, a_k sparse 371 x 371   (atilde.txt, 40 digits)
#pragma once
#include "diffeqs.hpp"
#include <fstream>
#include <map>
#include <tuple>

template<class Real>
struct HData {
    using cplx = std::complex<Real>;
    struct Poly { std::vector<int> i, j; std::vector<Real> c; };
    int n = 0;
    std::map<int, Poly> plain;
    std::vector<int> letter;                                        // letter ids, in file order
    std::vector<std::vector<std::tuple<int, int, cplx>>> ent;       // per letter: (row, col, value)
    // union sparsity pattern (CSR) and, per letter, positions of its entries in it
    std::vector<size_t> rowptr, col;
    std::vector<std::vector<std::pair<size_t, cplx>>> lpos;
    // same split into purely real entries (almost all) and the rest, for cheaper assembly
    std::vector<std::vector<std::pair<size_t, Real>>> lpos_re;
    std::vector<std::vector<std::pair<size_t, cplx>>> lpos_cx;

    static bool is_odd(int k) { return k == 7 || k == 8; }

    explicit HData(const std::string& dir) {
        {
            std::ifstream f(dir + "/letters.txt"); int np; f >> np;
            for (int q = 0; q < np; ++q) {
                int k, nt; f >> k >> nt; Poly p;
                for (int t = 0; t < nt; ++t) { int a, b; std::string c; f >> a >> b >> c; p.i.push_back(a); p.j.push_back(b); p.c.push_back(from_string<Real>(c)); }
                plain[k] = p;
            }
            if (!f) throw std::runtime_error("letters.txt");
        }
        std::ifstream f(dir + "/atilde.txt"); int nl; f >> n >> nl;
        std::map<std::pair<int, int>, size_t> uni;
        ent.resize(nl);
        for (int q = 0; q < nl; ++q) {
            int k, nnz; f >> k >> nnz; letter.push_back(k);
            for (int e = 0; e < nnz; ++e) {
                int i, j; std::string re, im; f >> i >> j >> re >> im;
                ent[q].emplace_back(i, j, cplx(from_string<Real>(re), from_string<Real>(im)));
                uni[{i, j}] = 0;
            }
        }
        if (!f) throw std::runtime_error("atilde.txt");
        rowptr.assign(n + 1, 0);
        size_t pos = 0;
        for (auto& kv : uni) { kv.second = pos++; col.push_back(kv.first.second); rowptr[kv.first.first + 1]++; }
        for (int i = 0; i < n; ++i) rowptr[i + 1] += rowptr[i];
        lpos.resize(nl);
        for (int q = 0; q < nl; ++q)
            for (auto& [i, j, v] : ent[q]) lpos[q].emplace_back(uni[{i, j}], v);
        lpos_re.resize(nl); lpos_cx.resize(nl);
        for (int q = 0; q < nl; ++q)
            for (auto& [p, v] : lpos[q]) {
                if (v.imag() == Real(0)) lpos_re[q].emplace_back(p, v.real());
                else lpos_cx[q].emplace_back(p, v);
            }
    }
    // out = a_q * in   (single letter)
    void apply(int q, const cplx* in, cplx* out) const {
        for (auto& [i, j, v] : ent[q]) out[i] += v * in[j];
    }
};
