// diffeqs.hpp -- solver for canonical differential equations dJ = eps A(x) J, following
// M. Czakon, L. Tancredi, "Solution of Canonical Differential Equations for Integrals on
// Arbitrary Geometries", arXiv:2606.30354, Sec. 2 and 4:
//
//  * ONE first-order system containing all eps orders f_{i,1..order} (f_{i,0} is constant and
//    is not integrated) PLUS the special functions (square roots, ...) with their own
//    differential equations  df_k = V_k(x, f)  (eq. 2.5/2.6: d sqrt(g) = sqrt(g) dg / (2 g)).
//  * integrated with boost::numeric::odeint::bulirsch_stoer along the complexified path
//      z_k(t) = x0_k + (1 + 4 i delta_k (1 - t)) (x_k(t) - x0_k),  x(t) = x0 + t (x1 - x0)   (eq. 4.1)
//  * default local-error estimate: maximal ABSOLUTE error of the functions
//    (bulirsch_stoer with eps_abs = error, eps_rel = 0).
//  * Real = double or dd_real (QD library).
//
// Public interface as in Sec. 4 (diffeqs_types, connection_t, vector_field_t, path_t,
// error_estimator_t, evaluate, change_path, set_log_stream, set_max_steps/evals/time, set_min_step).
// Deviations (documented): matrix is a CSR struct with a fixed sparsity pattern instead of a
// ublas::compressed_matrix (same content, avoids re-allocation per call); change_error_estimator
// is accepted but boost's bulirsch_stoer hard-wires its error checker, so only the default
// (max absolute error) is used.
#pragma once
#include <complex>
#include <vector>
#include <functional>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <chrono>
#include <boost/numeric/odeint.hpp>

#ifdef DIFFEQS_QD
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#endif

// ---------------------------------------------------------------- odeint support for complex states
namespace boost { namespace numeric { namespace odeint {
template<class R>
struct norm_result_type< std::vector< std::complex<R> > > { typedef R type; };
} } }

template<class Real> inline Real from_string(const std::string& s);
template<> inline double from_string<double>(const std::string& s) { return std::stod(s); }
template<class Real> inline std::string to_string(const Real& v);
template<> inline std::string to_string<double>(const double& v) {
    char b[64]; snprintf(b, sizeof b, "%.17e", v); return b; }
#ifdef DIFFEQS_QD
template<> inline dd_real from_string<dd_real>(const std::string& s) { return dd_real(s.c_str()); }
template<> inline std::string to_string<dd_real>(const dd_real& v) { return v.to_string(34, 0, std::ios_base::scientific); }
template<> inline qd_real from_string<qd_real>(const std::string& s) { return qd_real(s.c_str()); }
template<> inline std::string to_string<qd_real>(const qd_real& v) { return v.to_string(66, 0, std::ios_base::scientific); }
#endif

// ---------------------------------------------------------------- types (Table 2)
template<class Real>
struct diffeqs_types {
    using real = Real;
    using cplx = std::complex<Real>;
    using real_vec = std::vector<Real>;
    using cplx_vec = std::vector<cplx>;
    struct cplx_vec_range {                 // view into a complex vector
        cplx* p; size_t n;
        cplx& operator[](size_t i) const { return p[i]; }
        size_t size() const { return n; }
    };
    struct matrix {                         // CSR, fixed pattern, values filled by connection_t
        size_t n = 0;
        std::vector<size_t> rowptr, col;
        std::vector<cplx> val;
    };
    using connection_t = std::function<void(const cplx_vec& x, const cplx_vec& dx,
                                            const cplx_vec_range& f, matrix& A)>;
    using vector_field_t = std::function<void(const cplx_vec& x, const cplx_vec& dx,
                                              const cplx_vec_range& f, cplx_vec_range& df)>;
    using path_t = std::function<void(const real_vec& x0, const real_vec& x1, real t,
                                      real_vec& x, real_vec& dx)>;
    using error_estimator_t = std::function<real(const cplx_vec& f, const cplx_vec& f_err)>;
    // optional: compiled product out = A in (A filled by connection_t); default: generic CSR loop
    using matvec_t = std::function<void(const matrix& A, const cplx* in, cplx* out)>;
};

// ---------------------------------------------------------------- exceptions (Table 3)
struct diffeqs_exception : std::runtime_error { using std::runtime_error::runtime_error; };
struct diffeqs_max_steps_exception : diffeqs_exception { diffeqs_max_steps_exception() : diffeqs_exception("max steps") {} };
struct diffeqs_max_evals_exception : diffeqs_exception { diffeqs_max_evals_exception() : diffeqs_exception("max evals") {} };
struct diffeqs_max_time_exception  : diffeqs_exception { diffeqs_max_time_exception()  : diffeqs_exception("max time") {} };
struct diffeqs_min_step_exception  : diffeqs_exception { diffeqs_min_step_exception()  : diffeqs_exception("min step") {} };

// ---------------------------------------------------------------- solver
template<class Real = double>
class diffeqs {
public:
    using T = diffeqs_types<Real>;
    using real = typename T::real;  using cplx = typename T::cplx;
    using real_vec = typename T::real_vec;  using cplx_vec = typename T::cplx_vec;
    using cplx_vec_range = typename T::cplx_vec_range;  using matrix = typename T::matrix;
    using connection_t = typename T::connection_t;  using vector_field_t = typename T::vector_field_t;
    using path_t = typename T::path_t;  using error_estimator_t = typename T::error_estimator_t;
    using matvec_t = typename T::matvec_t;

    // A must come with its sparsity pattern (rowptr, col) set; connection_t fills A.val.
    diffeqs(connection_t connection, vector_field_t vector_field, std::istream& boundary, matrix pattern)
        : conn_(connection), vf_(vector_field), A_(pattern)
    {
        std::string s;
        boundary >> dimx_;
        x0_.resize(dimx_);
        for (auto& v : x0_) { boundary >> s; v = from_string<Real>(s); }
        boundary >> order_ >> dimbasis_ >> nspecial_;
        f0_.resize((order_ + 1) * dimbasis_);
        for (auto& v : f0_) { std::string a, b; boundary >> a >> b; v = cplx(from_string<Real>(a), from_string<Real>(b)); }
        sp0_.resize(nspecial_);
        for (auto& v : sp0_) { std::string a, b; boundary >> a >> b; v = cplx(from_string<Real>(a), from_string<Real>(b)); }
        if (!boundary) throw diffeqs_exception("bad boundary stream");
        path_ = [](const real_vec& a, const real_vec& b, real t, real_vec& x, real_vec& dx) {
            x.resize(a.size()); dx.resize(a.size());
            for (size_t k = 0; k < a.size(); ++k) { x[k] = a[k] + t * (b[k] - a[k]); dx[k] = b[k] - a[k]; }
        };
    }

    // integrate from the boundary point to x; result f = (f_.,1, ..., f_.,order, special)
    void evaluate(const real_vec& x, const real_vec& deformation, const real& error, cplx_vec& f)
    {
        x1_ = x; delta_ = deformation;
        const size_t N = order_ * dimbasis_ + nspecial_;
        cplx_vec state(N);
        for (size_t w = 1; w <= order_; ++w)
            for (size_t i = 0; i < dimbasis_; ++i) state[(w - 1) * dimbasis_ + i] = f0_[w * dimbasis_ + i];
        for (size_t k = 0; k < nspecial_; ++k) state[order_ * dimbasis_ + k] = sp0_[k];

        using namespace boost::numeric::odeint;
        // default (paper): max ABSOLUTE error.  set_relative_error(e): boost's error checker with
        // eps_rel = e, a_x = 1, i.e. |err_i| <= error + e |f_i|  (recommendation 1 of the paper: tailor the estimate)
        bulirsch_stoer<cplx_vec, real, cplx_vec, real> stepper(error, rel_, rel_ == real(0) ? real(0) : real(1), real(0));
        steps_ = 0; evals_ = 0; rejected_ = 0;
        auto t_start = std::chrono::steady_clock::now();
        auto sys = [this](const cplx_vec& s, cplx_vec& ds, const real t) { this->rhs(s, ds, t); };

        real t = 0, dt = real(1) / 10;
        const real tend = 1;
        while (t < tend) {
            if (t + dt > tend) dt = tend - t;
            controlled_step_result res;
            do {
                res = stepper.try_step(sys, state, t, dt);
                if (res == fail) ++rejected_;
                check_limits(dt, t_start);
            } while (res == fail);
            ++steps_;
            if (log_) *log_ << "t=" << to_string<Real>(t) << " dt=" << to_string<Real>(dt) << " evals=" << evals_ << "\n";
        }
        f = state;
        time_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    }

    void change_error_estimator(error_estimator_t e) { (void)e; /* not supported by boost bulirsch_stoer, see header */ }
    void change_path(path_t p) { path_ = p; }
    void set_matvec(matvec_t m) { mv_ = m; }
    void set_log_stream(std::ostream& out) { log_ = &out; }
    void set_max_steps(unsigned n) { max_steps_ = n; }
    void set_max_evals(unsigned n) { max_evals_ = n; }
    void set_max_time(unsigned n) { max_time_ = n; }
    void set_min_step(real d) { min_step_ = d; }
    void set_relative_error(real e) { rel_ = e; }

    // diagnostics
    size_t steps() const { return steps_; }
    size_t evals() const { return evals_; }
    size_t rejected() const { return rejected_; }
    double time() const { return time_; }
    size_t order() const { return order_; }
    size_t dimbasis() const { return dimbasis_; }
    const real_vec& x0() const { return x0_; }

private:
    void check_limits(const real& dt, std::chrono::steady_clock::time_point t0) {
        if (max_steps_ && steps_ > max_steps_) throw diffeqs_max_steps_exception();
        if (max_evals_ && evals_ > max_evals_) throw diffeqs_max_evals_exception();
        if (max_time_ && std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > max_time_)
            throw diffeqs_max_time_exception();
        if (min_step_ != real(0) && dt < min_step_) throw diffeqs_min_step_exception();
    }

    // complexified point and tangent, eq. (4.1)
    void point(real t, cplx_vec& z, cplx_vec& dz) {
        real_vec x, dx;
        path_(x0_, x1_, t, x, dx);
        z.resize(dimx_); dz.resize(dimx_);
        const cplx I(0, 1);
        for (size_t k = 0; k < dimx_; ++k) {
            cplx g = cplx(1) + real(4) * delta_[k] * (real(1) - t) * I;     // 1 + 4 i delta (1 - t)
            cplx dg = -real(4) * delta_[k] * I;
            z[k] = cplx(x0_[k]) + g * cplx(x[k] - x0_[k]);
            dz[k] = dg * cplx(x[k] - x0_[k]) + g * cplx(dx[k]);
        }
    }

    void rhs(const cplx_vec& s, cplx_vec& ds, const real t) {
        ++evals_;
        cplx_vec z, dz;
        point(t, z, dz);
        cplx_vec_range sp{const_cast<cplx*>(s.data()) + order_ * dimbasis_, nspecial_};
        conn_(z, dz, sp, A_);                                  // A = A_t (connection contracted with dz/dt)
        // d f_w / dt = A f_{w-1}
        for (size_t w = 1; w <= order_; ++w) {
            const cplx* in = (w == 1) ? f0_.data() : s.data() + (w - 2) * dimbasis_;
            cplx* out = ds.data() + (w - 1) * dimbasis_;
            if (mv_) { mv_(A_, in, out); continue; }
            for (size_t i = 0; i < dimbasis_; ++i) {
                cplx acc(0);
                for (size_t p = A_.rowptr[i]; p < A_.rowptr[i + 1]; ++p) acc += A_.val[p] * in[A_.col[p]];
                out[i] = acc;
            }
        }
        cplx_vec_range dsp{ds.data() + order_ * dimbasis_, nspecial_};
        vf_(z, dz, sp, dsp);
    }

    connection_t conn_; vector_field_t vf_; path_t path_; matrix A_; matvec_t mv_;
    size_t dimx_ = 0, order_ = 0, dimbasis_ = 0, nspecial_ = 0;
    real_vec x0_, x1_, delta_;
    cplx_vec f0_, sp0_;
    std::ostream* log_ = nullptr;
    unsigned max_steps_ = 0, max_evals_ = 0, max_time_ = 0;
    real min_step_ = 0, rel_ = 0;
    size_t steps_ = 0, evals_ = 0, rejected_ = 0; double time_ = 0;
};
