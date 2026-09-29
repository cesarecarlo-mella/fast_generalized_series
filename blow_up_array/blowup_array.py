"""
blowup_array.py -- the blow-up recursion done as ARRAY arithmetic on fixed
series "slots" instead of symbolic expressions (no Mathematica, no FORM).

Same mathematics as ../blow_up_form*/pipeline/01_recur_step*.wl:
    multT  = (atilde /. dt) . Jprev          multY = (atilde /. dy) . Jprev
    intt   = Int[multT, t]                    tointy = -D[intt, BIG] + multY
    toaddy = Int[tointy, BIG]                 J_w    = intt + toaddy + bdy_w
with the same truncation (t^e dropped for e > ORDt, BIG^e for e > ORDy).

SLOTS.  Every quantity is a set of dense 2D coefficient grids, one per
"tag" = (parity_t, parity_BIG, logpow_t, logpow_BIG, constant monomial):
    - exponents are half-integers:  e = n + parity/2,  grid index n+1
      (n = -1 .. ORD, so the 1/t, 1/BIG poles of the dlogs fit),
    - Log[t]^a Log[BIG]^b and the transcendental constant monomial
      (Pi^k Zeta[3]^m ... I^0/1) never interact with the (t,BIG) grid under
      multiplication by a dlog, so they are just labels riding along.
The dlog series carry no logs and no transcendentals, so the product is a
pure 2D truncated convolution of numbers (the "precomputed slots" idea).

REORGANIZATION.  rowVec.Jprev = sum_l dlog_l * (sum_j atilde_ij,l J_j):
first a plain linear combination over components (no convolution), then
ONE convolution per (letter, row) -- not one per nonzero matrix entry.
Both steps are single dense matrix products (BLAS):
    V_l      = A_l  (ncomp x ncomp)  @  J      (ncomp x tags*slots)
    mult    += V_l  (tags*ncomp x S) @  T_l    (S x S, the dlog's shift matrix)

ARITHMETIC.  All of it is done modulo word-size primes p < 2^22 -- the
common-denominator idea pushed to its limit: every rational (atilde
entries, dlog coefficients with their 2^85 denominators, 1/(a+1)^(j+1)
integration factors, boundary constants) becomes ONE machine integer, and
every sum is a plain integer multiply-add with no gcd anywhere.  With
p < 2^22 all BLAS partial sums stay < 2^53, so float64 BLAS is EXACT.
Exact rationals are recovered once at the end (CRT over several primes +
rational reconstruction), and verified against one extra prime.
"""
import math
import os
import time
import numpy as np
import gmpy2
from gmpy2 import mpq, mpz

import mparse as mp

# --------------------------------------------------------------------------
PRIMES_BITS = 22


def gen_primes(n, bits=PRIMES_BITS):
    out, c = [], (1 << bits) - 1
    while len(out) < n:
        c -= 2
        if gmpy2.is_prime(c):
            out.append(int(c))
    return out


def rmod(q, p):
    """rational -> residue mod p."""
    return int(mpz(q.numerator) * gmpy2.invert(mpz(q.denominator), p) % p)


def key_mul(a, b):
    k, s = mp._mulkey(a, b)
    return k, s


I_KEY = (('I', 2),)


# ==========================================================================
#   loading / slot decomposition (exact, prime independent)
# ==========================================================================
class Corner:
    """everything weight-independent for one corner at one order."""

    def __init__(self, ddir, corner, ORDt, ORDy, verbose=True):
        self.c, self.ORDt, self.ORDy = corner, ORDt, ORDy
        self.tv = 't'
        self.yv = 'z' if corner == 3 else 'y'
        self.ltv, self.lyv = 'Log[t]', 'Log[%s]' % self.yv
        self.Nt, self.Ny = ORDt + 2, ORDy + 2          # n = -1..ORD
        self.S = self.Nt * self.Ny
        t0 = time.time()
        A = mp.parse_file(os.path.join(ddir, 'atilde.m'))
        self.ncomp = len(A)
        # atilde_ij = sum_l (re + i im) l  ->  per letter: {(i,j): (re, im)}
        self.A = {}
        for i, row in enumerate(A):
            for j, e in enumerate(row):
                for k, v in e.items():
                    d = dict(k)
                    im = d.pop('I', 0)
                    (lt, _), = d.items()
                    slot = self.A.setdefault(lt, {}).setdefault((i, j), [mpq(0), mpq(0)])
                    slot[1 if im else 0] += v
        self.letters = sorted(self.A, key=lambda s: int(s[2:-1]))
        self.kern = {}
        for dname, f in (('t', 'dlog_dt_c%d.m' % corner), ('y', 'dlog_dy_c%d.m' % corner)):
            rules = mp.parse_file(os.path.join(ddir, f))
            for _, lhs, rhs in rules:
                (lk, _), = next(iter(lhs))  # lhs is {((l[k],2),):1}
                self.kern[(dname, lk)] = self._kernel(rhs)
        if verbose:
            print("[load c%d] %d comps, letters %s, parsed in %.1fs" %
                  (corner, self.ncomp, self.letters, time.time() - t0), flush=True)

    # exponent (exp2 units) -> (grid index, parity) or None if truncated
    def slot(self, e2, ORD):
        n, par = e2 >> 1, e2 & 1
        if e2 > 2 * ORD:
            return None
        if n < -1:
            raise ValueError("exponent below grid: %s/2" % e2)
        return n + 1, par

    def _kernel(self, poly):
        """dlog series -> {(par_t, par_y, iflag): [(it, iy, coef mpq), ...]}"""
        out = {}
        for k, v in poly.items():
            d = dict(k)
            fi = 1 if d.pop('I', 0) else 0
            et = d.pop(self.tv, 0); ey = d.pop(self.yv, 0)
            if d:
                raise ValueError("unexpected atoms in dlog: %r" % d)
            st, sy = self.slot(et, self.ORDt), self.slot(ey, self.ORDy)
            if st is None or sy is None:
                continue
            out.setdefault((st[1], sy[1], fi), []).append((st[0], sy[0], v))
        return out

    def split_series(self, poly):
        """series Poly -> {tag: [(it, iy, coef)]}, tag=(pt,py,lt,ly,ckey)"""
        out = {}
        for k, v in poly.items():
            d = dict(k)
            et = d.pop(self.tv, 0); ey = d.pop(self.yv, 0)
            lt = d.pop(self.ltv, 0) // 2; ly = d.pop(self.lyv, 0) // 2
            st, sy = self.slot(et, self.ORDt), self.slot(ey, self.ORDy)
            if st is None or sy is None:
                continue
            ckey = tuple(sorted(d.items()))
            out.setdefault((st[1], sy[1], lt, ly, ckey), []).append((st[0], sy[0], v))
        return out

    def load_vector(self, path):
        """exact_*.m / sol0.m / bdy_*.m  ->  {(tag, comp): [(it,iy,coef)]}"""
        vec = mp.parse_file(path)
        assert len(vec) == self.ncomp
        res = {}
        for j, p in enumerate(vec):
            for tag, ents in self.split_series(p).items():
                res[(tag, j)] = ents
        return res


# ==========================================================================
#   modular engine (one prime)
# ==========================================================================
class ModEngine:
    def __init__(self, C, p):
        self.C, self.p = C, p
        n = C.ncomp
        # dense atilde matrices (balanced residues, float64): letter -> (Re, Im)
        self.Am = {}
        for lt in C.letters:
            R = np.zeros((n, n)); Im = np.zeros((n, n))
            for (i, j), (re, im) in C.A[lt].items():
                if re:
                    R[i, j] = self.bal(rmod(re, p))
                if im:
                    Im[i, j] = self.bal(rmod(im, p))
            self.Am[lt] = (R if R.any() else None, Im if Im.any() else None)
        self.Tcache = {}
        self.inv = lambda a: int(gmpy2.invert(mpz(a % p), p))

    def bal(self, r):
        return r - self.p if r > self.p // 2 else r

    def vec_to_arrays(self, v):
        """{(tag,comp): ents} -> {tag: int64 (ncomp, Nt, Ny)}"""
        C, p = self.C, self.p
        out = {}
        for (tag, j), ents in v.items():
            a = out.get(tag)
            if a is None:
                a = out[tag] = np.zeros((C.ncomp, C.Nt, C.Ny), np.int64)
            for it, iy, c in ents:
                a[j, it, iy] = (a[j, it, iy] + rmod(c, p)) % p
        return out

    # ---- the shift matrix of one kernel class --------------------------
    def T(self, dname, lt, cls, ct, cy, rt, ry):
        key = (dname, lt, cls, ct, cy)
        T = self.Tcache.get(key)
        if T is not None:
            return T
        C, p = self.C, self.p
        Nt, Ny = C.Nt, C.Ny
        maxt = C.ORDt + 1 - rt          # max grid index allowed for result parity
        maxy = C.ORDy + 1 - ry
        T = np.zeros((Nt, Ny, Nt, Ny))
        for kt, ky, c in C.kern[(dname, lt)][cls]:
            # src index s (n1 = s-1), dst = s + (kt-1) + carry
            sh_t, sh_y = kt - 1 + ct, ky - 1 + cy
            s0t, s1t = max(0, -sh_t), min(Nt, maxt + 1 - sh_t)
            s0y, s1y = max(0, -sh_y), min(Ny, maxy + 1 - sh_y)
            if s0t >= s1t or s0y >= s1y:
                continue
            val = self.bal(rmod(c, p))
            it = np.arange(s0t, s1t); iy = np.arange(s0y, s1y)
            T[it[:, None], iy[None, :], (it + sh_t)[:, None], (iy + sh_y)[None, :]] += val
        T = T.reshape(C.S, C.S)
        T = np.fmod(T, p)   # sums of several terms hitting same cell
        self.Tcache[key] = T
        return T

    # ---- one multiplication  (atilde /. dlog_d<dname>) . J ---------------
    def multiply(self, J, tags, X2, dname):
        """X2: float64 (ncomp, ntag*S) balanced residues of J."""
        C, p = self.C, self.p
        n, S = C.ncomp, C.S
        ntag = len(tags)
        acc = {}

        def add(tag, arr):
            a = acc.get(tag)
            if a is None:
                acc[tag] = arr
            else:
                a += arr
                np.remainder(a, p, out=a)

        for lt in C.letters:
            kern = C.kern.get((dname, lt))
            if not kern:
                continue
            R, Im = self.Am[lt]
            for part, M in ((0, R), (1, Im)):
                if M is None:
                    continue
                V = M @ X2                                     # exact (< 2^53)
                V -= p * np.round(V / p)                       # balanced, |V| <= p/2
                V = V.reshape(n, ntag, S).transpose(1, 0, 2)
                # group tags by parity
                byp = {}
                for ti, tg in enumerate(tags):
                    byp.setdefault((tg[0], tg[1]), []).append(ti)
                for (pt, py), tis in byp.items():
                    Vp = np.ascontiguousarray(V[tis]).reshape(len(tis) * n, S)
                    for cls in kern:
                        qt, qy, ki = cls
                        rt, ct = (pt + qt) % 2, (pt + qt) // 2
                        ry, cy = (py + qy) % 2, (py + qy) // 2
                        T = self.T(dname, lt, cls, ct, cy, rt, ry)
                        W = Vp @ T                             # exact (< 2^53)
                        W = np.remainder(W, p).astype(np.int64).reshape(len(tis), n, C.Nt, C.Ny)
                        m = part + ki                  # total power of I picked up
                        for q, ti in enumerate(tis):
                            _, _, la, lb, ck = tags[ti]
                            if m == 0:
                                ck2, sign = ck, 1
                            elif m == 1:
                                ck2, sign = key_mul(ck, I_KEY)
                            else:                      # I^2 = -1
                                ck2, sign = ck, -1
                            arr = W[q] if sign == 1 else (p - W[q]) % p
                            add((rt, ry, la, lb, ck2), arr.copy())
        return acc

    # ---- closed-form termwise integration along one axis ----------------
    def integrate(self, acc, axis):
        C, p = self.C, self.p
        ORD = C.ORDt if axis == 1 else C.ORDy
        N = C.Nt if axis == 1 else C.Ny
        out = {}

        def add(tag, arr):
            a = out.get(tag)
            if a is None:
                out[tag] = arr
            else:
                a += arr; np.remainder(a, p, out=a)

        for tag, arr in acc.items():
            pt, py, la, lb, ck = tag
            par = pt if axis == 1 else py
            k = la if axis == 1 else lb
            # e+1 = (2n + par + 2)/2 with n = idx-1  ->  (2 idx + par)/2
            idx = np.arange(N)
            num = 2 * idx + par                    # 2(e+1)
            invf = np.array([0 if x == 0 else (2 * self.inv(int(x))) % p for x in num], np.int64)
            fall = math.factorial(k)
            powv = np.ones(N, np.int64)
            for j in range(k + 1):
                powv = (powv * invf) % p           # (1/(e+1))^(j+1)
                coef = ((-1) ** j) * (fall // math.factorial(k - j)) % p
                f = (powv * coef) % p
                res = np.zeros_like(arr)
                sl_src = [slice(None)] * 3; sl_dst = [slice(None)] * 3
                sl_src[axis] = slice(0, N - 1); sl_dst[axis] = slice(1, N)
                shp = [1, 1, 1]; shp[axis] = N - 1
                res[tuple(sl_dst)] = (arr[tuple(sl_src)] * f[:N - 1].reshape(shp)) % p
                if par == 1:     # e+1 half-integer: top slot n=ORD is beyond ORD
                    top = [slice(None)] * 3; top[axis] = N - 1
                    res[tuple(top)] = 0
                ntag = (pt, py, la - j, lb, ck) if axis == 1 else (pt, py, la, lb - j, ck)
                add(ntag, res)
            if par == 0:          # e = -1 slot (idx 0) -> Log^(k+1)/(k+1) at e=0 (idx 1)
                s0 = [slice(None)] * 3; s0[axis] = 0
                col = arr[tuple(s0)]
                if col.any():
                    res = np.zeros_like(arr)
                    s1 = [slice(None)] * 3; s1[axis] = 1
                    res[tuple(s1)] = (col * self.inv(k + 1)) % p
                    ntag = (pt, py, la + 1, lb, ck) if axis == 1 else (pt, py, la, lb + 1, ck)
                    add(ntag, res)
        return out

    def dy(self, acc):
        """-D[acc, BIG]  (D of y^e Log^b = e y^(e-1) Log^b + b y^(e-1) Log^(b-1))."""
        C, p = self.C, self.p
        N = C.Ny
        out = {}

        def add(tag, arr):
            a = out.get(tag)
            if a is None:
                out[tag] = arr
            else:
                a += arr; np.remainder(a, p, out=a)

        inv2 = self.inv(2)
        for tag, arr in acc.items():
            pt, py, la, lb, ck = tag
            assert not arr[:, :, 0].any(), "BIG^-1 in intt"
            idx = np.arange(N)
            e = ((2 * (idx - 1) + py) * inv2) % p           # exponent e mod p
            res = np.zeros_like(arr)
            res[:, :, :N - 1] = (arr[:, :, 1:] * ((p - e[1:]) % p)[None, None, :]) % p
            add(tag, res)
            if lb:
                res = np.zeros_like(arr)
                res[:, :, :N - 1] = (arr[:, :, 1:] * ((p - lb) % p)) % p
                add((pt, py, la, lb - 1, ck), res)
        return out

    # ---- one full recursion step ---------------------------------------
    def step(self, J, bdy):
        C, p = self.C, self.p
        tags = sorted(t for t, a in J.items() if a.any())
        X = np.stack([J[t] for t in tags])                  # (ntag, n, Nt, Ny)
        Xb = np.where(X > p // 2, X - p, X).astype(np.float64)
        X2 = np.ascontiguousarray(Xb.transpose(1, 0, 2, 3).reshape(C.ncomp, len(tags) * C.S))
        tm = time.time()
        multT = self.multiply(J, tags, X2, 't')
        self.Tcache.clear()                 # shift matrices are S x S: free them
        multY = self.multiply(J, tags, X2, 'y')
        self.Tcache.clear()
        tmul = time.time() - tm
        intt = self.integrate(multT, 1)
        tointy = self.dy(intt)
        for tg, a in multY.items():
            b = tointy.get(tg)
            if b is None:
                tointy[tg] = a
            else:
                b += a; np.remainder(b, p, out=b)
        toaddy = self.integrate(tointy, 2)
        res = intt
        for src in (toaddy, bdy):
            for tg, a in src.items():
                b = res.get(tg)
                if b is None:
                    res[tg] = a.copy()
                else:
                    b += a; np.remainder(b, p, out=b)
        res = {t: a for t, a in res.items() if a.any()}
        return res, tmul


# ==========================================================================
#   reconstruction
# ==========================================================================
def ratrec(a, m):
    """rational reconstruction of a mod m (Wang), returns mpq or None."""
    a = mpz(a) % m
    if a == 0:
        return mpq(0)
    bound = gmpy2.isqrt(m // 2)
    r0, r1, s0, s1 = m, a, mpz(0), mpz(1)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or gmpy2.gcd(r1, s1) != 1:
        return None
    return mpq(r1, s1)


NV = 2          # verification primes (false acceptance ~ p^-2 per coefficient)

_REC = {}       # fork-shared state for parallel reconstruction


def _rec_first_pass(tags):
    """rational reconstruction for a chunk of tags (runs in a worker)."""
    results, M, crt, PV = _REC['results'], _REC['M'], _REC['crt'], _REC['PV']
    out, failed, dens = {}, [], {}
    for tag in tags:
        arrs = [r.get(tag) for r in results]
        shape = next(a for a in arrs if a is not None).shape
        arrs = [a if a is not None else np.zeros(shape, np.int64) for a in arrs]
        st = np.stack(arrs)
        nz = np.argwhere(np.any(st != 0, axis=0))
        if not len(nz):
            continue
        vals = st[:, nz[:, 0], nz[:, 1], nz[:, 2]].T.tolist()
        for (j, it, iy), rs in zip(nz.tolist(), vals):
            v = mpz(0)
            for c, r in zip(crt, rs[:-NV]):
                v += c * r
            v %= M
            q_ = ratrec(v, M)
            if q_ is None or any(rmod(q_, pv) != r for pv, r in zip(PV, rs[-NV:])):
                failed.append((tag, j, it, iy, v, rs[-NV:]))
                continue
            if q_ != 0:
                out[(tag, j, it, iy)] = q_
                d = q_.denominator
                if d != 1:
                    dj = dens.get(j, mpz(1))
                    dens[j] = dj * d // gmpy2.gcd(dj, d)
    return out, failed, dens


def reconstruct(results, primes, jobs=1):
    """results: list (per prime, same order as primes) of {tag: int64 array}.
    Returns ({(tag, comp, it, iy): mpq}, n_unverified, n_second_pass).
    All but the last NV primes reconstruct, the last NV verify.
    Pass 1: rational reconstruction (numerator AND denominator unknown).
    Pass 2 (COMMON DENOMINATOR): for entries pass 1 could not fix, use the
    component's lcm of denominators D: the entry is N/D with N an integer,
    which plain CRT recovers with about half the bits."""
    import multiprocessing as mproc
    P = primes[:-NV]; PV = primes[-NV:]
    M = mpz(1)
    for q in P:
        M *= q
    crt = []
    for q in P:
        Mi = M // q
        crt.append(Mi * gmpy2.invert(Mi % q, q) % M)
    alltags = set()
    for r in results:
        alltags |= set(r)
    alltags = sorted(alltags)
    _REC.update(results=results, M=M, crt=crt, PV=PV)
    if jobs > 1 and len(alltags) > 1:
        chunks = [alltags[i::jobs * 4] for i in range(jobs * 4)]
        with mproc.get_context('fork').Pool(jobs) as pool:
            parts = pool.map(_rec_first_pass, [c for c in chunks if c])
    else:
        parts = [_rec_first_pass(alltags)]
    _REC.clear()
    out, failed, dens = {}, [], {}
    for o, f, d in parts:
        out.update(o); failed += f
        for j, dj in d.items():
            e = dens.get(j, mpz(1))
            dens[j] = e * dj // gmpy2.gcd(e, dj)
    bad = 0
    half = M // 2
    for tag, j, it, iy, v, rv in failed:
        D = dens.get(j, mpz(1))
        n = v * D % M
        if n > half:
            n -= M
        q_ = mpq(n, D)
        if any(rmod(q_, pv) != r for pv, r in zip(PV, rv)):
            bad += 1
            continue
        if q_ != 0:
            out[(tag, j, it, iy)] = q_
    return out, bad, len(failed)


def exact_to_slots(C, v):
    """{(tag,comp): ents} -> {(tag, comp, it, iy): mpq} (merged)."""
    out = {}
    for (tag, j), ents in v.items():
        for it, iy, c in ents:
            k = (tag, j, it, iy)
            out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c != 0}
