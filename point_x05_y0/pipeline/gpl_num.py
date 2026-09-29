"""
gpl_num.py -- numerical G(a1,...,an; x) for letters a_i in {0, 1} and 0 < x < 1
(harmonic polylogarithms), high precision with mpmath.

    G(a1, ..., an; x) = Int_0^x dt/(t - a1) G(a2, ..., an; t),   G(0,...,0; x) = Log[x]^n/n!

Built from the last letter outwards as  sum_{k,m} c[k,m] x^k Log[x]^m  (power
series truncated at x^K; at x = 1/2 the error is ~ 2^-K).

    value(word, x)            -> mpf
    atom_value("G[0, 1, 1/2]") -> mpf      (the atoms used in bdy_c0_w<w>.m)
"""
from functools import lru_cache
from fractions import Fraction
import math
import re

import mpmath as mpm

mpm.mp.dps = 40
K = 200                                   # series order: 2^-200 at x = 1/2


def _integrate(ser, a):
    """ser: {(k, m): c} (as a function of t);  returns Int_0^x dt/(t-a) ser."""
    out = {}

    def add(k, m, v):
        out[(k, m)] = out.get((k, m), 0) + v

    def int_pow_log(n, m, c):             # Int_0^x t^(n-1) Log[t]^m, n >= 1
        f = math.factorial(m)
        for j in range(m + 1):
            add(n, m - j, c * (-1) ** j * (f // math.factorial(m - j)) / mpm.mpf(n) ** (j + 1))

    for (k, m), c in ser.items():
        if a == 0:
            if k == 0:
                add(0, m + 1, c / (m + 1))
            else:
                int_pow_log(k, m, c)
        else:                             # 1/(t-1) = -sum_j t^j  ->  -t^(k+j)
            for j in range(K - k):
                int_pow_log(k + j + 1, m, -c)
    return {km: v for km, v in out.items() if km[0] <= K}


@lru_cache(maxsize=None)
def _series(word):
    if not word:
        return {(0, 0): mpm.mpf(1)}
    return _integrate(_series(word[1:]), word[0])


def value(word, x):
    x = mpm.mpf(x)
    lx = mpm.log(x)
    return sum(c * x ** k * lx ** m for (k, m), c in _series(tuple(word)).items())


_ATOM = re.compile(r'G\[(.*)\]$')


def atom_value(atom):
    """'G[0, 1, 1/2]' -> numeric value; other atoms: Pi, Zeta[n], Log[2], I."""
    m = _ATOM.match(atom)
    if m:
        parts = [s.strip() for s in m.group(1).split(',')]
        x = Fraction(parts[-1])
        return value(tuple(int(p) for p in parts[:-1]), mpm.mpf(x.numerator) / x.denominator)
    if atom == 'Pi':
        return +mpm.pi
    if atom == 'I':
        return mpm.mpc(0, 1)
    m = re.match(r'Zeta\[(\d+)\]$', atom)
    if m:
        return mpm.zeta(int(m.group(1)))
    m = re.match(r'Log\[(\d+)\]$', atom)
    if m:
        return mpm.log(int(m.group(1)))
    raise ValueError("no numerical value for atom %r" % atom)


def monomial_value(key):
    """mparse key ((atom, exp2), ...) -> value (exp2 even for constants)."""
    v = mpm.mpf(1)
    for at, e in key:
        v *= atom_value(at) ** (mpm.mpf(e) / 2)
    return v


if __name__ == '__main__':
    x = mpm.mpf(1) / 2
    checks = [((1,), mpm.log(1 - x)), ((0,), mpm.log(x)), ((0, 1), -mpm.polylog(2, x)),
              ((0, 0, 1), -mpm.polylog(3, x)), ((1, 1), mpm.log(1 - x) ** 2 / 2),
              ((1, 0), mpm.log(x) * mpm.log(1 - x) + mpm.polylog(2, x))]
    for w, ref in checks:
        print(w, mpm.nstr(value(w, x), 25), mpm.nstr(abs(value(w, x) - ref), 3))
