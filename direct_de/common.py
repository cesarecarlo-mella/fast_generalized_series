"""common data for the direct DE solver (H family, Euclidean region)."""
import os, re
import numpy as np
import mparse as mp

ACTIVE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16]
ZETA = {3: 1.2020569031595942853997, 5: 1.0369277551433699263314}


def const_value(key):
    v = 1.0 + 0j
    for at, e in key:
        if at == 'Pi': v *= np.pi ** (e // 2)
        elif at == 'I': v *= (1j) ** (e // 2)
        elif at.startswith('Zeta['): v *= ZETA[int(at[5:-1])] ** (e // 2)
        elif at == 'Log[2]': v *= np.log(2) ** (e // 2)
        else: raise ValueError(at)
    return v


def load_atilde(path):
    """-> {k: dense complex (n, n)}"""
    A = mp.parse_file(path)
    n = len(A)
    M = {k: np.zeros((n, n), complex) for k in ACTIVE}
    for i, row in enumerate(A):
        for j, e in enumerate(row):
            for key, v in e.items():
                d = dict(key); im = d.pop('I', 0); (lt, _), = d.items()
                M[int(lt[2:-1])][i, j] += float(v) * (1j if im else 1)
    return M


def load_vector(path):
    vec = mp.parse_file(path)
    return np.array([sum(float(c) * const_value(k) for k, c in p.items()) for p in vec], complex)


# ------------------------------------------------------------- exact helpers
from fractions import Fraction


def load_atilde_exact(path):
    """-> {k: {(i, j): (Fraction re, Fraction im)}}  (sparse, exact)"""
    A = mp.parse_file(path)
    M = {k: {} for k in ACTIVE}
    for i, row in enumerate(A):
        for j, e in enumerate(row):
            for key, v in e.items():
                d = dict(key); im = d.pop('I', 0); (lt, _), = d.items()
                k = int(lt[2:-1]); q = Fraction(int(v.numerator), int(v.denominator))
                re0, im0 = M[k].get((i, j), (Fraction(0), Fraction(0)))
                M[k][(i, j)] = (re0, im0 + q) if im else (re0 + q, im0)
    return M


def load_vector_exact(path):
    """-> {constant key (without I): [complex Fraction pairs per component]} as {ckey: {j: (re, im)}}"""
    vec = mp.parse_file(path)
    out = {}
    for j, p in enumerate(vec):
        for key, v in p.items():
            d = dict(key); im = d.pop('I', 0) // 2
            ck = tuple(sorted(d.items())); q = Fraction(int(v.numerator), int(v.denominator))
            sign = 1 if im % 4 in (0, 1) else -1
            re_, im_ = (q * sign, Fraction(0)) if im % 2 == 0 else (Fraction(0), q * sign)
            slot = out.setdefault(ck, {})
            a, b = slot.get(j, (Fraction(0), Fraction(0)))
            slot[j] = (a + re_, b + im_)
    return out


def matvec_exact(M, vec):
    """sparse exact complex matrix times exact vector {ckey: {j: (re, im)}}"""
    out = {}
    for ck, comp in vec.items():
        r = {}
        for (i, j), (mr, mi) in M.items():
            if j in comp:
                vr, vi = comp[j]
                a, b = r.get(i, (Fraction(0), Fraction(0)))
                r[i] = (a + mr * vr - mi * vi, b + mr * vi + mi * vr)
        r = {i: z for i, z in r.items() if z[0] or z[1]}
        if r:
            out[ck] = r
    return out


def to_numeric(vec, n=371):
    out = np.zeros(n, complex)
    for ck, comp in vec.items():
        cv = const_value(ck)
        for j, (a, b) in comp.items():
            out[j] += (float(a) + 1j * float(b)) * cv
    return out
