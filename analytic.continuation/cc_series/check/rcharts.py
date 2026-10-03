"""helpers for the ratio-line charts CC2r / CC3r (exact files in t, p)."""
import numpy as np
from overlap import load_terms


def tp_of(region, pts):
    u, w, v = pts.T
    if region == 'CC2r':          # x = -1/p, y = -1 - p + t
        p = v / w; t = 1 - u / v + p
    else:                         # y = -1/p, x = -1 - p + t
        p = v / u; t = 1 - w / v + p
    return t, p


def uwv_of(region, t, p):
    if region == 'CC2r':
        x, y = -1 / p, -1 - p + t
    else:
        y, x = -1 / p, -1 - p + t
    v = 1 / (1 - x - y)
    return (-y * v, -x * v, v)


def eval_exact(path, region, pts):
    t, p = tp_of(region, pts)
    comps = load_terms(path)
    A = {'t': t + 0j, 'p': p + 0j, 'Log[t]': np.log(t + 0j), 'Log[p]': np.log(p + 0j)}
    out = np.zeros((len(comps), len(t)), complex)
    for j, terms in enumerate(comps):
        acc = np.zeros(len(t), complex)
        for c, key in terms:
            val = c
            for at, e in key:
                if at in ('t', 'p'):
                    val = val * A[at] ** (e / 2)
                elif at in ('Log[t]', 'Log[p]'):
                    val = val * A[at] ** (e // 2)
                elif at == 'Pi':
                    val = val * np.pi ** (e // 2)
                elif at == 'I':
                    val = val * (1j) ** (e // 2)
                elif at == 'Log[2]':
                    val = val * np.log(2) ** (e // 2)
                elif at == 'Zeta[3]':
                    val = val * 1.2020569031595942 ** (e // 2)
                else:
                    raise ValueError(at)
            acc = acc + val
        out[j] = acc
    return out
