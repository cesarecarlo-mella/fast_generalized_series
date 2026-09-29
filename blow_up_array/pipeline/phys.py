"""
phys.py -- undo the blow-up on reconstructed slots (same map as
03_to_physical_vars.wl) and compare against phys_c<c>_w<w>_ord<ORD>.m files
(e.g. the cluster results), restricted to a window SMALL <= X, BIG <= Y.

    t = SMALL/BIG^2 :  t^e BIG^m Log[t]^a Log[BIG]^b
        -> SMALL^e BIG^(m-2e) (Log[SMALL] - 2 Log[BIG])^a Log[BIG]^b
"""
import math
import re
import mparse as mp

MAP = {1: ('x', 'y'), 2: ('z', 'y'), 3: ('x', 'z')}


def slots_to_phys(C, slots, X, Y):
    """{(tag, j, it, iy): q} -> {(j, ex2, ey2, lx, ly, ckey): q} in window."""
    sm, bg = MAP[C.c]
    out = {}
    for (tag, j, it, iy), q in slots.items():
        pt, py, la, lb, ck = tag
        e2t = 2 * (it - 1) + pt
        e2y = 2 * (iy - 1) + py - 2 * e2t
        if e2t > 2 * X or e2y > 2 * Y:
            continue
        for k in range(la + 1):
            c = q * math.comb(la, k) * (-2) ** k
            key = (j, e2t, e2y, la - k, lb + k, ck)
            v = out.get(key, 0) + c
            if v:
                out[key] = v
            else:
                out.pop(key, None)
    return out


_SPLIT = re.compile(r'(\s[+-]\s)')


def _terms(comp):
    """split a component string into top-level signed terms (fast)."""
    parts = _SPLIT.split(comp)
    terms, cur, bal = [], parts[0], parts[0].count('(') - parts[0].count(')')
    for i in range(1, len(parts), 2):
        sep, piece = parts[i], parts[i + 1]
        if bal != 0:
            cur += sep + piece
        else:
            terms.append(cur)
            cur = sep.strip() + piece
        bal += piece.count('(') - piece.count(')')
    terms.append(cur)
    return terms


def _deg(term, v):
    """degree (exp2) of symbol v in a monomial term string; None = unknown."""
    s = term.replace('Log[%s]' % v, '')
    d = 0
    for m in re.finditer(r'Sqrt\[%s\]|(?<![A-Za-z\[])%s(?![A-Za-z\]])(\^\((-?\d+)(?:/(\d+))?\)|\^(\d+))?' % (v, v), s):
        if m.group(0).startswith('Sqrt'):
            d += 1
        elif m.group(1) is None:
            d += 2
        elif m.group(4) is not None:
            d += 2 * int(m.group(4))
        else:
            num = int(m.group(2)); den = int(m.group(3) or 1)
            if den not in (1, 2):
                return None
            d += num * 2 // den
    # a symbol in a denominator ("/(x*..." or "/x") flips sign -> unknown
    if re.search(r'/\(?[^()]*\b%s\b' % v, s):
        return None
    return d


def load_phys_ref(path, corner, X, Y):
    sm, bg = MAP[corner]
    lx, ly = 'Log[%s]' % sm, 'Log[%s]' % bg
    with open(path) as f:
        s = f.read()
    s = re.sub(r'\\\r?\n\s*', '', s)
    out = {}
    kept = total = 0
    for j, comp in enumerate(iter_top_list(s)):
        keep = []
        for tm in _terms(comp.strip()):
            total += 1
            dx, dy = _deg(tm, sm), _deg(tm, bg)
            if (dx is not None and dx > 2 * X) or (dy is not None and dy > 2 * Y):
                continue
            keep.append(tm)
        kept += len(keep)
        if not keep:
            continue
        expr = " + ".join("(" + t + ")" for t in keep)
        poly = mp.parse_string(expr)
        for k, v in poly.items():
            d = dict(k)
            ex = d.pop(sm, 0); ey = d.pop(bg, 0)
            a = d.pop(lx, 0) // 2; b = d.pop(ly, 0) // 2
            if ex > 2 * X or ey > 2 * Y:
                continue
            key = (j, ex, ey, a, b, tuple(sorted(d.items())))
            out[key] = out.get(key, 0) + v
    return {k: v for k, v in out.items() if v != 0}, kept, total


def compare(mine, ref):
    onlyr = [k for k in ref if k not in mine]
    onlym = [k for k in mine if k not in ref]
    diff = [k for k in ref if k in mine and ref[k] != mine[k]]
    return onlyr, onlym, diff


def iter_top_list(s):
    """yield the elements of a top-level Mathematica list one at a time
    (keeps memory flat on the multi-hundred-MB phys/exact files)."""
    s = s.strip()
    assert s[0] == '{' and s[-1] == '}', "not a list"
    depth, start, n = 0, 1, len(s) - 1
    find = re.compile(r'[(){}\[\],]').finditer
    for m in find(s, 1, n):
        c = m.group(0)
        if c in '({[':
            depth += 1
        elif c in ')}]':
            depth -= 1
        elif depth == 0:
            yield s[start:m.start()]
            start = m.end()
    yield s[start:n]


def write_phys(C, physd, path):
    """{(j, ex2, ey2, lx, ly, ck): q} -> phys_c<c>_w<w>_ord<ORD>.m (Get[]-able)."""
    sm, bg = MAP[C.c]
    comps = [dict() for _ in range(C.ncomp)]
    for (j, ex2, ey2, lx, ly, ck), q in physd.items():
        key = list(ck)
        if ex2: key.append((sm, ex2))
        if ey2: key.append((bg, ey2))
        if lx: key.append(('Log[%s]' % sm, 2 * lx))
        if ly: key.append(('Log[%s]' % bg, 2 * ly))
        comps[j][tuple(key)] = q
    with open(path, 'w') as f:
        f.write("{" + ",\n ".join(mp.fmt_poly(c) for c in comps) + "}\n")
