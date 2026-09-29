"""shared paths / helpers for the blow_up_array cluster pipeline."""
import os
import pickle

import blowup_array as ba

PRIME_TABLE = ba.gen_primes(400)      # deterministic: prime index k -> PRIME_TABLE[k]


def setup_path(wd, c):
    return os.path.join(wd, 'setup_c%d.pkl' % c)


def res_dir(wd):
    d = os.path.join(wd, 'residues')
    os.makedirs(d, exist_ok=True)
    return d


def res_path(wd, c, w, k):
    return os.path.join(res_dir(wd), 'J_c%d_w%d_k%03d.pkl' % (c, w, k))


def rec_path(wd, c, w):
    return os.path.join(wd, 'rec_c%d_w%d.pkl' % (c, w))


def exact_path(wd, c, w):
    return os.path.join(wd, 'exact_c%d_w%d.m' % (c, w))


def load_pickle(p):
    with open(p, 'rb') as f:
        return pickle.load(f)


def save_pickle(obj, p):
    tmp = p + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f, protocol=4)
    os.replace(tmp, p)                  # atomic: a killed job never leaves a half file


def save_residues(J, p):
    """residues are < 2^22: store as int32, zlib-compressed (the grids are
    mostly zeros, typically 10-30x smaller than raw int64)."""
    import zlib, numpy as np
    blob = {t: (a.shape, zlib.compress(np.ascontiguousarray(a, dtype=np.int32).tobytes(), 1))
            for t, a in J.items()}
    save_pickle(blob, p)


def load_residues(p):
    import zlib, numpy as np
    blob = load_pickle(p)
    return {t: np.frombuffer(zlib.decompress(b), dtype=np.int32).astype(np.int64).reshape(shp)
            for t, (shp, b) in blob.items()}


def primes_with_residues(wd, c, w):
    d = res_dir(wd)
    pre = 'J_c%d_w%d_k' % (c, w)
    ks = sorted(int(f[len(pre):len(pre) + 3]) for f in os.listdir(d)
                if f.startswith(pre) and f.endswith('.pkl'))
    return ks


def write_exact(C, slots, path):
    """reconstructed slots -> exact_c<c>_w<w>.m (same content as the FORM
    pipeline's file; term order differs, which Get[] does not care about)."""
    import mparse as mp
    comps = [dict() for _ in range(C.ncomp)]
    for (tag, j, it, iy), q in slots.items():
        pt, py, la, lb, ck = tag
        e2t = 2 * (it - 1) + pt
        e2y = 2 * (iy - 1) + py
        key = list(ck)
        if e2t: key.append((C.tv, e2t))
        if e2y: key.append((C.yv, e2y))
        if la: key.append((C.ltv, 2 * la))
        if lb: key.append((C.lyv, 2 * lb))
        comps[j][tuple(key)] = q
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write("{" + ",\n ".join(mp.fmt_poly(c) for c in comps) + "}\n")
    os.replace(tmp, path)
