#!/usr/bin/env python3
"""
compare_petr.py -- exact, coefficient-by-coefficient comparison of Petr's physical
series (JSON, schema 1/3, rectangular truncation u,v <= 20) with ours
(phys_c<c>_w<w>_ord<N>.m from the array pipeline), both truncated to
u, v <= K (half-integer exponents included: e <= K).

Petr's JSON can be read straight from the compressed archive (no unpacking to
disk), or from an unpacked results/ directory:

    zstd -d --long=31 -c series_results.tar.zst | \
        python3 compare_petr.py --tar - --ours ../cluster_results/summary --order 15 --out report_ord15.txt
    python3 compare_petr.py --dir results --ours <dir with <F>/phys_*.m> --order 20 --families H ...

Per (family, corner, weight): number of masters identical / different, and for
the first few different masters what differs (only in ours, only in Petr's,
different coefficient; overall sign flips are recognised).
Pure python, no packages; streams one file at a time.
"""
import argparse
import json
import os
import re
import sys
import tarfile
import time

VARS = {"c1": ("x", "y"), "c2": ("y", "z"), "c3": ("x", "z")}
CONST = {"Pi": 4, "Zeta[3]": 5, "Zeta[5]": 6}     # index in the key


def norm_q(num, den):
    return (int(num), int(den))


def petr_master(row, tokens, half, K):
    """-> {key: (num, den)}, key = (2u, 2v, logu, logv, pi, z3, z5), truncated."""
    out = {}
    for tid, n, d in row:
        u, v, lu, lv, pi, z3, z5 = tokens[tid]
        e2u, e2v = (2 * u + 1, 2 * v + 1) if half else (2 * u, 2 * v)
        if e2u > 2 * K or e2v > 2 * K:
            continue
        out[(e2u, e2v, lu, lv, pi, z3, z5)] = norm_q(n, d)
    return out


_ATOM = re.compile(r'^([A-Za-z]+(?:\[[^\]]*\])?)(?:\^(?:\((-?\d+)(?:/(\d+))?\)|(\d+)))?$')


def our_master(s, corner, K, unknown):
    """one component string written by fmt_poly: '(q)*a*b^2 + (q)*...' or '0'."""
    s = s.strip()
    if s == "0" or not s:
        return {}
    u, v = VARS[corner]
    lu, lv = "Log[%s]" % u, "Log[%s]" % v
    out = {}
    for term in s.split(" + "):
        term = term.strip()
        assert term.startswith("("), term[:80]
        j = term.index(")")
        q = term[1:j]
        num, den = (q.split("/") + ["1"])[:2]
        key = [0] * 7
        rest = term[j + 1:]
        if rest:
            assert rest[0] == "*", term[:80]
            for f in rest[1:].split("*"):
                m = _ATOM.match(f)
                if not m:
                    unknown.add(f); key = None; break
                at = m.group(1)
                if m.group(4):                    # a^k
                    e2 = 2 * int(m.group(4))
                elif m.group(2):                  # a^(n) or a^(n/2)
                    e2 = int(m.group(2)) if m.group(3) == "2" else 2 * int(m.group(2))
                else:                             # a
                    e2 = 2
                if at == u:
                    key[0] += e2
                elif at == v:
                    key[1] += e2
                elif at == lu:
                    key[2] += e2 // 2
                elif at == lv:
                    key[3] += e2 // 2
                elif at in CONST:
                    key[CONST[at]] += e2 // 2
                else:
                    unknown.add(at); key = None; break
        if key is None:
            continue
        if key[0] > 2 * K or key[1] > 2 * K:
            continue
        out[tuple(key)] = norm_q(num, den)
    return out


def iter_top_list(s):
    s = s.strip()
    assert s[0] == "{" and s[-1] == "}"
    depth, start = 0, 1
    for m in re.finditer(r"[(){}\[\],]", s):
        c = m.group(0)
        if m.start() == 0 or m.start() == len(s) - 1:
            continue
        if c in "({[":
            depth += 1
        elif c in ")}]":
            depth -= 1
        elif depth == 0:
            yield s[start:m.start()]
            start = m.end()
    yield s[start:len(s) - 1]


def compare(data, ours_path, K, nshow=3):
    corner, w = data["corner"], data["weight"]
    tokens = [tuple(t) for t in data["tokens"]]
    pref = data.get("master_prefactors", ["1"] * data["dimension"])
    s = open(ours_path).read()
    s = re.sub(r"\\\r?\n\s*", "", s)
    ours = list(iter_top_list(s))
    rows = data["rows"]
    if len(ours) != len(rows):
        return "DIMENSION MISMATCH ours %d vs Petr %d" % (len(ours), len(rows)), []
    unknown = set()
    nsame = ndiff = nsign = 0
    show = []
    for m, (row, os_) in enumerate(zip(rows, ours)):
        P = petr_master(row, tokens, pref[m] == "sqrt_uv", K)
        O = our_master(os_, corner, K, unknown)
        if P == O:
            nsame += 1
            continue
        if P and set(P) == set(O) and all(O[k] == (-P[k][0], P[k][1]) for k in P):
            nsign += 1
            if len(show) < nshow:
                show.append("  master %d (prefactor %s): OVERALL SIGN FLIP (%d terms)" % (m + 1, pref[m], len(P)))
            continue
        ndiff += 1
        if len(show) < nshow:
            onlyP = [k for k in P if k not in O]
            onlyO = [k for k in O if k not in P]
            dq = [k for k in P if k in O and P[k] != O[k]]
            show.append("  master %d (prefactor %s): %d terms Petr, %d ours | only Petr %d, only ours %d, "
                        "different coefficient %d  e.g. %s"
                        % (m + 1, pref[m], len(P), len(O), len(onlyP), len(onlyO), len(dq),
                           (("Petr-only %s=%s" % (onlyP[0], P[onlyP[0]])) if onlyP else
                            ("ours-only %s=%s" % (onlyO[0], O[onlyO[0]])) if onlyO else
                            ("%s: Petr %s ours %s" % (dq[0], P[dq[0]], O[dq[0]])))))
    status = "IDENTICAL" if ndiff == 0 and nsign == 0 else "DIFFERENT"
    msg = "%s  (%d masters: %d identical, %d sign-flipped, %d different)" % (status, len(rows), nsame, nsign, ndiff)
    if unknown:
        msg += "  !! unknown atoms in ours: %s" % sorted(unknown)[:5]
    return msg, show


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--tar", help="tar stream ('-' = stdin) of Petr's results/")
    src.add_argument("--dir", help="unpacked Petr results/ directory")
    ap.add_argument("--ours", required=True, help="dir with <FAMILY>/phys_c<c>_w<w>_ord<N>.m")
    ap.add_argument("--ours-order", type=int, default=None, help="N in our file names (default = --order)")
    ap.add_argument("--order", type=int, required=True, help="compare up to u, v <= K")
    ap.add_argument("--families", default=None, help="comma list (default: all found)")
    ap.add_argument("--corners", default="c1,c2,c3")
    ap.add_argument("--weights", default="1,2,3,4,5,6")
    ap.add_argument("--out", default=None)
    ap.add_argument("--budget", type=float, default=None,
                    help="stop (cleanly, after the current file) once this many seconds have passed")
    ap.add_argument("--resume-after", default=None,
                    help="skip archive members up to and including this file name (e.g. C_c2_w4.json)")
    a = ap.parse_args()
    N = a.ours_order if a.ours_order is not None else a.order
    fams = set(a.families.split(",")) if a.families else None
    out = open(a.out, "a") if a.out else None

    def log(t):
        print(t, flush=True)
        if out:
            out.write(t + "\n"); out.flush()

    def handle(name, fileobj):
        base = os.path.basename(name)
        m = re.match(r"([A-Z0-9]+)_(c[123])_w(\d)\.json$", base)
        if not m:
            return
        F, c, w = m.group(1), m.group(2), int(m.group(3))
        if (fams and F not in fams) or w == 0 or c not in a.corners.split(",") \
                or str(w) not in a.weights.split(","):
            return
        ours = os.path.join(a.ours, F, "phys_%s_w%d_ord%d.m" % (c, w, N))
        if not os.path.isfile(ours):
            log("%-3s %s w%d: no file of ours (%s)" % (F, c, w, ours)); return
        t0 = time.time()
        data = json.load(fileobj)
        msg, show = compare(data, ours, a.order)
        log("%-3s %s w%d: %s  [%.0fs]" % (F, c, w, msg, time.time() - t0))
        for sline in show:
            log(sline)

    T0 = time.time()
    log("=== compare_petr  order <= %d  ours: %s ===" % (a.order, a.ours))
    if a.tar:
        stream = sys.stdin.buffer if a.tar == "-" else open(a.tar, "rb")
        skipping = a.resume_after is not None
        last = None
        with tarfile.open(fileobj=stream, mode="r|") as tf:
            for ti in tf:
                if not (ti.isfile() and ti.name.endswith(".json")):
                    continue
                base = os.path.basename(ti.name)
                if skipping:
                    if base == a.resume_after:
                        skipping = False
                    continue
                handle(ti.name, tf.extractfile(ti))
                last = base
                if a.budget and time.time() - T0 > a.budget:
                    log("RESUME_AFTER %s" % last)
                    return
        if last is None and skipping:
            log("resume point %s not found" % a.resume_after)
    else:
        for root, _, files in sorted(os.walk(a.dir)):
            for f in sorted(files):
                if f.endswith(".json"):
                    with open(os.path.join(root, f)) as fh:
                        handle(f, fh)
    log("=== done ===")


if __name__ == "__main__":
    main()
