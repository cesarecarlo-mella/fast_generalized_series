#!/usr/bin/env python3
"""
check_map.py -- exact cross-check of the CC2/CC3 letter series against the
OLD corner-2/3 master letters (Mathematica-built, 00_build_letters.wl),
through the alphabet map of ../../alphabet_map/ (times --odd-sign for l7, l8,
since the old letters use the principal root):

    dlog l_k(x(w,u), y(w,u)) = sum_j M_kj dlog l_j(w, u)     (exact, dlog level)

and the relabelling (x, y, z) -> (w, u, v) that turns old c2 (t, y) into
CC2 (t, u) and old c3 (t, z) into CC3 (t, v). Every coefficient up to
t^ORDT B^ORDB must agree exactly.

    python3 check_map.py --new ../inputs --old <dir with dlog_d{t,y}_c2.m, dlog_d{t,y}_c3.m> [--ordt 20 --ordb 60]
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mparse as mp

M = {1: {1: 1, 3: -1}, 2: {2: 1, 3: -1}, 3: {3: -1}, 4: {5: 1, 3: -1}, 5: {4: 1, 3: -1}, 6: {6: 1, 3: -1},
     7: {7: -1, 8: -1}, 8: {8: 1},
     9: {13: 1, 3: -2}, 10: {16: 1, 3: -2}, 11: {12: 1, 3: -2}, 12: {11: 1, 3: -2}, 13: {9: 1, 3: -2},
     16: {10: 1, 3: -2}}

ap = argparse.ArgumentParser()
ap.add_argument('--new', required=True)
ap.add_argument('--old', required=True)
ap.add_argument('--ordt', type=int, default=20)
ap.add_argument('--ordb', type=int, default=60)
ap.add_argument('--odd-sign', type=int, default=-1, choices=[1, -1],
                help='sign used in build_cc_letters.py; the old letters use the principal root (+1)')
a = ap.parse_args()

def load(path, rename):
    out = {}
    for _, lhs, rhs in mp.parse_file(path):
        (lk, _), = next(iter(lhs))
        k = int(lk[2:-1])
        p = {}
        for key, v in rhs.items():
            d = {rename.get(at, at): e for at, e in key}
            if d.get('t', 0) > 2 * a.ordt or max([e for at, e in d.items() if at not in ('t', 'I')] or [0]) > 2 * a.ordb:
                continue
            p[tuple(sorted(d.items()))] = v
        out[k] = p
    return out

ok = True
# CC2 now uses the chart u = t v^2 (u << v), which is NOT the old-c2 shape, so
# only CC3 (= old c3 shape, relabelled) can be checked this way.
for CC, old_c, oldbig, big in (('CC3', 3, 'z', 'v'),):
    for d_new, d_old in (('t', 't'), (big, 'y')):
        new = load(os.path.join(a.new, 'dlog_d%s_%s.m' % (d_new, CC)), {})
        old = load(os.path.join(a.old, 'dlog_d%s_c%d.m' % (d_old, old_c)), {oldbig: big})
        nbad = 0
        for k, combo in M.items():
            pred = {}
            for j, c in combo.items():
                for key, v in old[j].items():
                    pred[key] = pred.get(key, 0) + c * v
            if k in (7, 8):              # odd letters: r -> odd_sign * r flips their dlog
                pred = {kk: a.odd_sign * v for kk, v in pred.items()}
            pred = {kk: v for kk, v in pred.items() if v}
            if pred != new[k]:
                nbad += 1; ok = False
                diff = set(pred) ^ set(new[k]) or {kk for kk in pred if pred[kk] != new[k][kk]}
                print("  %s d%s l%d: MISMATCH (%d keys differ, e.g. %s)" % (CC, d_new, k, len(diff), sorted(diff)[:2]))
        print("%s  d/d%-2s : %d/%d letters identical to sum_j M_kj (old c%d letters)  [%d terms total]"
              % (CC, d_new, len(M) - nbad, len(M), old_c, sum(len(v) for v in new.values())))
print("ALL IDENTICAL" if ok else "MISMATCHES FOUND")
sys.exit(0 if ok else 1)
