#!/usr/bin/env python3
"""
run.py -- drive the array-based blow-up recursion for one corner, weights
1..W, modulo several primes in parallel, reconstruct exact rationals, and
(optionally) compare against a reference set of exact_c<c>_w<w>.m files.

    python3 run.py --data DIR --corner 1 --ordt 10 --ordy 30 --wmax 4 \
                   --nprimes 6 --jobs 4 [--ref DIR] [--write OUTDIR]

--data : directory with atilde.m, dlog_dt/dy_c<c>.m, sol0.m, bdy_c<c>_w<w>.m
         (the dlog files may be at a HIGHER order than --ordt/--ordy: series
          terms beyond the requested order are simply dropped on load)
--ref  : directory with reference exact_c<c>_w<w>.m (any order >= ours;
         truncated on load to ours) -- compared exactly, slot by slot
--write: write exact_c<c>_w<w>.m in the pipeline's own format
"""
import argparse
import os
import sys
import time
import pickle
import multiprocessing as mproc

ap = argparse.ArgumentParser()
ap.add_argument('--data', required=True)
ap.add_argument('--corner', type=int, default=1)
ap.add_argument('--ordt', type=int, default=10)
ap.add_argument('--ordy', type=int, default=30)
ap.add_argument('--wmax', type=int, default=4)
ap.add_argument('--nprimes', type=int, default=6, help='initial number of primes (one is kept for verification)')
ap.add_argument('--jobs', type=int, default=os.cpu_count())
ap.add_argument('--ref', default=None)
ap.add_argument('--write', default=None)
ap.add_argument('--scratch', default='scratch')
ap.add_argument('--addprimes', type=int, default=4)
ap.add_argument('--phys-ref', default=None, help='dir with phys_c<c>_w<w>_ord<N>.m to compare against')
ap.add_argument('--phys-ord', type=int, default=20)
ap.add_argument('--phys-win', type=int, nargs=2, default=None, help='SMALL and BIG max degree for the comparison')
args = ap.parse_args()
if args.phys_win is None:
    # largest window fully determined by a (ORDt, ORDy) rectangle after t = SMALL/BIG^2
    args.phys_win = [args.ordt, args.ordy - 2 * args.ordt]

# one BLAS thread per worker when running several primes side by side
if args.jobs > 1:
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    os.environ.setdefault('MKL_NUM_THREADS', '1')

import numpy as np
import blowup_array as ba
import mparse as mp

C = None
SOL0 = None
BDY = None


def worker(p):
    """full weight chain 1..wmax modulo one prime; results pickled to disk."""
    t0 = time.time()
    E = ba.ModEngine(C, p)
    J = E.vec_to_arrays(SOL0)
    tinit = time.time() - t0
    times = []
    for w in range(1, args.wmax + 1):
        ts = time.time()
        J, tmul = E.step(J, E.vec_to_arrays(BDY[w]))
        times.append((w, time.time() - ts, tmul, len(J)))
        with open(os.path.join(args.scratch, 'J_c%d_w%d_p%d.pkl' % (C.c, w, p)), 'wb') as f:
            pickle.dump(J, f, protocol=4)
    return p, tinit, times


def main():
    global C, SOL0, BDY
    os.makedirs(args.scratch, exist_ok=True)
    T0 = time.time()
    C = ba.Corner(args.data, args.corner, args.ordt, args.ordy)
    SOL0 = C.load_vector(os.path.join(args.data, 'sol0.m'))
    BDY = {w: C.load_vector(os.path.join(args.data, 'bdy_c%d_w%d.m' % (args.corner, w)))
           for w in range(1, args.wmax + 1)}
    tload = time.time() - T0
    print("[run] corner %d  order (%d,%d)  weights 1..%d  load %.1fs" %
          (args.corner, args.ordt, args.ordy, args.wmax, tload), flush=True)

    primes = ba.gen_primes(200)
    used = []
    ctx = mproc.get_context('fork')
    batch = primes[:args.nprimes]
    per_prime = {}
    recon, recon_info = {}, {}
    todo = list(range(1, args.wmax + 1))
    twall_comp = trec = 0.0
    while True:
        tc = time.time()
        with ctx.Pool(min(args.jobs, len(batch))) as pool:
            for p, tinit, times in pool.imap_unordered(worker, batch):
                per_prime[p] = times
                print("  prime %d: " % p + "  ".join("w%d %.1fs" % t[:2] for t in times), flush=True)
        twall_comp += time.time() - tc
        used += batch
        # reconstruct the not-yet-verified weights (2 primes kept for verification)
        tr = time.time()
        for w in list(todo):
            res = []
            for p in used:
                with open(os.path.join(args.scratch, 'J_c%d_w%d_p%d.pkl' % (C.c, w, p)), 'rb') as f:
                    res.append(pickle.load(f))
            out, bad, ncd = ba.reconstruct(res, used)
            print("  reconstruct w%d with %d(+2) primes: %d coefficients, %d via common denominator, %d unverified"
                  % (w, len(used) - 2, len(out), ncd - bad, bad), flush=True)
            if not bad:
                recon[w] = out; recon_info[w] = (len(used), ncd)
                todo.remove(w)
                with open(os.path.join(args.scratch, 'recon_c%d_w%d.pkl' % (C.c, w)), 'wb') as f:
                    pickle.dump(out, f, protocol=4)
        trec += time.time() - tr
        if not todo:
            break
        batch = primes[len(used):len(used) + args.addprimes]
        print("  -> weights %s need more primes, adding %d" % (todo, len(batch)), flush=True)

    Ttot = time.time() - T0
    print("\n=== TIMING (corner %d, order %d/%d, %d primes of %d bits) ===" %
          (args.corner, args.ordt, args.ordy, len(used), ba.PRIMES_BITS))
    print("  load/parse inputs          : %.1fs" % tload)
    tot_cpu = 0
    for w in range(1, args.wmax + 1):
        tt = [dict((x[0], x[1]) for x in per_prime[p])[w] for p in used]
        tot_cpu += sum(tt)
        print("  weight %d : %6.2fs per prime  x %d primes = %7.1fs CPU   (verified with %d primes)"
              % (w, sum(tt) / len(tt), len(used), sum(tt), recon_info[w][0]))
    print("  recursion CPU total        : %.1fs   (wall %.1fs on %d jobs)" % (tot_cpu, twall_comp, args.jobs))
    print("  reconstruction             : %.1fs" % trec)
    print("  TOTAL wall                 : %.1fs" % Ttot)

    if args.phys_ref:
        import phys
        X, Y = args.phys_win
        print("\n=== PHYSICAL-VARIABLE COMPARISON vs %s  (window SMALL<=%d, BIG<=%d) ===" % (args.phys_ref, X, Y))
        for w in range(1, args.wmax + 1):
            f = os.path.join(args.phys_ref, 'phys_c%d_w%d_ord%d.m' % (args.corner, w, args.phys_ord))
            if not os.path.isfile(f):
                print("  w%d: no file %s" % (w, f)); continue
            tr = time.time()
            ref, kept, total = phys.load_phys_ref(f, args.corner, X, Y)
            mine = phys.slots_to_phys(C, recon[w], X, Y)
            onlyr, onlym, diff = phys.compare(mine, ref)
            bad = sorted(set(k[0] + 1 for k in onlyr + onlym + diff))
            print("  w%d: ref %d coeffs (%d/%d terms in window), mine %d | missing %d, extra %d, different %d | bad comps %s (%.0fs)"
                  % (w, len(ref), kept, total, len(mine), len(onlyr), len(onlym), len(diff), bad[:20], time.time() - tr), flush=True)
            for k in (onlyr + onlym + diff)[:5]:
                print("     e.g.", k, ref.get(k), mine.get(k))

    # --- compare with reference
    if args.ref:
        print("\n=== EXACT COMPARISON vs reference %s ===" % args.ref)
        for w in range(1, args.wmax + 1):
            f = os.path.join(args.ref, 'exact_c%d_w%d.m' % (args.corner, w))
            if not os.path.isfile(f):
                print("  w%d: no reference file" % w); continue
            tr = time.time()
            ref = ba.exact_to_slots(C, C.load_vector(f))
            mine = recon[w]
            onlyr = [k for k in ref if k not in mine]
            onlym = [k for k in mine if k not in ref]
            diff = [k for k in ref if k in mine and ref[k] != mine[k]]
            comps_bad = sorted(set(k[1] + 1 for k in onlyr + onlym + diff))
            print("  w%d: ref %d coeffs, mine %d | missing %d, extra %d, different %d | bad comps %s  (%.1fs)"
                  % (w, len(ref), len(mine), len(onlyr), len(onlym), len(diff),
                     comps_bad[:20], time.time() - tr), flush=True)
            for k in (onlyr + onlym + diff)[:5]:
                print("     e.g.", k, ref.get(k), mine.get(k))

    if args.write:
        os.makedirs(args.write, exist_ok=True)
        for w in range(1, args.wmax + 1):
            write_exact(recon[w], os.path.join(args.write, 'exact_c%d_w%d.m' % (args.corner, w)))
        print("wrote exact files to", args.write)


def write_exact(slots, path):
    """slots -> Mathematica list of 371 expressions (same content as the
    pipeline's exact_c<c>_w<w>.m; term order differs, which Get[] ignores)."""
    comps = [dict() for _ in range(C.ncomp)]
    for (tag, j, it, iy), q in slots.items():
        pt, py, la, lb, ck = tag
        e2t = 2 * (it - 1) + pt; e2y = 2 * (iy - 1) + py
        key = list(ck)
        if e2t: key.append((C.tv, e2t))
        if e2y: key.append((C.yv, e2y))
        if la: key.append((C.ltv, 2 * la))
        if lb: key.append((C.lyv, 2 * lb))
        comps[j][tuple(key)] = q
    with open(path, 'w') as f:
        f.write("{" + ",\n ".join(mp.fmt_poly(c) for c in comps) + "}\n")


if __name__ == '__main__':
    main()
