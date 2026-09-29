#!/usr/bin/env python3
"""
build_bdy_from_bc.py -- rebuild a family's sol0.m and bdy_c<c>_w<w>.m from
all.ints/boundaries_all/<Family>_bc.m.  This is the Python equivalent of the
boundary part of all_families/<F>/pipeline/00_build_inputs.wl, using the same
recipe:

    <Family>_bc.m = { BC[F,"",1] -> <eps series>, ... }  with 6*Nf entries
    batch b = entries (b-1)Nf+1 .. b Nf
    corner 1 <- batch 1,  corner 2 <- batch 3,  corner 3 <- batch 2
    sol0          = Coefficient[batch1, ep, 0]
    bdy_c<c>_w<w> = Coefficient[batch(c), ep, w]      w = 1..6

    python3 build_bdy_from_bc.py --bc .../I3_bc.m --nf 333 --out I3/inputs [--check DIR]

--check DIR : compare every rebuilt file with the one already in DIR and report
              which files differ; nothing is written unless --out is given.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mparse as mp

BATCH_OF_CORNER = {1: 1, 2: 3, 3: 2}

ap = argparse.ArgumentParser()
ap.add_argument('--bc', required=True)
ap.add_argument('--nf', type=int, default=None, help='component count (default: from --atilde)')
ap.add_argument('--atilde', default=None)
ap.add_argument('--wmax', type=int, default=6)
ap.add_argument('--out', default=None)
ap.add_argument('--check', default=None)
a = ap.parse_args()

nf = a.nf
if nf is None:
    if not a.atilde:
        sys.exit("give --nf or --atilde")
    nf = len(mp.split_top_list(open(a.atilde).read()))
rules = mp.parse_file(a.bc)
if len(rules) != 6 * nf:
    sys.exit("%s has %d entries, expected 6*Nf = %d" % (a.bc, len(rules), 6 * nf))
vals = [r[2] for r in rules]


def coeff(p, w):
    out = {}
    for k, v in p.items():
        d = dict(k)
        if d.pop('ep', 0) != 2 * w:
            continue
        out[tuple(sorted(d.items()))] = v
    return out


files = {'sol0.m': [coeff(vals[i], 0) for i in range(nf)]}
for c in (1, 2, 3):
    b = BATCH_OF_CORNER[c]
    batch = vals[(b - 1) * nf: b * nf]
    for w in range(1, a.wmax + 1):
        files['bdy_c%d_w%d.m' % (c, w)] = [coeff(e, w) for e in batch]

bad_atoms = sorted({at for comps in files.values() for p in comps for k in p for at, _ in k
                    if not (at in ('Pi', 'I') or at.startswith('Zeta['))})
if bad_atoms:
    print("WARNING: unexpected symbols in the boundary constants:", bad_atoms[:10])

changed = []
if a.check:
    for name, comps in files.items():
        f = os.path.join(a.check, name)
        if not os.path.isfile(f):
            print("  %-14s  (no existing file)" % name); changed.append(name); continue
        old = mp.parse_file(f)
        diff = [i + 1 for i in range(nf) if old[i] != comps[i]]
        if diff:
            changed.append(name)
        print("  %-14s  %s" % (name, "identical" if not diff else
                               "DIFFERS in %d components, e.g. %s" % (len(diff), diff[:8])))
if a.out:
    os.makedirs(a.out, exist_ok=True)
    for name, comps in files.items():
        with open(os.path.join(a.out, name), 'w') as fh:
            fh.write("{" + ",\n ".join(mp.fmt_poly(p) for p in comps) + "}\n")
    print("wrote %d files to %s" % (len(files), a.out))
print("files that differ from --check:", changed if a.check else "(not checked)")
