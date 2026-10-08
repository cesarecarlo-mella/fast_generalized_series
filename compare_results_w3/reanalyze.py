#!/usr/bin/env python3
"""
reanalyze.py -- recompute compare_<tag>.m from existing rawdata_<tag>.m with a
different chop threshold (no re-evaluation of any series).

rawdata rows: {x, y, z, corner, diffN[], diffB[], ref[]}. Per component:
    digits = 16                              if |diff| < 1e-16
           = min(16, -log10(|diff|/|ref|))   otherwise
components with |ref| < CHOP are dropped; MIN / MEDIAN over the rest.
(Same formula as 03_analyze.wl / fast_selfconv.py; only CHOP changes.)

    python3 reanalyze.py --chop 1e-6 --out results_w2_chop1e-6 rawdata_*.m
"""
import argparse
import os
import re

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--chop', type=float, default=1e-6)
ap.add_argument('--out', required=True)
ap.add_argument('files', nargs='+')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)


def load_raw(path):
    s = open(path).read()
    head = s[:s.index('"data" -> ')]
    s = s[s.index('"data" -> ') + 10:s.rindex('|>')]
    s = re.sub(r'`[0-9.]*', '', s).replace('*^', 'e').replace('{', '[').replace('}', ']')
    s = re.sub(r'\bI\b', '1j', s.replace('Indeterminate', 'float("nan")'))
    return head, eval(s)


def digits(d, r):
    d, r = np.abs(np.asarray(d, complex)), np.abs(np.asarray(r, complex))
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.minimum(16.0, -np.log10(d / r + 1e-30))
    return np.where(d < 1e-16, 16.0, g)


def num(v):
    return "Indeterminate" if not np.isfinite(v) else repr(float(v)).replace('e', '*^')


for f in a.files:
    head, rows = load_raw(f)
    meta = dict(re.findall(r'"(mode|weight|ordp|ordb)" -> "?([A-Za-z0-9]+)"?', head))
    out, dropped, total = [], 0, 0
    for x, y, z, c, dN, dB, ref in rows:
        keep = np.abs(np.asarray(ref, complex)) >= a.chop
        total += len(ref); dropped += int((~keep).sum())
        if not keep.any():
            continue
        gN = digits(np.asarray(dN, complex)[keep], np.asarray(ref, complex)[keep])
        gB = digits(np.asarray(dB, complex)[keep], np.asarray(ref, complex)[keep])
        out.append("{%s, %s, %s, %d, %s, %s, %s, %s}" % (num(x), num(y), num(z), int(c), num(gN.min()),
                   num(np.median(gN)), num(gB.min()), num(np.median(gB))))
    tag = os.path.basename(f)[len('rawdata_'):-2]
    path = os.path.join(a.out, 'compare_%s.m' % tag)
    with open(path, 'w') as fh:
        fh.write('<|"mode" -> "%s", "weight" -> %s, "ordp" -> %s, "ordb" -> %s, "chopThreshold" -> %s, '
                 '"header" -> {"x", "y", "z", "corner", "digits_native_min", "digits_native_median", '
                 '"digits_bernoulli_min", "digits_bernoulli_median"}, "data" -> {%s}|>\n'
                 % (meta['mode'], meta['weight'], meta['ordp'], meta['ordb'], num(a.chop), ",\n ".join(out)))
    print("%s: %d points, dropped %.1f%% of component values (|ref| < %g)"
          % (os.path.basename(path), len(out), 100 * dropped / total, a.chop), flush=True)
