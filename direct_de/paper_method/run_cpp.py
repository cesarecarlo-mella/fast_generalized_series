"""Python driver for the C++ paper-method solver (h_family_{double,dd})."""
import subprocess, os
import numpy as np
import mpmath as mpm

HERE = os.path.dirname(os.path.abspath(__file__))


def run(prec, bdy, error, pts, delta=(0.1, 0.2), timeout=36000, max_time=120):
    exe = os.path.join(HERE, 'h_family_%s' % prec)
    inp = ''.join('%s %s\n' % (x, y) for x, y in pts)     # pass strings; reference must use mpf(string)
    out = subprocess.run(['nice', exe, os.path.join(HERE, 'data'), os.path.join(HERE, bdy), '%g' % error,
                          str(delta[0]), str(delta[1])], input=inp, capture_output=True, text=True,
                         timeout=timeout, check=True, env=dict(os.environ, DIFFEQS_MAX_TIME=str(max_time))).stdout
    res, cur = [], None
    for line in out.splitlines():
        if line.startswith('#'):
            p = line.split()
            if 'FAILED' in line:
                cur = None; res.append(dict(x=p[1], y=p[2], failed=' '.join(p[4:]))); continue
            cur = dict(x=p[1], y=p[2], steps=int(p[3]), evals=int(p[4]), rejected=int(p[5]),
                       time=float(p[6]), vals=[])
            res.append(cur)
        elif cur is not None:
            a, b = line.split()
            cur['vals'].append(mpm.mpc(a, b))
    for r in res:
        if 'vals' in r:
            v = r.pop('vals'); n = 371
            w = (len(v) - 1) // n
            r['J'] = [v[k * n:(k + 1) * n] for k in range(w)]     # J_1 ... J_w
            r['r'] = v[-1]
    return res


def digits(a, b, chop=1e-4):
    """a, b: lists of mpc. digits = -log10 |a-b|/|b| over components with |b| >= chop (cap 32)."""
    out = []
    for u, v in zip(a, b):
        if abs(v) < chop: continue
        d = abs(u - v)
        out.append(32.0 if d == 0 else min(32.0, float(-mpm.log10(d / abs(v)))))
    return np.array(out)
