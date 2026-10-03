"""exact weight-1,2 H masters (eps_{1,2}.m) in the Euclidean triangle, closed-form GPLs (mpmath):
   G(0;x)=log x, G(a;x)=log(1-x/a), G(0,0;x)=log^2 x/2, G(0,b;x)=-Li2(x/b),
   G(a,0;x)=log x log(1-x/a)+Li2(x/a), G(a,a;x)=log^2(1-x/a)/2,
   G(a,b;x)=Li2((b-x)/(b-a))-Li2(b/(b-a))+log(1-x/b) log((x-a)/(b-a))   (a != b, a,b != 0)"""
import re, os
import numpy as np, mpmath as mpm
import mparse as mp
EX = os.path.dirname(os.path.abspath(__file__))     # eps_1.m, eps_2.m live next to this file


def G(args, x):
    if len(args) == 1:
        a, = args
        return mpm.log(x) if a == 0 else mpm.log(1 - x / a)
    a, b = args
    if a == 0 and b == 0: return mpm.log(x) ** 2 / 2
    if a == 0: return -mpm.polylog(2, x / b)
    if b == 0: return mpm.log(x) * mpm.log(1 - x / a) + mpm.polylog(2, x / a)
    if a == b: return mpm.log(1 - x / a) ** 2 / 2
    return (mpm.polylog(2, (b - x) / (b - a)) - mpm.polylog(2, b / (b - a))
            + mpm.log(1 - x / b) * mpm.log((x - a) / (b - a)))


def G_quad(args, x):
    if len(args) == 1: return G(args, x)
    a, b = args
    inner = (lambda s: mpm.log(s) ** 2 / 2 * 0) if False else None
    f = (lambda s: mpm.log(s) / (s - a)) if b == 0 else (lambda s: mpm.log(1 - s / b) / (s - a))
    if a == 0 and b == 0: return mpm.log(x) ** 2 / 2
    return mpm.quad(f, [0, x])


class Exact:
    def __init__(self, w, dps=30):
        self.vec = mp.parse_file('%s/eps_%d.m' % (EX, w)); self.dps = dps

    def __call__(self, x0, y0, Gf=G):
        with mpm.workdps(self.dps):
            env = {'x': mpm.mpf(x0), 'y': mpm.mpf(y0)}; cache = {}
            def atom(at):
                if at not in cache:
                    m = re.fullmatch(r'G\[(.*)\]', at)
                    if m:
                        parts = [q.strip() for q in m.group(1).split(',')]
                        cache[at] = Gf([mpm.mpf(eval(q, {}, env)) for q in parts[:-1]], env[parts[-1]])
                    else:
                        cache[at] = mpm.pi
                return cache[at]
            out = []
            for p in self.vec:
                s = mpm.mpf(0)
                for key, v in p.items():
                    t = mpm.mpf(v.numerator) / v.denominator
                    for at, e in key:
                        t *= (mpm.mpc(0, 1) ** (e // 2)) if at == 'I' else atom(at) ** (e // 2)
                    s += t
                out.append(s)
            return out
