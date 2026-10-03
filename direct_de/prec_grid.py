"""precision of the direct solver at weight <= 2 on the paper grid, vs the exact GPL solution.
   python3 prec_grid.py GRIDFILE OUT.pkl [i0 i1]"""
import sys, re, time, pickle; sys.path.insert(0, '.')
import numpy as np, mpmath as mpm
import mparse as mp
from solver import Solver
mpm.mp.dps = 20
EX = '/mnt/user-data/uploads/series.all/compare_results/exact_solutions'


def G(args, z):
    if all(a == 0 for a in args):
        return mpm.log(z) ** len(args) / mpm.factorial(len(args))
    if len(args) == 1:
        return mpm.log(1 - z / args[0])
    a1, rest = args[0], args[1:]
    return mpm.quad(lambda s: G(rest, s) / (s - a1), [0, z])


class Exact:
    def __init__(self, w):
        self.vec = mp.parse_file('%s/eps_%d.m' % (EX, w))

    def __call__(self, x0, y0):
        env = {'x': mpm.mpf(x0), 'y': mpm.mpf(y0)}; cache = {}
        def atom(at):
            if at not in cache:
                m = re.fullmatch(r'G\[(.*)\]', at)
                if m:
                    parts = [q.strip() for q in m.group(1).split(',')]
                    cache[at] = G([mpm.mpf(eval(q, {}, env)) for q in parts[:-1]], env[parts[-1]])
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
            out.append(complex(s))
        return np.array(out)


def load_grid(path):
    s = open(path).read().replace('\n', ' ')
    return [tuple(float(eval(q)) for q in t.split(',')) for t in re.findall(r'\{([^{}]*)\}', s)]


if __name__ == '__main__':
    grid = load_grid(sys.argv[1]); out = sys.argv[2]
    i0, i1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (0, len(grid))
    S = Solver('../ccpipe/inputs', '/mnt/user-data/uploads/series.all/blow_up_array/pipeline/inputs')
    E1, E2 = Exact(1), Exact(2)
    res = []
    t0 = time.time()
    for i in range(i0, i1):
        x, y, z = grid[i]
        J, info = S.evaluate(x, y, wmax=2)
        res.append(dict(x=x, y=y, J1=J[1], J2=J[2], ref1=E1(x, y), ref2=E2(x, y), panels=info['panels'], err=info['err']))
        if (i - i0) % 20 == 0:
            print(i, "%.0fs" % (time.time() - t0), flush=True)
            pickle.dump(res, open(out, 'wb'))
    pickle.dump(res, open(out, 'wb'))
    print("done", time.time() - t0)
