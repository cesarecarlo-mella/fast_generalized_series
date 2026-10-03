"""c_w = f_w(P) - J_w(P; c_w = 0) at two matching points; recognise exactly."""
import sys, pickle, numpy as np
from fractions import Fraction
from overlap import load_terms, evaluate
gt = pickle.load(open('gt_match.pkl', 'rb'))
wt = int(sys.argv[1]); workdir = sys.argv[2]; outdir = sys.argv[3]
def recog(z, kind):
    """weight 1: z = q*I*Pi ; weight 2: z = q*Pi^2 (+ check)."""
    if kind == 1:
        q = Fraction(z.imag / np.pi).limit_denominator(20000); approx = complex(0, float(q) * np.pi)
        s = "0" if q == 0 else "(%s)*I*Pi" % q
    else:
        q = Fraction(z.real / np.pi ** 2).limit_denominator(20000); approx = float(q) * np.pi ** 2
        s = "0" if q == 0 else "(%s)*Pi^2" % q
    return s, abs(z - approx)
for c in ('CC2e', 'CC3e'):
    pts = [p for (cc, p) in gt if cc == c]
    J = evaluate(load_terms('%s/phys_%s_w%d_ord12.m' % (workdir, c, wt)), c, np.array(pts))
    cs = [gt[(c, p)][wt - 1] - J[:, i] for i, p in enumerate(pts)]
    agree = np.max(np.abs(cs[0] - cs[1]))
    strs, res = zip(*[recog(z, wt) for z in cs[0]])
    print("%s w%d: constants from the two points agree to %.1e ; recognition residual max %.1e ; nonzero %d/371"
          % (c, wt, agree, max(res), sum(s != '0' for s in strs)))
    open('%s/bdy_%s_w%d.m' % (outdir, c, wt), 'w').write('{' + ', '.join(strs) + '}\n')
