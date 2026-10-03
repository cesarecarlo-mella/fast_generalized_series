"""T3: two independent starts (different ray and s0 out of corner 1) must give the
same J at far points, all weights 1-6 (checks Frobenius + path integration +
the regularised log structure, which enters differently along the two rays)."""
import sys, time; sys.path.insert(0,'.')
import numpy as np
from solver import Solver
I='/mnt/user-data/uploads/series.all/blow_up_array/pipeline/inputs'
S1=Solver('../ccpipe/inputs',I,ray=(1.0,0.7),s0=0.1)
S2=Solver('../ccpipe/inputs',I,ray=(0.4,1.0),s0=0.12)
def digits(a,b,chop=1e-4):
    d=np.abs(a-b); r=np.abs(b)
    with np.errstate(divide='ignore',invalid='ignore'): g=np.minimum(16,-np.log10(d/r+1e-300))
    g=np.where(d<1e-300,16,g); return np.where(r>=chop,g,np.nan)
for P in ((0.5,0.3),(0.15,0.8),(0.85,0.1),(0.33,0.33)):
    a,ia=S1.evaluate(*P); b,ib=S2.evaluate(*P)
    print(P, " ".join("w%d %.1f/%.1f"%(w,np.nanmin(digits(a[w],b[w])),np.nanmedian(digits(a[w],b[w]))) for w in range(1,7)), "panels",ia['panels'],ib['panels'])
