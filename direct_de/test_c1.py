"""tests of the corner-1 direct solver:
 T1 near corner 1 vs the old corner-1 series (weights 1-4)
 T2 exact weight 1, 2 solution (GPLs, eps_{1,2}.m) at points across the triangle"""
import sys, time, pickle; sys.path.insert(0,'.'); sys.path.insert(0,'../ccpipe/check')
import numpy as np
from solver import Solver
from overlap import load_terms
import exact_continued as E

def digits(a,b,chop=1e-4):
    d=np.abs(a-b); r=np.abs(b)
    with np.errstate(divide='ignore'): g=np.minimum(16,-np.log10(d/r+1e-300))
    g=np.where(d<1e-300,16,g); return np.where(r>=chop,g,np.nan)

def eval_c1(path, pts):
    comps=load_terms(path); out=np.zeros((len(comps),len(pts)),complex)
    X=np.array([p[0] for p in pts]); Y=np.array([p[1] for p in pts])
    A={'x':X+0j,'y':Y+0j,'Log[x]':np.log(X)+0j,'Log[y]':np.log(Y)+0j}
    for j,t in enumerate(comps):
        acc=0
        for c,k in t:
            v=c
            for at,e in k:
                if at in ('x','y'): v=v*A[at]**(e/2)
                elif at.startswith('Log'): v=v*A[at]**(e//2)
                elif at=='Pi': v=v*np.pi**(e//2)
                elif at=='Zeta[3]': v=v*1.2020569031595942**(e//2)
                elif at=='I': v=v*(1j)**(e//2)
                else: raise ValueError(at)
            acc=acc+v
        out[j]=acc
    return out

S=Solver('../ccpipe/inputs','/mnt/user-data/uploads/series.all/blow_up_array/pipeline/inputs')
C1='/mnt/user-data/uploads/series.all/cluster_results/blow_up_form_2_order_20/out'
pts1=[(0.12,0.08),(0.05,0.2),(0.2,0.05),(0.25,0.25)]
sol={P:S.evaluate(*P) for P in pts1}
print("T1: direct solver vs corner-1 series (order 20)")
for w in (1,2,3,4):
    ser=eval_c1(f'{C1}/phys_c1_w{w}_ord20.m',pts1)
    for i,P in enumerate(pts1):
        d=digits(sol[P][0][w],ser[:,i]); print("  w%d %s: min %.1f median %.1f digits  (panels %d)"%(w,P,np.nanmin(d),np.nanmedian(d),sol[P][1]['panels']))
print("T2: direct solver vs exact GPL solution")
ex='/mnt/user-data/uploads/series.all/compare_results/exact_solutions'
for P in ((0.3,0.2),(0.6,0.3),(0.1,0.85),(0.45,0.45),(0.8,0.1)):
    J,info=S.evaluate(*P)
    for w in (1,2):
        ref=E.eval_exact(f'{ex}/eps_{w}.m',*P)
        d=digits(J[w],ref); print("  w%d %s: min %.1f median %.1f digits (panels %d, est. rel err %.0e)"%(w,P,np.nanmin(d),np.nanmedian(d),info['panels'],info['err']))
