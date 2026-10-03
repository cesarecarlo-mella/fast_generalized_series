import sys, pickle; sys.path.insert(0,'.')
import numpy as np
exec(open('test_c1.py').read().split("S=Solver")[0])
from common import *
from corner import Corner1
I='/mnt/user-data/uploads/series.all/blow_up_array/pipeline/inputs'
Aex=load_atilde_exact('../ccpipe/inputs/atilde.m'); A=load_atilde('../ccpipe/inputs/atilde.m')
C=Corner1(Aex,A,'../ccpipe/inputs/sol0.m',[f'{I}/bdy_c1_w{w}.m' for w in range(1,7)])
pickle.dump(C.C,open('Cab.pkl','wb'))
C1='/mnt/user-data/uploads/series.all/cluster_results/blow_up_form_2_order_20/out'
ser={w:None for w in (2,4)}
for (X,Y,s0,N) in ((1,0.7,0.1,60),(1,0.7,0.05,60),(1,0.7,0.1,100),(1,0.5,0.1,60)):
    J0,_=C.frobenius(X,Y,s0,N)
    out=[]
    for w in (2,4):
        r=eval_c1(f'{C1}/phys_c1_w{w}_ord20.m',[(s0*X,s0*Y)])[:,0]
        d=digits(J0[w],r); out.append("w%d min %.1f med %.1f"%(w,np.nanmin(d),np.nanmedian(d)))
    print((X,Y,s0,N),out)
