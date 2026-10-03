import sys, pickle, numpy as np
from fractions import Fraction
sys.path.insert(0,'.')
import exact_continued as E
from rcharts import uwv_of, eval_exact
base='/tmp/claude-0/-home-claude/62cee7e5-47c1-5a82-adde-578e1e20317b/scratchpad/ccpipe'; ex='/mnt/user-data/uploads/series.all/compare_results/exact_solutions'
wt=int(sys.argv[1]); workdir=sys.argv[2]; outdir=sys.argv[3]
try: gt=pickle.load(open('gt_r.pkl','rb'))
except Exception:
    C=E.Continuer(base+'/inputs/atilde.m',base+'/inputs/sol0.m',ex+'/eps_1.m',ex+'/eps_2.m',M=1000)
    gt={}
    for reg in ('CC2r','CC3r'):
        for (t,p) in ((0.02,0.03),(-0.03,0.02)):
            u,w,v=uwv_of(reg,t,p); r=E.evaluate_adaptive(C,-w/v,-u/v); gt[(reg,(t,p))]=(r[1],r[2],(u,w,v)); print(reg,(t,p),r[3])
    pickle.dump(gt,open('gt_r.pkl','wb'))
def recog(z, kind):
    if kind==1:
        q=Fraction(z.imag/np.pi).limit_denominator(20000); approx=complex(0,float(q)*np.pi); s="0" if q==0 else "(%s)*I*Pi"%q
    else:
        # weight 2 at the ratio point: a Pi^2 + b Log[2]^2 + c I Pi Log[2]
        import mpmath as mpm
        L2=np.log(2)
        c=Fraction(z.imag/(np.pi*L2)).limit_denominator(20000)
        if abs(z.real)<1e-9: a=b=Fraction(0)
        else:
            rel=mpm.pslq([z.real,float(mpm.pi**2),L2**2],tol=1e-9,maxcoeff=10**5,maxsteps=10**6)
            if rel is None or rel[0]==0: return "FAIL", 1e9
            a=Fraction(-rel[1],rel[0]); b=Fraction(-rel[2],rel[0])
        approx=float(a)*np.pi**2+float(b)*L2**2+1j*float(c)*np.pi*L2
        parts=[t for t,q in (("(%s)*Pi^2"%a,a),("(%s)*Log[2]^2"%b,b),("(%s)*I*Pi*Log[2]"%c,c)) if q!=0]
        s=" + ".join(parts) if parts else "0"
    return s, abs(z-approx)
for reg in ('CC2r','CC3r'):
    keys=[k for k in gt if k[0]==reg]
    pts=np.array([gt[k][2] for k in keys])
    J=eval_exact('%s/exact_%s_w%d.m'%(workdir,reg,wt),reg,pts)
    cs=[gt[k][wt-1]-J[:,i] for i,k in enumerate(keys)]
    strs,res=zip(*[recog(z,wt) for z in cs[0]])
    print("%s w%d: two points agree to %.1e ; recognition residual %.1e ; nonzero %d"%(reg,wt,np.max(np.abs(cs[0]-cs[1])),max(res),sum(s!='0' for s in strs)))
    open('%s/bdy_%s_w%d.m'%(outdir,reg,wt),'w').write('{'+', '.join(strs)+'}\n')
