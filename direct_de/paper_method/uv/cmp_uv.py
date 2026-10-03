import sys; sys.path.insert(0,'../..')
import mpmath as mpm
from exact_w2 import Exact, G_quad
mpm.mp.dps=30
def rd(fn):
    out=[];cur=None
    for l in open(fn):
        if l.startswith('#'): cur=[]; out.append((l.split()[1:3],cur))
        else: cur.append(mpm.mpc(*l.split()))
    return out
def dig(a,b):
    e=[abs(p-q)/abs(q) for p,q in zip(a,b) if abs(q)>1e-4]
    return float(-mpm.log10(max(e))) if e else 99
n=371
A=rd(sys.argv[1]); B=rd(sys.argv[2])
for (p,a),(_,b) in zip(A,B):
    u,v=map(mpm.mpf,p); x,y=-(1-u-v)/v,-u/v
    E=[Exact(w,30)(x,y,Gf=G_quad) for w in (1,2)]
    print('u,v=',p,' uv vs exact w1,w2: %.1f %.1f'%(dig(a[0:n],E[0]),dig(a[n:2*n],E[1])),
          '| uv-path vs xy-path w1..6:',' '.join('%.1f'%dig(a[w*n:(w+1)*n],b[w*n:(w+1)*n]) for w in range(6)), '| r',mpm.nstr(a[-1].real,10),mpm.nstr(b[-1].real,10))
