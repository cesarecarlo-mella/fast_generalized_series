import mpmath as mp, random
from sq import lets, newl, refl
mp.mp.dps=30
M={1:{1:1,3:-1},2:{2:1,3:-1},3:{3:-1},4:{5:1,3:-1},5:{4:1,3:-1},6:{6:1,3:-1},
   7:{7:-1,8:-1},8:{8:1},
   9:{13:1,3:-2},10:{16:1,3:-2},11:{12:1,3:-2},12:{11:1,3:-2},13:{9:1,3:-2},14:{18:1,3:-2},
   15:{17:1,3:-2},16:{10:1,3:-2},17:{15:1,3:-2},18:{14:1,3:-2},19:{20:1,3:-2},20:{19:1,3:-2}}
random.seed(7)
pts=[(mp.mpf(random.uniform(.02,.6)),mp.mpf(random.uniform(.02,.6))) for _ in range(400)]
pts=[p for p in pts if p[0]+p[1]<0.97][:200]
for k in range(1,21):
    offs=set()
    for u,w in pts:
        N=newl(u,w);R=refl(u,w)
        d=mp.log(N[k])-sum(c*mp.log(R[j]) for j,c in M[k].items())
        offs.add((mp.nstr(mp.chop(mp.re(d),1e-20),6), mp.nstr(mp.chop(mp.im(d)/mp.pi,1e-20),6)))
    print(f"l{k:<2} log(new) - combo = {sorted(offs)}  (re, im/pi) over {len(pts)} pts")
# also: third sqrt letter X in terms of l7, l8 ?
offs=set()
for u,w in pts:
    R=refl(u,w); d=mp.log(R['X'])+mp.log(R[7])+mp.log(R[8])
    offs.add((mp.nstr(mp.chop(mp.re(d),1e-20),6), mp.nstr(mp.chop(mp.im(d)/mp.pi,1e-20),6)))
print("X + l7 + l8 =",sorted(offs))
