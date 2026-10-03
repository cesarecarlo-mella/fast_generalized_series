import sympy as sp, itertools
x,y,u,w=sp.symbols('x y u w')
v=1-u-w
L={1:x,2:y,3:1-x-y,4:1-x,5:1-y,6:x+y,9:1-(1+x)*y,10:1-x*(1+y),11:x-y+y**2,12:-x+x**2+y,
   13:(x-1)**2-y,14:-y+(x+y)**2,15:-x+(x+y)**2,16:-x+(1-y)**2,17:x+x*y+y**2,18:x**2+y+x*y,
   19:(y-1)**2+x*y,20:(x-1)**2+x*y}
sub={x:-w/v,y:-u/v}
def canon(p):
    p=sp.Poly(sp.expand(p),u,w); c=p.LC(); return sp.expand(p.as_expr()/c)
new={}
for k,p in L.items():
    e=sp.factor(sp.together(p.subs(sub)))
    n,d=sp.fraction(e)
    fs=[]
    for part,s in ((n,1),(d,-1)):
        c,fl=sp.factor_list(part)
        for f,m in fl:
            if f.free_symbols: fs.append((canon(f),s*m))
    new[k]=fs
    print(f"l{k}: {sp.factor(p)}  ->  {e}")
# original polys under relabelings (x,y,z)->perm of (w,u,v)
orig={k:canon(sp.expand(p)) for k,p in L.items()}
newset=set(f for fs in new.values() for f,_ in fs)
print("\nirreducible factors after map:",len(newset))
for perm in itertools.permutations([w,u,v]):
    X,Y,Z=perm
    mapped={k:canon(sp.expand(L[k].subs({x:X,y:Y},simultaneous=True))) for k in L}
    # also factor mapped (some may factor)
    mset=set()
    for k in L:
        c,fl=sp.factor_list(sp.expand(L[k].subs({x:X,y:Y},simultaneous=True)))
        for f,m in fl:
            if f.free_symbols: mset.add(canon(f))
    print("perm x,y,z ->",perm," new==orig set:",newset==mset," missing:",[str(f) for f in newset-mset]," unused:",len(mset-newset))
