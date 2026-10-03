import sympy as sp
exec(open('rat.py').read().split("orig=")[0].replace('print(f"l{k}','pass#'))
# reference letters in new variables: relabel x->w, y->u (z->v)
ref={}
for k,p in L.items():
    c,fl=sp.factor_list(sp.expand(p.subs({x:w,y:u},simultaneous=True)))
    assert len(fl)==1 and fl[0][1]==1
    ref[canon(fl[0][0])]=k
def conic(k):
    p=sp.Poly(L[k],x,y)
    A=p.coeff_monomial(x**2);B=p.coeff_monomial(x*y);C=p.coeff_monomial(y**2)
    if A==B==C==0: return "lin"
    d=B**2-4*A*C; return "parab" if d==0 else ("hyperb" if d>0 else "ellip")
for k in L:
    terms=" ".join(("%+d"%m if m!=1 else "+")+f"l{ref[f]}" for f,m in new[k])
    print(f"l{k:<2} [{conic(k):6}] -> {terms.lstrip('+')}" + ("" if len(new[k])>1 or True else ""))
