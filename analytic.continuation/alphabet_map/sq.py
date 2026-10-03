import mpmath as mp, random
mp.mp.dps=40
def lets(x,y):
    z=1-x-y; s=mp.sqrt(x*z*y); I=mp.mpc(0,1)
    d={1:x,2:y,3:z,4:1-x,5:1-y,6:x+y,
       7:(x*z-I*s)/(x*z+I*s), 8:(x*y-I*s)/(x*y+I*s),
       9:1-(1+x)*y,10:1-x*(1+y),11:x-y+y**2,12:-x+x**2+y,13:(x-1)**2-y,
       14:-y+(x+y)**2,15:-x+(x+y)**2,16:-x+(1-y)**2,17:x+x*y+y**2,18:x**2+y+x*y,
       19:(y-1)**2+x*y,20:(x-1)**2+x*y,
       'X':(y*z-I*s)/(y*z+I*s)}   # third sqrt letter (a = yz), not in alphabet
    return d
def newl(u,w):           # original letters after x=-w/v, y=-u/v
    v=1-u-w; return lets(-w/v,-u/v)
def refl(u,w):           # reference letters, relabel x->w, y->u
    return lets(w,u)
def dlog(f,u,w,k):
    g=lambda a,b: mp.log(f(a,b)[k])
    return (mp.diff(lambda a: g(a,w),u), mp.diff(lambda b: g(u,b),w))
if __name__=="__main__":
    random.seed(1)
    pts=[(mp.mpf(random.uniform(.05,.45)),mp.mpf(random.uniform(.05,.45))) for _ in range(14)]
    for basis in (list(range(1,21)), list(range(1,21))+['X']):
        print("basis:", "alphabet" if len(basis)==20 else "alphabet + third sqrt letter")
        rows=[];  # each point gives 2 equations
        B=[[dlog(refl,u,w,j) for j in basis] for (u,w) in pts]
        for k in (7,8):
            A=[];b=[]
            for i,(u,w) in enumerate(pts):
                t=dlog(newl,u,w,k)
                for c in (0,1):
                    A.append([B[i][j][c] for j in range(len(basis))]); b.append(t[c])
            A=mp.matrix(A); b=mp.matrix(b)
            sol=mp.lu_solve(A.T*A, A.T*b) if True else None
            res=mp.norm(A*sol-b)
            coef={basis[j]:mp.nstr(mp.chop(sol[j],1e-20),6) for j in range(len(basis)) if abs(sol[j])>1e-15}
            print(f"  new l{k}: residual {mp.nstr(res,3)}  coeffs {coef}")
