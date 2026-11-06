import itertools, numpy as np
th=np.pi/8
def ph(L,x): return 2*th*(sum(x[q] for q in L)%2)
kept=[(0,)+S for r in range(4) for S in itertools.combinations((1,2,3),r)]
rem=[S for r in range(1,4) for S in itertools.combinations((1,2,3),r)]
# full BI check
xs=list(itertools.product((0,1),repeat=4))
full=[sum(ph(L,x) for L in kept+rem)%(2*np.pi) for x in xs]
print("full BI phases distinct:",sorted(set(np.round(full,9))))
# output after projecting q0 on |+>
out={}
for x in itertools.product((0,1),repeat=3):
    a=sum(np.exp(1j*sum(ph(L,(x0,)+x) for L in kept)) for x0 in (0,1))/2
    out[x]=a
a0=out[(0,0,0)]
for x,a in out.items(): print(x, "|amp|=%.3f"%abs(a), "phase/theta=%.3f"%(np.angle(a/a0)/th))
# compare to CCZ|+> up to Clifford: target = -pi on x!=0 = Z1Z2Z3 CZ12CZ13CZ23 CCZ
for x in out:
    w=sum(x); p=np.pi*(w - (x[0]*x[1]+x[0]*x[2]+x[1]*x[2]) + x[0]*x[1]*x[2])
    assert abs(np.exp(1j*p)-out[x]/a0)<1e-9
print("output = Z^3 CZ^3 CCZ |+>^3 : OK; 8 of 8 basis states")
# accept prob
print("acceptance norm^2:",sum(abs(a)**2 for a in out.values())/8)
