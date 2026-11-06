# Exact phase-table simulation (units of pi/4, mod 8) of v1 Fig.6 Quirk parents and the current odd-weight parent.
import itertools
def rank(vs):
    vs=[v for v in vs]; r=0; rows=[]
    for v in vs:
        for b in rows: v=min(v,v^b)
        if v: rows.append(v); r+=1
    return r
def factory(rots,N,checks):
    outs=[q for q in range(N) if q not in checks]
    keep=[R for R in rots if set(R)&set(checks)]
    rem=[R for R in rots if not set(R)&set(checks)]
    # postselect checks on |+>: amplitude on output basis y = sum_xc w^{phase(keep)}; check |.| const and residual phase
    amps={}
    for y in itertools.product((0,1),repeat=len(outs)):
        s=0
        import cmath
        for xc in itertools.product((0,1),repeat=len(checks)):
            x=dict(zip(checks,xc)); x.update(zip(outs,y))
            ph=sum(1 for R in keep if sum(x[q] for q in R)%2)
            s+=cmath.exp(1j*cmath.pi/4*ph)
        amps[y]=s
    mags=set(round(abs(a),6) for a in amps.values())
    # residual phase in units pi/4
    import math
    res={y:round((cmath.phase(a)/(math.pi/4)))%8 for y,a in amps.items()} if len(mags)==1 else None
    # output parts rank (affine) of kept rotations
    def mask(R): return sum(1<<outs.index(q) for q in R if q in outs)
    ms=[mask(R) for R in keep]
    aff=rank([m^ms[0] for m in ms[1:]])
    return len(keep),len(outs),mags,res,aff,rem
def anf(res,k):
    # Moebius transform of residual phase (units pi/4 mod 8)
    a={}
    for r in range(k+1):
        for T in itertools.combinations(range(k),r):
            v=0
            for rr in range(r+1):
                for U in itertools.combinations(T,rr):
                    y=tuple(1 if i in U else 0 for i in range(k)); v+=(-1)**(r-rr)*res[y]
            a[T]=v%8
    return {T:v for T,v in a.items() if v}
# v1 six-qubit parent (link 39): qubits 0..5
six=[(0,1,3,4,5),(0,2,3,4,5),(0,1,5),(0,2,5),(0,1,4),(0,2,4),(0,1,3),(0,2,3),(1,2,3,4,5),(1,2,5),(1,2,4),(1,2,3),(3,4,5)]
for checks in ([0],[0,1,2]):
    n,k,m,res,aff,rem=factory(six,6,checks)
    print('six-qubit parent cut',checks,'n=',n,'k=',k,'|amp| set',m,'affine rank of output parts',aff,'removed',rem)
    if res: print('   residual ANF (units pi/4):',anf(res,k))
# BI check of full six-qubit parent
import cmath
ok=len(set(sum(1 for R in six if sum(x[q] for q in R)%2)%8 for x in itertools.product((0,1),repeat=6)))
print('six-qubit parent: distinct phases over basis states (BI iff 1):',ok)
# current odd-weight five-qubit parent
five=[S for w in (5,3,1) for S in itertools.combinations(range(5),w)]
print('odd five-qubit parent: distinct phases:',len(set(sum(1 for R in five if sum(x[q] for q in R)%2)%8 for x in itertools.product((0,1),repeat=5))))
for checks in ([0],[0,1],[0,1,2]):
    n,k,m,res,aff,rem=factory(five,5,checks)
    print('odd parent cut',checks,'n=',n,'k=',k,'|amp|',m,'aff rank',aff)
    if res: print('   residual ANF:',anf(res,k))
from collections import Counter
print(Counter(sum(1 for R in six if sum(x[q] for q in R)%2)%8 for x in itertools.product((0,1),repeat=6)))
