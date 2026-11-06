# [[8,3,2]] cube factory vs Clifford-padded [[8,4,2]]; Z-fault model on the 8 T inputs, exact statevector.
import numpy as np, itertools
def run(N, outs, rots, Mfun):
    # rots: list of (support set, coeff in units of pi/4 applied as exp(i*pi/4*c*parity))
    dim=1<<N; xs=np.arange(dim)
    bit=lambda x,q:(x>>q)&1
    par=lambda S:np.bitwise_xor.reduce([bit(xs,q) for q in S]) if S else 0*xs
    checks=[q for q in range(N) if q not in outs]
    def final(faults):
        ph=np.zeros(dim)
        for i,(S,c) in enumerate(rots): ph+= np.pi/4*c*par(S)
        psi=np.exp(1j*ph)/np.sqrt(dim)
        for i in faults:
            psi=psi*(1-2*par(rots[i][0]))
        # project checks on |+>
        k=len(outs); out=np.zeros(1<<k,complex)
        for x in range(dim):
            y=sum(bit(x,q)<<j for j,q in enumerate(outs)); out[y]+=psi[x]
        out/=np.sqrt(1<<len(checks))
        return out
    ideal=final([]); p0=np.vdot(ideal,ideal).real; ideal/=np.sqrt(p0)
    n=len(rots); res={}
    for w in (1,2):
        bad=0; acc=0
        for F in itertools.combinations(range(n),w):
            o=final(list(F)); pa=np.vdot(o,o).real/p0
            if pa>1e-9:
                acc+=1
                f=abs(np.vdot(ideal,o/np.sqrt(np.vdot(o,o).real)))**2
                if f<1-1e-9: bad+=1
        res[w]=(acc,bad,len(list(itertools.combinations(range(n),w))))
    return p0,ideal,res
# cube: check qubit 3, outputs 0,1,2 ; rotation on parity x3^x_T, coefficient (-1)^{|T|} (units pi/4 -> T gates)
rots=[(set(T)|{3},(-1)**len(T)) for r in range(4) for T in itertools.combinations(range(3),r)]
p0,ideal,res=run(4,[0,1,2],rots,None)
# check output is CCZ up to diagonal Clifford: phases
ph=np.angle(ideal/ideal[0])/(np.pi/4); print("cube p0",round(p0,4),"output phases (pi/4)",np.round(ph,3)); print("cube faults (accepted, accepted-and-wrong, total) w1,w2:",res)
# padded: add output qubit 4 in |+>, apply CNOT circuit M on outputs: map support S -> image. Use x-substitution:
# target CCZ123 CCZ134 CCZ234 (0-indexed outputs 0,1,2,3 with qubit 4 being check). Build by CNOT on output parities:
# rotation parity over (x_check + a.x_out); CNOT circuit acts as x_out -> A x_out ; parity a.x -> (A^T a).x
A=np.array([[1,0,0,1],[0,1,0,1],[0,0,1,0],[0,0,0,1]])  # new vars: u0=x0+x3,u1=x1+x3,u2=x2,u3=x3
rots2=[]
for S,c in rots:
    a=np.zeros(4,int)
    for q in S:
        if q<3: a[q]=1
    b=(A.T@a)%2
    rots2.append((set(np.nonzero(b)[0].tolist())|{4},c))
p0b,ideal2,res2=run(5,[0,1,2,3],rots2,None)
print("padded rotations",rots2); print("padded p0",round(p0b,4),"faults:",res2)
ph2=np.angle(ideal2/ideal2[0])/(np.pi/4); print("padded phases",np.round(ph2,3))
# can the padded output detect anything more? stabilizer of the |+> factor in new frame = X on qubits mapped; test: for every weight-2 accepted wrong fault, is the
# extra-qubit stabilizer measurement (in original frame: X on padded qubit) ever flipped? equivalently compare res with res2.
# kernel direction v of the padded target: phase invariant under x->x^v ?
for v in range(1,16):
    if all(abs(ph2[x]-ph2[x^v])<1e-9 or abs(abs(ph2[x]-ph2[x^v])-8)<1e-9 for x in range(16)):
        vv=[(v>>j)&1 for j in range(4)]
        print("free X direction v",vv,"all rotation output parts orthogonal:",all(sum(vv[q] for q in S if q<4)%2==0 for S,c in rots2))
