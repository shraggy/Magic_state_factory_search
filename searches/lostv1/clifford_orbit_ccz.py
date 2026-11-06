# BFS over the 3-qubit Clifford orbit of |CCZ>; check whether CS12CS13|+++> and others lie in it.
import numpy as np, itertools, sys
n=3; D=8
H=np.array([[1,1],[1,-1]])/np.sqrt(2); S=np.diag([1,1j]); I=np.eye(2)
def kron(*a):
    r=np.array([[1]]); 
    for x in a: r=np.kron(r,x)
    return r
def one(g,q): return kron(*[g if i==q else I for i in range(n)])
def cnot(c,t):
    M=np.zeros((D,D))
    for x in range(D):
        b=[(x>>(n-1-i))&1 for i in range(n)]
        if b[c]: b[t]^=1
        y=sum(b[i]<<(n-1-i) for i in range(n)); M[y,x]=1
    return M
gens=[one(H,q) for q in range(n)]+[one(S,q) for q in range(n)]+[cnot(a,b) for a in range(n) for b in range(n) if a!=b]
def canon(v):
    k=np.argmax(np.abs(v)>1e-9); v=v*np.conj(v[k])/abs(v[k])
    return tuple(np.round(v.real,6))+tuple(np.round(v.imag,6))
def diagstate(f):
    v=np.array([f([(x>>(n-1-i))&1 for i in range(n)]) for x in range(D)],dtype=complex); return v/np.linalg.norm(v)
ccz=diagstate(lambda b:(-1)**(b[0]*b[1]*b[2]))
targets={'CS12CS13':diagstate(lambda b:1j**(b[0]*b[1]+b[0]*b[2])),
 'CCZ123CS12':diagstate(lambda b:(-1)**(b[0]*b[1]*b[2])*1j**(b[0]*b[1])),
 'CS12 x |+>':diagstate(lambda b:1j**(b[0]*b[1])),
 'T1T2T3':diagstate(lambda b:np.exp(1j*np.pi/4*(b[0]+b[1]+b[2]))),
 'CS12CS13CS23':diagstate(lambda b:1j**(b[0]*b[1]+b[0]*b[2]+b[1]*b[2]))}
tk={k:canon(v) for k,v in targets.items()}
seen={canon(ccz)}; frontier=[ccz]
while frontier:
    nf=[]
    for v in frontier:
        for g in gens:
            w=g@v; c=canon(w)
            if c not in seen: seen.add(c); nf.append(w)
    frontier=nf
print('orbit size of |CCZ> under 3-qubit Cliffords:',len(seen))
for k,c in tk.items(): print(k,'in Clifford orbit of CCZ:',c in seen)
import cmath
r12=diagstate(lambda b: cmath.exp(1j*cmath.pi/4*(6*(b[0]+b[1]+b[2])+2*(b[0]*b[1]+b[0]*b[2]+b[1]*b[2])+4*b[0]*b[1]*b[2])))
print('[[12,3,2]] residual output of odd 5-qubit parent in Clifford orbit of CCZ:',canon(r12) in seen)
