import itertools
from math import comb
from d_mod import dist
def bi(N,W,l):
    return all(sum(comb(N-i,w-i) for w in W if w>=i)%(2**(l-i+1))==0 for i in range(1,l+1))
def pred(N,W,k):
    for wS in range(1,N-k+1):
        s=[wO for wO in range(0,k+1) if wO+wS in W]
        if sum(comb(k,w) for w in s)>1: return False
    return True
tot=agree=d3=0
for l in (2,3):
  for N in range(2,8):
    for r in range(1,N+1):
      for W in itertools.combinations(range(1,N+1),r):
        if not bi(N,W,l): continue
        for k in range(1,N):
          # target nonClifford? skip check; require at least one input
          d,n=dist(N,set(W),k,dmax=3)
          p=pred(N,set(W),k); real=(d!=1 and d!=2)
          tot+=1; agree+=(p==real); d3+=real
          if p!=real: print('MISMATCH',l,N,W,k,d)
          if real: print('d>=3',l,N,W,k,n,d)
print(agree,'of',tot,'agree;',d3,'with d>=3')
