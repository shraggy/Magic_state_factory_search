from math import comb
ok=0;tot=0
for N in range(2,14):
  for i in range(1,N):
    tot+=1; ok+= sum(comb(N-i,w-i) for w in range(2,N+1,2) if w>=i)==2**(N-i-1)
print('even-weight closed form holds',ok,'of',tot)
def bi(N,W,l): return all(sum(comb(N-i,w-i) for w in W if w>=i) % 2**(l-i+1)==0 for i in range(1,l+1))
for l in range(2,7):
  Nmin=min(N for N in range(2,20) if bi(N,range(2,N+1,2),l)); print('l',l,'min N even',Nmin,'rotations',sum(comb(Nmin,w) for w in range(2,Nmin+1,2)))
