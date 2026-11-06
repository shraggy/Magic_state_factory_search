import itertools,sys
def dist(N,W,k=1,dmax=4):
    O=set(range(k))
    rots=[frozenset(c) for w in W for c in itertools.combinations(range(N),w)]
    inp=[r for r in rots if not r<=O]
    vec=[(sum(1<<q for q in r if q>=k), sum(1<<q for q in r if q<k)) for r in inp]
    for d in range(1,dmax+1):
        for S in itertools.combinations(range(len(vec)),d):
            c=o=0
            for i in S: c^=vec[i][0]; o^=vec[i][1]
            if c==0 and o!=0: return d,len(inp)
    return '>%d'%dmax,len(inp)
if __name__=='__main__':
  for l in (2,3):
    N=l+1; print('l',l,'all',N,dist(N,range(1,N+1)))
    N=l+2; print('l',l,'odd',N,dist(N,range(1,N+1,2)))
