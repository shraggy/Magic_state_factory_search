import itertools,random
random.seed(1);N=5;l=3;M=2**l
subs=[frozenset(s) for r in range(1,N+1) for s in itertools.combinations(range(N),r)]
def phi(c,Q): return sum(v for L,v in c.items() if Q<=L)%M
bad=0;tot=0
for t in range(200):
    c={L:random.randrange(M) for L in subs}
    p=list(range(N));random.shuffle(p)
    cp={frozenset(p[i] for i in L):v for L,v in c.items()}
    for Q in subs:
        if len(Q)<=l:
            tot+=1; bad+= phi(c,Q)!=phi(cp,frozenset(p[i] for i in Q))
print(bad,"of",tot,"mismatches")
