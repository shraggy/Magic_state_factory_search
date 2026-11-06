import itertools
def run(N,O):
    subs=[frozenset(L) for r in range(1,N+1) for L in itertools.combinations(range(N),r)]
    Qs=[frozenset(Q) for r in range(1,4) for Q in itertools.combinations(range(N),r) if not set(Q)<=O]
    best={}
    for m in range(1,1<<len(subs)):
        S=[subs[j] for j in range(len(subs)) if m>>j&1]
        if any(sum(1 for L in S if Q<=L)%2 for Q in Qs): continue
        out=tuple(sorted(tuple(sorted(L)) for L in S if L<=O))
        n=sum(1 for L in S if not L<=O)
        if n==0: continue
        best[out]=min(best.get(out,99),n)
    return best
for N in (3,4):
    for k in (1,2):
        print(N,k,run(N,frozenset(range(k))))
