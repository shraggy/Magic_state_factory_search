"""Exact number of classes of CCZ layers on k outputs, by Burnside's lemma over GL(k,2).

Up to CNOT circuits, Pauli X and diagonal Clifford gates on the outputs, a CCZ layer (a product of CCZ gates, phase
(-1)^f with f a cubic Boolean form) is a GL(k,2)-orbit of cubic forms modulo quadratic ones, RM(3,k)/RM(2,k), which is
the space Lambda^3 of alternating trilinear forms (X gates and diagonal Cliffords change only terms of degree <= 2).
Burnside's lemma, summed over the conjugacy classes of GL(k,2) (rational canonical forms from irreducible polynomials
and partitions, centralizer orders from the standard formula; the class sizes are checked to sum to |GL(k,2)|),
gives the number of orbits including the zero form. A form whose radical is nonzero depends on fewer than k
variables after a CNOT circuit, and those forms are exactly the orbits at k-1, so the classes that need all k outputs
number orb(k) - orb(k-1). Expected (counting zero): 2, 2, 3, 6, 12, 32 orbits for k = 3..8, so 1, 0, 1, 3, 6, 20
classes that need all k outputs. Known answers: k <= 6 agree with the genuine pure-CCZ classes of
outputs/class_table_l3.json (asserted); the nonzero counts 5, 11, 31 at k = 6, 7, 8 agree with Table 1 of
arXiv:2606.28473, which attributes them to X.-D. Hou (1996).
It also recomputes the GL(k,2)-invariant fingerprints of searches/k78_orbits.py on the target list
outputs/k78_pipeline/targets_all.json: distinct fingerprints are distinct orbits, so if their number equals the
Burnside count, the fingerprint separates every class and the list contains every class.
Deterministic, exact (fractions).
    python searches/k78_class_count.py   -> outputs/checks/k78_class_count.json   (about a minute)
"""
import os, sys, json, itertools, re, collections
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from fractions import Fraction
def pmul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1;b>>=1
    return r
def pmod(a,m):
    dm=m.bit_length()
    while a.bit_length()>=dm: a^=m<<(a.bit_length()-dm)
    return a
def irreducibles(d):
    out=[]
    for p in range(1<<d,1<<(d+1)):
        if all(pmod(p,q) for e in range(1,d//2+1) for q in range(1<<e,1<<(e+1))): out.append(p)
    return out
def partitions(n,mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    for i in range(min(n,mx),0,-1):
        for p in partitions(n-i,i): yield (i,)+p
def companion(p):
    d=p.bit_length()-1; M=[[0]*d for _ in range(d)]
    for i in range(1,d): M[i][i-1]=1
    for i in range(d): M[i][d-1]=(p>>i)&1
    return M
def blockdiag(bs):
    n=sum(len(b) for b in bs); M=[[0]*n for _ in range(n)]; o=0
    for b in bs:
        for i in range(len(b)):
            for j in range(len(b)): M[o+i][o+j]=b[i][j]
        o+=len(b)
    return M
def rank(rows):
    B=[]
    for v in rows:
        for b in B: v=min(v,v^b)
        if v: B.append(v)
    return len(B)
def fixdim(A,m):
    # action on Lambda^3: e_i^e_j^e_k -> Ae_i ^ Ae_j ^ Ae_k
    tr=list(itertools.combinations(range(m),3)); idx={t:i for i,t in enumerate(tr)}
    col=[[r for r in range(m) if A[r][c]] for c in range(m)]
    rows=[]
    for t in tr:
        v=0
        for a in col[t[0]]:
            for b in col[t[1]]:
                for c in col[t[2]]:
                    if len({a,b,c})==3: v^=1<<idx[tuple(sorted((a,b,c)))]
        v^=1<<idx[t]; rows.append(v)
    return len(tr)-rank(rows)
def phi(m,Q): 
    r=Fraction(1)
    for j in range(1,m+1): r*=1-Fraction(1,Q**j)
    return r
def cent(lam,Q):
    conj=[sum(1 for x in lam if x>i) for i in range(lam[0])] if lam else []
    r=Fraction(Q)**sum(c*c for c in conj)
    for v in set(lam): r*=phi(lam.count(v),Q)
    return r
def classes(m):
    irr=[(p,p.bit_length()-1) for d in range(1,m+1) for p in irreducibles(d) if p!=2]
    def rec(i,rem):
        if rem==0: yield []; return
        if i==len(irr): return
        p,d=irr[i]
        for s in range(0,rem//d+1):
            for lam in (partitions(s) if s else [()]):
                for rest in rec(i+1,rem-s*d): yield ([(p,d,lam)] if s else [])+rest
    yield from rec(0,m)
def G(m):
    r=1
    for i in range(m): r*=(2**m-2**i)
    return r



def orbits(m):
    tot = Fraction(0); S = Fraction(0); nc = 0
    for cl in classes(m):
        blocks = []; c = Fraction(1)
        for p, d, lam in cl:
            c *= cent(lam, 2 ** d)
            for part in lam:
                fp = 1
                for _ in range(part):
                    fp = pmul(fp, p)
                blocks.append(companion(fp))
        A = blockdiag(blocks); sz = Fraction(G(m)) / c
        tot += sz; S += sz * 2 ** fixdim(A, m); nc += 1
    assert tot == G(m), m
    S = S / G(m)
    assert S.denominator == 1
    return int(S), nc


def main():
    orb = {}; conj = {}
    for m in range(2, 9):
        orb[m], conj[m] = orbits(m) if m >= 3 else (1, 0)
    need_all = {m: orb[m] - orb[m - 1] for m in range(3, 9)}
    tab = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))
    known = collections.Counter(r["k"] for r in tab if all(g.startswith("CCZ") for g in r["spec"].split(",")))
    for m in range(3, 7):
        assert need_all[m] == known.get(m, 0), (m, need_all[m], known.get(m, 0))
    assert [orb[m] - 1 for m in (6, 7, 8)] == [5, 11, 31]
    import k78_orbits as K
    t = json.load(open(os.path.join(ROOT, "outputs", "k78_pipeline", "targets_all.json")))
    fp = collections.defaultdict(set); n = collections.Counter(); bad = 0; stored = collections.defaultdict(set)
    for r in t["targets"]:
        k = r["k"]; L = [tuple(x) for x in r["ccz"]]
        bad += not K.genuine(L, k); fp[k].add(K.fingerprint(L, k)); n[k] += 1; stored[k].add(r["fp"])
    lists = {k: dict(layers=n[k], non_genuine=0 if not bad else bad, fingerprints=len(fp[k]),
                     stored_fingerprint_labels=len(stored[k]), classes=need_all[k],
                     list_has_every_class=len(fp[k]) == need_all[k]) for k in sorted(n)}
    res = dict(orbits_including_zero={m: orb[m] for m in range(3, 9)}, conjugacy_classes={m: conj[m] for m in range(3, 9)},
               classes_needing_all_k_outputs=need_all, known_answer_k_le_6=dict(known), target_lists=lists)
    os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
    json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "k78_class_count.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
