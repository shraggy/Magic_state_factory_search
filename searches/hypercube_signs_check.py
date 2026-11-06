"""Lemma (hypercube factories), exact check of its proof at l = 2..5, in exact rational arithmetic.

(1) The free lattice of a degree-t term at level l is 2^(e_t - l) 2 pi, e_t = min(l, max(l-3+t, t)) (App. S1.1), which
    contains fewer free gates than 2 pi / 2^(l-t) for l >= 4 (a single T would be free at l = 4, 5 under the latter).
(2) Signs that meet the bound tau(C^{l-1}Z) = 2^l - 1: c_Lambda = (-1)^(|Lambda|+1) on every nonempty subset of the l
    outputs gives C^{l-1}Z exactly at every l; all coefficients +1 give it up to a free diagonal gate at l = 2, 3 only.
Prints one line per level and writes outputs/checks/hypercube_signs_check.json. Asserted: alternating signs exact at
l = 2..5; all +1 free at l = 2, 3 and not free at l = 4, 5; T not free under the true lattice at any l.
    python searches/hypercube_signs_check.py
"""
import os, json, itertools
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
# phases in units of 2pi, exact fractions
def monos(l, phase):  # phase: dict x-tuple -> Fraction (units of 2pi); return Mobius coefficients
    co={}
    for Q in itertools.product([0,1],repeat=l):
        s=F(0)
        for x in itertools.product([0,1],repeat=l):
            if all(x[i]<=Q[i] for i in range(l)):
                s+=(-1)**(sum(Q)-sum(x))*phase[x]
        co[Q]=s
    return co
def e(l,t): return min(l,max(l-3+t,t))
def is_free(l,co,lattice):
    for Q,v in co.items():
        t=sum(Q)
        if t==0: continue
        g=lattice(l,t)
        if (v/g).denominator!=1: return False
    return True
true_lat=lambda l,t: F(2**e(l,t),2**l)          # 2^{e_t-l} (units of 2pi)
paper_lat=lambda l,t: F(1,2**(l-t))             # 2pi/2^{l-|Q|}
dclif_lat=lambda l,t: {1:F(1,4),2:F(1,2)}.get(t,F(1))  # diagonal Clifford pi/2, pi, 2pi
def prep(l,c):
    subs=[L for L in itertools.product([0,1],repeat=l) if any(L)]
    ph={}
    for x in itertools.product([0,1],repeat=l):
        ph[x]=sum(F(c[L],2**l)*(sum(a*b for a,b in zip(L,x))%2) for L in subs)
    return subs,ph
res = {}
for l in [2,3,4,5]:
    print('l',l,'true lattice',[str(true_lat(l,t)) for t in range(1,l+1)],'old',[str(paper_lat(l,t)) for t in range(1,l+1)],'diagClif',[str(dclif_lat(l,t)) for t in range(1,l+1)])
    subs=[L for L in itertools.product([0,1],repeat=l) if any(L)]
    tgt={x:F(1,2) if all(x) else F(0) for x in itertools.product([0,1],repeat=l)}
    for name,c in [('all+1',{L:1 for L in subs}),('alt',{L:(-1)**(sum(L)+1) for L in subs})]:
        _,ph=prep(l,c)
        diff={x:ph[x]-tgt[x] for x in ph}
        co=monos(l,diff)
        r = dict(free_true=is_free(l,co,true_lat), free_old_lattice=is_free(l,co,paper_lat), exact=all(v.denominator==1 for v in diff.values()))
        res.setdefault(l, {})[name] = r
        print('  ',name,r)
    # T on qubit 1 under paper lattice at l=4
    T={x:F(x[0],8) for x in itertools.product([0,1],repeat=l)}
    res[l]['T_free_old_lattice'] = is_free(l,monos(l,T),paper_lat); res[l]['T_free_true'] = is_free(l,monos(l,T),true_lat)
    print('   T free under old lattice:', res[l]['T_free_old_lattice'], ' under true:', res[l]['T_free_true'])
for l in res:
    assert res[l]['alt']['exact'] and not res[l]['T_free_true']
    assert res[l]['all+1']['free_true'] == (l <= 3)
os.makedirs(os.path.join(ROOT, 'outputs', 'checks'), exist_ok=True)
json.dump({str(l): v for l, v in res.items()}, open(os.path.join(ROOT, 'outputs', 'checks', 'hypercube_signs_check.json'), 'w'), indent=1)
