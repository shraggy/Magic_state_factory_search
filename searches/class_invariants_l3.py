"""Class invariants of the 766 output classes at l=3 with k<=6: does a class contain an output of one gate type alone,
and does it contain a product of outputs on disjoint qubits?

Counted on the class, not on the gate string of its representative. A class is 'single-type' if its orbit (under
CNOT circuits, Pauli X and diagonal Cliffords on the outputs) contains an output built from one gate type alone (T only,
CS only or CCZ only); it is a 'product' if its orbit contains A (x) B with A, B genuine classes on k1 + k2 = k disjoint
qubits (complete, since the class of A (x) B depends only on the classes of A and B). Every candidate is mapped to its
exact class with the repository's orbit tables (searches/enumerate_targets.py for k<=5, enumerate_targets_fiber.py for
k=6). For T and CS every subset of supports is enumerated; for CCZ at k=6 every orbit in the fiber over the zero
quotient. Against the classes with fewer inputs than synthillation (n < syn_strongest_n in outputs/l3_catalogue_k6.csv)
it reports how many are entangled (not products) as classes and as labels.
Known-answer check (asserted): each representative maps back to its own class (1, 2, 6, 13, 65, 679 at k=1..6).
Deterministic; about 15 min on one core.
    nice -n 19 python searches/class_invariants_l3.py   -> outputs/checks/class_invariants_l3.json
"""
import os,sys,pickle,itertools,csv,collections,json
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0,HERE); sys.path.insert(0,os.path.join(ROOT,'classification'))
from enumerate_targets import enumerate_targets, index_of_poly
from coset_certify import parse_target
csv.field_size_limit(10**9)
lab={}; cls={}
for k in range(1,6):
    res,(C,off,L)=enumerate_targets(3,k,log=lambda *a:None)
    cls[k]={str(c['orbit']):c['spec'] for c in res['classes'] if c['genuine']}
    lab[k]=(lambda C,off,L,k: (lambda P: str(int(L[index_of_poly(P,k,3,C,off)]))))(C,off,L,k)
from enumerate_targets_fiber import enumerate_targets_fiber
res6,F=enumerate_targets_fiber(3,6,log=lambda *a:None)
print('fiber built',flush=True)
cls[6]={c['orbit']:c['spec'] for c in res6['classes']}
def l6(P):
    r=F.label(F.digits_of_poly(P)); return f"{r[0]}:{r[1]}"
lab[6]=l6
spec2orb={}
for k in cls:
    for o,s in cls[k].items(): spec2orb[(k,s)]=o
    assert all(lab[k](parse_target(s,3))==o for o,s in cls[k].items()), k
print('labels reproduce representatives', {k:len(cls[k]) for k in cls})
single={k:set() for k in cls}
for k in cls:
    for gtype,size,coef in (('T',1,1),('CS',2,2),('CCZ',3,4)):
        sups=list(itertools.combinations(range(k),size))
        if k==6 and size==3:
            import numpy as np
            r0,_=F.label(np.zeros(len(F.C),dtype=np.int64))
            for v in set(int(x) for x in np.unique(F.flab[r0])):
                o=f"{r0}:{v}"
                if o in cls[6]: single[6].add((o,'CCZ'))
            continue
        cap=len(sups)
        for r_ in range(1,cap+1):
          for sub in itertools.combinations(range(len(sups)),r_):
            P={frozenset(sups[i]):coef for i in sub}
            o=lab[k](P)
            if o in cls[k]: single[k].add((o,gtype))
    print(k,'single-type classes',len({o for o,_ in single[k]}))
prod={k:set() for k in cls}
def shift(P,s): return {frozenset(q+s for q in S):c for S,c in P.items()}
for k in range(2,7):
    for k1 in range(1,k//2+1):
        k2=k-k1
        for s1 in cls[k1].values():
            for s2 in cls[k2].values():
                P=dict(parse_target(s1,3)); P.update(shift(parse_target(s2,3),k1))
                o=lab[k](P); assert o in cls[k]; prod[k].add(o)
    print(k,'product classes',len(prod[k]),'of',len(cls[k]))
rows=list(csv.DictReader(open(os.path.join(ROOT,'outputs','l3_catalogue_k6.csv'))))
def types(s): return {g.rstrip('0123456789') for g in s.split(',')}
def conn(s):
    gs=[set(int(c) for c in g.lstrip('TCSZ')) for g in s.split(',')]
    qs=set().union(*gs); comp=[gs[0]]
    ch=True; cur=set(gs[0])
    while ch:
        ch=False
        for g in gs:
            if g&cur and not g<=cur: cur|=g; ch=True
    return cur==qs
out=[]
for r in rows:
    k=int(r['k']); o=spec2orb[(k,r['spec'])]
    st=sorted({t for oo,t in single[k] if oo==o})
    out.append(dict(k=k,spec=r['spec'],n=int(r['n']),syn=int(r['syn_strongest_n']),tau=int(r['tcount']),
      single=st,product=o in prod[k],label_types=len(types(r['spec'])),label_conn=conn(r['spec'])))
few=[x for x in out if x['n']<x['syn']]
print('rows',len(out),'fewer',len(few))
print('mixed by label',sum(x['label_types']>1 for x in out),'| classes with no single-type element',sum(not x['single'] for x in out))
print('label-mixed but class single-type:',[(x['spec'],x['single']) for x in out if x['label_types']>1 and x['single']])
print('label-pure classes also of other pure type:',[(x['spec'],x['single']) for x in out if len(x['single'])>1])
print('among 295: label-connected',sum(x['label_conn'] for x in few),'| class-entangled (not product)',sum(not x['product'] for x in few))
print('disagree:',[(x['spec'],x['label_conn'],x['product']) for x in out if x['label_conn']==x['product']])
print('T-only classes:',[(x['k'],x['spec'],x['n']) for x in out if 'T' in x['single']])
summ=dict(classes=len(out),fewer_than_syn=len(few),
  product_classes={k:len(prod[k]) for k in prod},genuine_classes={k:len(cls[k]) for k in cls},
  single_type_classes={k:len({o for o,_ in single[k]}) for k in single},
  single_type_by_gate=dict(collections.Counter(t for x in out for t in x['single'])),
  no_single_type_output=sum(not x['single'] for x in out),
  label_mixed=sum(x['label_types']>1 for x in out),
  fewer_entangled_class=sum(not x['product'] for x in few), fewer_product_class=sum(x['product'] for x in few),
  fewer_product_by_k=dict(collections.Counter(x['k'] for x in few if x['product'])),
  fewer_label_connected=sum(x['label_conn'] for x in few),
  t_only_classes=[(x['k'],x['spec'],x['n']) for x in out if 'T' in x['single']])
os.makedirs(os.path.join(ROOT,'outputs','checks'),exist_ok=True)
json.dump(dict(source='outputs/l3_catalogue_k6.csv',summary=summ,rows=out),open(os.path.join(ROOT,'outputs','checks','class_invariants_l3.json'),'w'),indent=0)
print(json.dumps(summ))
