"""Independent re-verification of every per-class best factory of a k78 pass file (us.is_factory monomial and brute
force, us.distance == 2, entangled_gluing.deposits against the class's layer, n recomputed), and the list of classes
not on the bound written as XOR targets.   python searches/k78_verify_pass.py PASS.json [XOR_TARGETS_OUT.json]
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import unified_search as us
from entangled_gluing import deposits
d = json.load(open(sys.argv[1])); k = d["k"]; ok = 0; todo = []
for o, b in sorted(d["best"].items()):
    g = [(frozenset(S), c) for S, c in b["gates"]]; f = dict(k=k, r=b["N"] - k, gates=g)
    v = bool(us.is_factory(f, 3) and us.is_factory(f, 3, brute_force=True) and us.distance(f, 3) == 2
             and deposits(g, k, {frozenset(x): 4 for x in b["ccz"]}) and sum(us.inputs(c, 3) for _, c in g) == b["n"])
    ok += v
    if not v:
        print("FAILED", o)
    if not b.get("class_on_bound"):
        todo.append(dict(id=b.get("id") or ("B:" + o), k=k, ccz=b["ccz"], fp=o))
print(f"k={k}: {ok} of {len(d['best'])} class-best factories re-verified; {len(todo)} classes not on the bound")
if len(sys.argv) > 2:
    old = json.load(open(sys.argv[2])) if os.path.exists(sys.argv[2]) else []
    json.dump(old + todo, open(sys.argv[2], "w"))
