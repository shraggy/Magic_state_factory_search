"""Acceptance test of the XOR free-SAT pass against sweep A's free SAT (iii_f), on every l=3 class with k <= 6.

Both passes start from the same seeds: the published iii_r best per class (outputs/staged_pc/l3_k{1..4}/iii_r.json,
outputs/staged_A/l3_k{5,6}/iii_r.json). Per class:
  A_n    = n after sweep A's iii_f (iii_f.json best: min of seed, CP-SAT bool:w5:d600 and z3 coef:1.5e9:unseeded)
  X_n    = n after the XOR pass = min(seed, every factory the chosen XOR variants credit to the class)
  A_cpu  = sweep A free-SAT CPU seconds on the class (bool + coef, per-class `seconds`)
  X_cpu  = XOR pass CPU seconds on the class (sum over the chosen variants)
Counts better / equal / worse (X_n vs A_n). Also reported per single variant.
    python searches/xor_acceptance.py OUT [--vars x_w5_c2e6,x_wNone_c2e6]  -> OUT/acceptance_iii_f.json
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CELLS = [(1, "staged_pc"), (2, "staged_pc"), (3, "staged_pc"), (4, "staged_pc"), (5, "staged_A"), (6, "staged_A")]


def secs(pat):
    s = {}
    for f in glob.glob(pat):
        for o, p in json.load(open(f))["per_class"].items():
            s[o] = s.get(o, 0.0) + p.get("seconds", 0.0)
    return s


def xor_of(out, k):
    """{variant: (own best n per class, seconds per class, classes done)}"""
    res = {}
    for vd in sorted(glob.glob(os.path.join(out, f"l3_k{k}", "xor_x_*"))):
        v = os.path.basename(vd)[4:]
        if "seed" in v:
            continue                                     # gluing-seeded runs are a separate comparison
        own, sec, done = {}, {}, set()
        for f in glob.glob(os.path.join(vd, "*.json")):
            x = json.load(open(f)); done |= set(x["classes"])
            for o, p in x["per_class"].items():
                sec[o] = p.get("seconds", 0.0)
            for r in x["results"]:
                o = str(r["orbit"])
                if o not in own or r["n"] < own[o]:
                    own[o] = r["n"]
        res[v] = (own, sec, done)
    return res


def main(out, chosen=None):
    rows, per_var = [], {}
    for k, src in CELLS:
        d = os.path.join(ROOT, "outputs", src, f"l3_k{k}")
        seed = json.load(open(os.path.join(d, "iii_r.json")))
        A = json.load(open(os.path.join(d, "iii_f.json")))["best"]
        acpu = secs(os.path.join(d, "free_*", "*.json"))
        X = xor_of(out, k)
        use = [v for v in X if chosen is None or v in chosen]
        for o in sorted(seed["specs"]):
            s = seed["best"].get(o, {}).get("n")
            xs = [X[v][0][o] for v in use if o in X[v][0]]
            xn = min([n for n in [s] + xs if n is not None], default=None)
            done = all(o in X[v][2] for v in use)
            an = A.get(o, {}).get("n")
            cmp = None if xn is None or an is None else ("better" if xn < an else "equal" if xn == an else "worse")
            rows.append(dict(k=k, orbit=o, spec=seed["specs"][o], seed=s, A_n=an, X_n=xn, result=cmp, done=done,
                             A_cpu=round(acpu.get(o, 0.0), 1), X_cpu=round(sum(X[v][1].get(o, 0.0) for v in use), 1)))
            for v in X:
                xv = min([n for n in [s, X[v][0].get(o)] if n is not None], default=None)
                c = per_var.setdefault(v, dict(better=0, equal=0, worse=0, not_done=0, cpu=0.0))
                if o not in X[v][2]:
                    c["not_done"] += 1; continue
                c["cpu"] += X[v][1].get(o, 0.0)
                c["better" if xv < an else "equal" if xv == an else "worse"] += 1
    cnt = {r: sum(1 for x in rows if x["result"] == r) for r in ("better", "equal", "worse")}
    cnt["undecided"] = sum(1 for x in rows if x["result"] is None); cnt["not_done"] = sum(1 for x in rows if not x["done"])
    rec = dict(chosen=sorted(chosen) if chosen else "all xor variants", classes=len(rows), counts=cnt,
               A_cpu_seconds=round(sum(x["A_cpu"] for x in rows), 1), X_cpu_seconds=round(sum(x["X_cpu"] for x in rows), 1),
               per_variant={v: dict(c, cpu=round(c["cpu"], 1)) for v, c in per_var.items()},
               worse=[x for x in rows if x["result"] == "worse"], better=[x for x in rows if x["result"] == "better"],
               rows=rows)
    json.dump(rec, open(os.path.join(out, "acceptance_iii_f.json"), "w"), indent=0)
    print(f"{len(rows)} classes; {cnt}; CPU-s sweep A free SAT {rec['A_cpu_seconds']}, XOR {rec['X_cpu_seconds']}")
    for v, c in rec["per_variant"].items():
        print(" ", v, c)
    for x in rec["worse"]:
        print("  worse:", x["k"], x["orbit"], "A", x["A_n"], "X", x["X_n"])


if __name__ == "__main__":
    a = sys.argv[1:]; ch = None
    if "--vars" in a:
        i = a.index("--vars"); ch = set(a[i + 1].split(",")); a = a[:i] + a[i + 2:]
    main(a[0], ch)
