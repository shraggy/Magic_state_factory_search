"""Every l=3 output class with k <= 6 (766 classes): our smallest factory against synthillation (Fig. 1).

For each class in outputs/cascade/l3_k{1..6}.json (best factory per class: n, N, stage, target spec) this computes
  tcount  exact T-count of the target (Amy-Mosca decoder), so n - tcount - 1 >= 0 is the gap to the bound of App. S4;
  syn_n   inputs of the synthillation construction (Campbell-Howard, arXiv:1606.01904, Thm 1 + Table I), built with
          searches/synthillation_compare.py (optimal A, smallest Lempel B, mu and mu+1) and checked with the same
          borrowed-identity test; None when no valid construction results (the Table I case does not apply);
  syn_N   its qubit count k + rows of the Table I case.
    python searches/class_table_l3.py --workers 8      -> outputs/class_table_l3.json
"""
import os, sys, json, argparse
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def parse(spec):
    P = {}
    for g in spec.split(","):
        if g.startswith("CCZ"):
            P[frozenset(int(c) for c in g[3:])] = 4
        elif g.startswith("CS"):
            P[frozenset(int(c) for c in g[2:])] = 2
        elif g.startswith("T"):
            P[frozenset([int(g[1:])])] = 1
    return P


def synth_n(k, target):
    import synthillation_compare as sc
    A = sc.opt_A(k, target); tau = len(A)
    lin = any(len(T) == 1 for T in target)
    mu0 = 0
    while sc.opt_B(k, target, mu0) is None:
        mu0 += 1
    best = None
    for mu in (mu0, mu0 + 1):
        B = sc.opt_B(k, target, mu)
        if B is None:
            continue
        if lin:
            tm = sc.C14[{(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(tau % 2, mu % 2)]]
            c = sum(1 << i for i in range(k) if frozenset([i]) in target); blocks = [A, B, B] + [[c]] * 4
        else:
            tm = sc.T5[(tau % 2, mu % 2)]; blocks = [A, B, B]
        blocks += [[0]] * (len(tm[0]) - len(blocks)); R = len(tm)
        cols = [v | sum(1 << (k + r) for r in range(R) if tm[r][bi]) for bi, blk in enumerate(blocks) for v in blk]
        if sc.coeffs_ok(cols, k + R, k, target, 3) and all(y >> k for y in cols):
            if best is None or len(cols) < best[0]:
                best = (len(cols), k + R, mu)
    return tau, mu0, best


def one(item):
    k, orbit, b = item
    P = parse(b["spec"])
    tau, mu0, best = synth_n(k, P)
    return dict(k=k, orbit=orbit, spec=b["spec"], n=b["n"], N=b["N"], stage=b["stage"], tcount=tau, mu=mu0,
                gap=b["n"] - tau - 1, syn_n=best and best[0], syn_N=best and best[1], syn_mu=best and best[2],
                source=b.get("source", "this repo"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=1); ap.add_argument("--kmax", type=int, default=6)
    a = ap.parse_args()
    items = []
    # one row per l=3 class; tcount, mu and synthillation depend only on the class, n and N are the staged cascade's
    for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]:
        if r["l"] == 3 and r["k"] <= a.kmax:
            items.append((r["k"], r["orbit"], dict(spec=r["spec"], n=r["n"], N=r["N"], stage=r["stage"])))
    with Pool(a.workers) as p:
        rows = []
        for i, r in enumerate(p.imap_unordered(one, items, chunksize=1)):
            rows.append(r)
            if i % 50 == 0:
                print(i, len(items), flush=True)
    rows.sort(key=lambda r: (r["k"], r["n"], r["spec"]))
    json.dump(rows, open(os.path.join(ROOT, "outputs", "class_table_l3.json"), "w"), indent=0)
    print("done", len(rows), flush=True)
