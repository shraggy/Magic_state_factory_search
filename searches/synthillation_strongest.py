"""Strongest synthillation construction per l=3 class (Campbell-Howard arXiv:1606.01904, Thm 1 and Table I).

Thm 1 states n = col(A) + 2 col(B) + Delta, Delta fixed by the parities of col(A), col(B) and whether the target has
linear terms (Table I, cases 1-11). The construction holds for any synthesis matrices A, B, not only optimal ones
("if suboptimal matrices are used, the resource cost is n = col(A) + 2col(B) + Delta"), so a zero column appended to
A or B (which changes nothing in |A^T x| or |B^T x|) flips a parity and can lower Delta or make the template valid.
For each class we build every variant A in {A*, A*+0}, B in {B*, B*+0, B(mu*+1)} with A* optimal (tau columns) and
B* a smallest Lempel factor (mu* columns), on its Table I case (cases 9/11 with one check row for pure-CCZ targets),
check each with the borrowed-identity condition, and keep the smallest valid n.
    python searches/synthillation_strongest.py --workers 8   -> outputs/synthillation_strongest_l3.json
"""
import os, sys, json, argparse
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def build(k, target, A, B, lin, c):
    import synthillation_compare as sc
    tau, mu = len(A), len(B)
    if lin:
        tm = sc.C14[{(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(tau % 2, mu % 2)]]
        blocks = [A, B, B] + [[c]] * 4
    elif mu == 0:                                         # cases 9 (tau even) and 11 (tau odd): one check row
        cols = [a | 1 << k for a in A] + ([1 << k] if tau % 2 else [])
        ok = sc.coeffs_ok(cols, k + 1, k, target, 3) and all(y >> k for y in cols)
        return (len(cols), k + 1, "9" if tau % 2 == 0 else "11") if ok else None
    else:
        tm = sc.T5[(tau % 2, mu % 2)]; blocks = [A, B, B]
    blocks += [[0]] * (len(tm[0]) - len(blocks)); R = len(tm)
    cols = [v | sum(1 << (k + r) for r in range(R) if tm[r][bi]) for bi, blk in enumerate(blocks) for v in blk]
    ok = sc.coeffs_ok(cols, k + R, k, target, 3) and all(y >> k for y in cols)
    case = ({(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4} if lin else {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8})[(tau % 2, mu % 2)]
    return (len(cols), k + R, str(case)) if ok else None


def one(r):
    import class_table_l3 as CT, synthillation_compare as sc
    k = r["k"]; P = CT.parse(r["spec"]); A = sc.opt_A(k, P)
    mu = 0
    while sc.opt_B(k, P, mu) is None:
        mu += 1
    B0 = sc.opt_B(k, P, mu); B1 = sc.opt_B(k, P, mu + 1)
    lin = any(len(T) == 1 for T in P); c = sum(1 << i for i in range(k) if frozenset([i]) in P)
    variants = []
    for an, Av in (("A", A), ("A+0", A + [0])):
        for bn, Bv in (("B", B0), ("B+0", B0 + [0]), ("B(mu+1)", B1)):
            if Bv is None:
                continue
            if not lin and len(Bv) > 0 and mu == 0 and bn != "B":
                pass                                      # pure CCZ with a padded B: generic case 5-8 template
            res = build(k, P, Av, Bv, lin, c)
            variants.append(dict(A=an, B=bn, tau=len(Av), mu=len(Bv), valid=res is not None,
                                 n=res and res[0], N=res and res[1], case=res and res[2]))
    valid = [v for v in variants if v["valid"]]
    best = min(valid, key=lambda v: (v["n"], v["N"])) if valid else None
    thm1 = [v for v in variants if v["A"] == "A" and v["B"] == "B"][0]
    return dict(k=k, spec=r["spec"], tau=len(A), mu=mu, thm1_n=thm1["n"] if thm1["valid"] else None,
                thm1_formula=len(A) + 2 * mu, best_n=best and best["n"], best_N=best and best["N"],
                best_variant=best and f'{best["A"]},{best["B"]},case {best["case"]}', variants=variants)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=1); a = ap.parse_args()
    R = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))
    with Pool(a.workers) as p:
        out = list(p.imap(one, R, chunksize=1))
    json.dump(out, open(os.path.join(ROOT, "outputs", "synthillation_strongest_l3.json"), "w"))
    M = {(x["k"], x["spec"]): x for x in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]
         if x["l"] == 3}                                  # smallest factory per class, staged cascade
    fewer = eq = more = fails = 0
    for o in out:
        n = M[(o["k"], o["spec"])]["n"]
        if o["best_n"] is None:
            fails += 1
        elif n < o["best_n"]:
            fewer += 1
        elif n == o["best_n"]:
            eq += 1
        else:
            more += 1
    print(f"vs strongest synthillation variant: ours fewer {fewer}, equal {eq}, more {more}, no valid variant {fails}")
    print("classes where a padded variant beats Thm 1 as stated:",
          sum(1 for o in out if o["best_n"] and (o["thm1_n"] is None or o["best_n"] < o["thm1_n"])))
