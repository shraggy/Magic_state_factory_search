"""How strong the synthillation benchmark is: the counts its Theorem 1 allows, case 2 of its Table I, and margins.

Campbell-Howard (arXiv:1606.01904v5, Thm 1 and Table I) give n = col(A) + 2 col(B) + Delta for any valid synthesis
matrix A (col(A) >= tau) and Lempel factor B (col(B) >= mu), with Delta fixed by the parities of col(A), col(B) and
whether the target has T (linear) terms. The class enters only through tau and mu, so every count the theorem allows
is reached with col(A) in {tau, tau+1} and col(B) in {mu, mu+1}; two more columns only add 2 or 4 inputs. This script
reads outputs/synthillation_strongest_l3.json (the variants searches/synthillation_strongest.py built and checked
with Eq. (4)) and outputs/staged_best.json (the cascade), and reports, over the 766 l=3 classes with k<=6:
  (a) mu: the Lempel rank (rank Q, +1 if diag Q = 0 and Q != 0, Q the quadratic form of the target) equals the stored
      mu on 766 of 766 (known-answer check, asserted);
  (b) validity per Table-I case over all built variants; case 2 (T terms, col(A) even, col(B) odd) fails on every
      instance. Structural reason: its second constant row as printed has a 1 in 2 col(B) + 5 columns (the two B
      blocks, the c columns marked in that row and three padding columns), an odd number, so the linear term on that
      check is a T rotation (computed from searches/synthillation_compare.C14[2]);
  (c) benchmark = min Theorem-1 count with the Delta realised by the printed matrices (case 2 excluded, case 6 at
      Delta = 2 since its matrix has two zero columns), on N of 766;
  (d) a generous bound: Table I's printed Delta for every case, case 2 counted as valid, and n >= tau+1; the benchmark
      against it and the cascade against it (by how many inputs below, on the classes with n < n_syn);
  (e) the synthillation paper's own instances with k<=6, by their class representatives here: CCZ 8, CCZ^2 14
      (CCZ045,CCZ123), Tof# 12 (CCZ014,CCZ023), CS 12, CS^2 18 (CS03,CS12), CS^3 26 (CS05,CS14,CS23);
  (f) exact coset minima on 6 qubits (searches/coset_certify.py) for T_1 CCZ_123 and T_3 CS_12, which the case-2
      formula would give at 17 and 19 on k+3 = 6 qubits;
  (g) the case-2 classes (T terms, tau even, mu odd): the benchmark, the formula tau + 2 mu + 9 and the cascade.
Deterministic. (f) takes about a minute.
    python searches/synthillation_bound.py   -> outputs/synthillation_bound_l3.json
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))

DL_PRINTED = {(0, 0): 8, (0, 1): 9, (1, 0): 9, (1, 1): 11}     # cases 1-4 (T terms), Table I as printed
DN_PRINTED = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 3}      # cases 5-8 (no T terms), as printed
DL_REAL = {(0, 0): 8, (1, 0): 9, (1, 1): 11}                    # case 2 removed (never valid)
DN_REAL = {(0, 0): 0, (0, 1): 2, (1, 0): 1, (1, 1): 3}          # case 6 has two zero columns: Delta = 2


def parse(spec):
    P = {}
    for g in spec.split(","):
        if g.startswith("CCZ"):
            P[tuple(int(c) for c in g[3:])] = 4
        elif g.startswith("CS"):
            P[tuple(int(c) for c in g[2:])] = 2
        else:
            P[(int(g[1:]),)] = 1
    return P


def rank2(rows):
    r = 0; rows = list(rows)
    for bit in range(8):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def lempel_mu(spec, k):
    P = parse(spec); lin = [P.get((i,), 0) % 2 for i in range(k)]
    Q = [lin[i] << i for i in range(k)]
    for T, v in P.items():
        if len(T) == 2 and (v // 2) % 2:
            i, j = T; Q[i] |= 1 << j; Q[j] |= 1 << i
    rk = rank2(Q)
    return 0 if rk == 0 else (rk if any(lin) else rk + 1)


def thm1(tau, mu, lin, DL, DN):
    D = DL if lin else DN
    v = min(a + 2 * b + D[(a % 2, b % 2)] for a in (tau, tau + 1) for b in (mu, mu + 1) if (a % 2, b % 2) in D)
    if not lin and mu == 0:
        v = min(v, tau if tau % 2 == 0 else tau + 1)             # cases 9 (tau even) and 11 (tau odd)
    return v


def main():
    import synthillation_compare as sc
    S = json.load(open(os.path.join(ROOT, "outputs", "synthillation_strongest_l3.json")))
    M = {(x["k"], x["spec"]): x for x in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]
         if x["l"] == 3}
    res = {}
    # (a)
    mu_ok = sum(lempel_mu(o["spec"], o["k"]) == o["mu"] for o in S)
    assert mu_ok == len(S) == 766, mu_ok
    res["mu_equals_lempel_rank"] = f"{mu_ok} of {len(S)}"
    # (b)
    val = collections.defaultdict(collections.Counter)
    for o in S:
        lin = any(g.startswith("T") for g in o["spec"].split(","))
        for v in o["variants"]:
            if v["valid"]:
                val[v["case"]]["valid"] += 1
            else:
                # case of an invalid variant from the parities (cases 1-8; 9/11 when mu = 0 without T terms)
                t, m = v["tau"], v["mu"]
                if not lin and m == 0:
                    c = "9" if t % 2 == 0 else "11"
                else:
                    c = str(({(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4} if lin else
                             {(0, 0): 5, (0, 1): 6, (1, 0): 7, (1, 1): 8})[(t % 2, m % 2)])
                val[c]["invalid"] += 1
    res["validity_by_case"] = {c: dict(val[c]) for c in sorted(val, key=int)}
    row2 = sc.C14[2][1]
    res["case2_row2_weight"] = "2 col(B) + %d (odd for every target)" % (sum(row2) - 2)
    assert (sum(row2) - 2) % 2 == 1 and sc.C14[2][1][1] == 1 and sc.C14[2][1][2] == 1
    res["case2_valid_instances"] = val["2"].get("valid", 0)
    # (c)
    c = collections.Counter(o["best_n"] - thm1(o["tau"], o["mu"], any(g.startswith("T") for g in o["spec"].split(",")),
                                               DL_REAL, DN_REAL) for o in S)
    res["benchmark_minus_realised_thm1"] = dict(sorted(c.items()))
    # (d)
    bg = collections.Counter(); cg = collections.Counter(); fewer = 0
    for o in S:
        lin = any(g.startswith("T") for g in o["spec"].split(","))
        lb = max(thm1(o["tau"], o["mu"], lin, DL_PRINTED, DN_PRINTED), o["tau"] + 1)
        bg[o["best_n"] - lb] += 1
        n = M[(o["k"], o["spec"])]["n"]
        if n < o["best_n"]:
            fewer += 1; cg[lb - n] += 1
    res["generous_bound"] = dict(benchmark_minus_bound=dict(sorted(bg.items())), classes_fewer_than_syn=fewer,
                                 cascade_below_bound_by=dict(sorted(cg.items())),
                                 cascade_below_bound_on=sum(v for g, v in cg.items() if g > 0))
    # (e)
    by = {o["spec"]: o for o in S}
    inst = {"CCZ012": 8, "CCZ045,CCZ123": 14, "CCZ014,CCZ023": 12, "CS01": 12, "CS03,CS12": 18, "CS05,CS14,CS23": 26}
    got = {}
    for spec, want in inst.items():
        o = by.get(spec)
        got[spec] = dict(paper=want, benchmark=o and o["best_n"])
    res["paper_instances"] = got
    res["paper_instances_reproduced"] = f'{sum(v["paper"] == v["benchmark"] for v in got.values())} of {len(got)}'
    # (f)
    import coset_certify as cc
    cert = {}
    for name, t, formula in (("T_1 CCZ_123", {frozenset([0]): 1, frozenset([0, 1, 2]): 4}, 17),
                             ("T_3 CS_12", {frozenset([2]): 1, frozenset([0, 1]): 2}, 19)):
        mins = {}
        for N in (4, 5, 6):
            r = cc.affine_space(t, 3, N)
            mins[N] = None if r is None else cc.min_weight(*r)
        cert[name] = dict(case2_formula=formula, exact_min_n_by_N=mins)
    res["coset_minima"] = cert
    # (g)
    c2 = [o for o in S if any(g.startswith("T") for g in o["spec"].split(",")) and o["tau"] % 2 == 0 and o["mu"] % 2 == 1]
    res["case2_classes"] = dict(
        classes=len(c2), best_variant=dict(collections.Counter(o["best_variant"] for o in c2)),
        unpadded_invalid_all_classes=sum(o["thm1_n"] is None for o in S),
        benchmark_minus_formula=dict(collections.Counter(o["best_n"] - (o["tau"] + 2 * o["mu"] + 9) for o in c2)),
        cascade_minus_formula=dict(sorted(collections.Counter(M[(o["k"], o["spec"])]["n"] - (o["tau"] + 2 * o["mu"] + 9)
                                                          for o in c2).items())),
        by_k=dict(sorted(collections.Counter(o["k"] for o in c2).items())))
    json.dump(res, open(os.path.join(ROOT, "outputs", "synthillation_bound_l3.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
