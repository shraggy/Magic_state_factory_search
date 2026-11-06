"""Output-error coefficient c of the cascade's factory against that of the synthillation construction, all l=3 classes.

At distance 2 the output error is c p^2 + O(p^3), where c counts the unordered pairs of faulty input rotations whose
check parts are equal (so the pair is not detected) and whose output parts differ (so the pair changes the output).
c is the total over all k outputs of the factory. For each of the 766 classes in outputs/class_table_l3.json:
  - n and c of the factory in outputs/staged_best.json (recounted from its gate list and checked against the stored
    n and c);
  - every padded synthillation variant that searches/synthillation_strongest.py builds (A in {A*, A*+0},
    B in {B*, B*+0, B(mu*+1)}, Campbell-Howard Thm 1 / Table I), with n, N and c, valid ones only;
  - n_syn = smallest valid n (checked against best_n in outputs/synthillation_strongest_l3.json) and c_syn = the
    smallest c among the valid variants with n = n_syn; c_syn_any = smallest c over all valid variants at any n.
The synthillation c is that of the construction built here (one optimal A and one Lempel factor B per class, as in
synthillation_strongest.py), not a minimum over all syntheses.
Deterministic (z3 models of synthillation_compare.opt_B with fixed constraints; no randomness).
    nice -n 19 python searches/syn_c_compare.py --workers 2   -> outputs/checks/syn_c_compare.json  (about 1 h)
Known-answer checks (asserted): c and n recount equal staged_best.json on 766 of 766; n_syn equals best_n of
synthillation_strongest_l3.json on 766 of 766.
"""
import os, sys, json, argparse, collections
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def c_of(faults):
    """faults: list of (check part, output part) of the input rotations; returns c."""
    g = collections.defaultdict(list)
    for ch, o in faults:
        g[ch].append(o)
    bad = 0
    for L in g.values():
        m = len(L); C = collections.Counter(L)
        bad += m * (m - 1) // 2 - sum(v * (v - 1) // 2 for v in C.values())
    return bad


def factory_nc(row):
    k = row["k"]; G = json.loads(row["gates"]) if isinstance(row["gates"], str) else row["gates"]
    f = [(sum(1 << (q - k) for q in lam if q >= k), sum(1 << q for q in lam if q < k)) for lam, c in G if c % 2]
    return len(f), c_of(f)


def cols_of(k, target, A, B, lin, c):
    import synthillation_compare as sc
    tau, mu = len(A), len(B)
    if lin:
        tm = sc.C14[{(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(tau % 2, mu % 2)]]
        blocks = [A, B, B] + [[c]] * 4
    elif mu == 0:
        cols = [a | 1 << k for a in A] + ([1 << k] if tau % 2 else [])
        ok = sc.coeffs_ok(cols, k + 1, k, target, 3) and all(y >> k for y in cols)
        return (cols, k + 1) if ok else None
    else:
        tm = sc.T5[(tau % 2, mu % 2)]; blocks = [A, B, B]
    blocks += [[0]] * (len(tm[0]) - len(blocks)); R = len(tm)
    cols = [v | sum(1 << (k + r) for r in range(R) if tm[r][bi]) for bi, blk in enumerate(blocks) for v in blk]
    ok = sc.coeffs_ok(cols, k + R, k, target, 3) and all(y >> k for y in cols)
    return (cols, k + R) if ok else None


def one(r):
    import synthillation_compare as sc, class_table_l3 as CT
    k = r["k"]; P = CT.parse(r["spec"]); A = sc.opt_A(k, P); mu = 0
    while sc.opt_B(k, P, mu) is None:
        mu += 1
    B0 = sc.opt_B(k, P, mu); B1 = sc.opt_B(k, P, mu + 1)
    lin = any(len(T) == 1 for T in P); c = sum(1 << i for i in range(k) if frozenset([i]) in P)
    res = []
    for an, Av in (("A", A), ("A+0", A + [0])):
        for bn, Bv in (("B", B0), ("B+0", B0 + [0]), ("B(mu+1)", B1)):
            if Bv is None:
                continue
            x = cols_of(k, P, Av, Bv, lin, c)
            if x:
                cols, N = x; mk = (1 << k) - 1
                res.append(dict(A=an, B=bn, n=len(cols), N=N, c=c_of([(y >> k, y & mk) for y in cols])))
    return dict(k=k, spec=r["spec"], variants=res)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=1); a = ap.parse_args()
    R = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))
    SB = {(x["k"], x["spec"]): x for x in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]
          if x["l"] == 3}
    SS = {(x["k"], x["spec"]): x for x in json.load(open(os.path.join(ROOT, "outputs", "synthillation_strongest_l3.json")))}
    with Pool(a.workers) as p:
        out = list(p.imap(one, R, chunksize=2))
    rows = []; recount_ok = nsyn_ok = 0; S = collections.Counter(); ratios = []; cfr = []; csfr = []
    for o in out:
        key = (o["k"], o["spec"]); fr = SB[key]; n, c = factory_nc(fr)
        recount_ok += (n == fr["n"] and c == fr["c"])
        V = o["variants"]; nsyn = min(v["n"] for v in V)
        nsyn_ok += (nsyn == SS[key]["best_n"])
        csyn = min(v["c"] for v in V if v["n"] == nsyn); cany = min(v["c"] for v in V)
        if n < nsyn:
            kind = "fewer_n_lower_c" if c < csyn else ("fewer_n_equal_c" if c == csyn else "fewer_n_higher_c")
            ratios.append(c / csyn); cfr.append(c); csfr.append(csyn)
        elif n == nsyn:
            kind = "same_n_lower_c" if c < csyn else ("same_n_equal_c" if c == csyn else "same_n_higher_c")
        else:
            kind = "more_n"
        S[kind] += 1
        S["any_variant_lower_c"] += (cany < c)
        S["c_equals_n_choose_2"] += (c == n * (n - 1) // 2)
        rows.append(dict(k=o["k"], orbit=fr["orbit"], spec=o["spec"], n=n, N=fr["N"], c=c, tau=fr["tau"],
                         n_syn=nsyn, c_syn=csyn, c_syn_any=cany, kind=kind, variants=V))
    assert recount_ok == len(R), recount_ok
    assert nsyn_ok == len(R), nsyn_ok
    ratios.sort()
    summary = dict(classes=len(R), recount_ok=recount_ok, n_syn_matches_strongest=nsyn_ok, **S,
                   fewer_ratio_c_over_csyn=dict(min=round(ratios[0], 3), median=round(ratios[len(ratios) // 2], 3),
                                                max=round(ratios[-1], 3)) if ratios else None,
                   fewer_c_range=[min(cfr), max(cfr)] if cfr else None,
                   fewer_csyn_range=[min(csfr), max(csfr)] if csfr else None,
                   on_bound_c_equals_n_choose_2=sum(1 for r in rows if r["n"] == r["tau"] + 1
                                                    and r["c"] == r["n"] * (r["n"] - 1) // 2),
                   on_bound=sum(1 for r in rows if r["n"] == r["tau"] + 1))
    os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
    json.dump(dict(source="outputs/staged_best.json, outputs/synthillation_strongest_l3.json", summary=summary,
                   rows=rows), open(os.path.join(ROOT, "outputs", "checks", "syn_c_compare.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
