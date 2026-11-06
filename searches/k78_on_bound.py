"""Which CCZ layers on the k=7, 8 target list have a factory on the bound n = tau(U)+1.

The pass files (outputs/k78_pipeline/full/l3_k{7,8}/ii_eg.json) mark a layer on the bound when its T-count is decided:
exactly at k=7 (tcount7.exists_below2), and at k=8 when the isotropy floor equals the upper bound of the best
compilation. A factory also bounds the T-count itself: by Corollary 3, n >= tau(U) + r, with r the GF(2) rank of the
check parts of its input rotations, so tau(U) <= n - r. Where that meets the isotropy floor, tau is decided and the
factory is on the bound. This adds the layers whose compilation upper bound was looser than the factory's own.
    python searches/k78_on_bound.py   -> outputs/k78_pipeline/full/on_bound.json
"""
import os, json, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = os.path.join(ROOT, "outputs", "k78_pipeline")


def rank_of_checks(gates, k):
    basis = []
    for S, c in gates:
        if c % 2 == 0:                                   # even coefficient: a Clifford rotation, consumes no input
            continue
        v = sum(1 << (q - k) for q in S if q >= k)
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def main():
    res = {}
    T = {}
    for f in glob.glob(os.path.join(P, "all_isd", "*.json")):
        for t in json.load(open(f))["results"]:
            T[t["id"]] = t
    for k in (7, 8):
        rows = json.load(open(os.path.join(P, "full", f"l3_k{k}", "ii_eg.json")))["rows"]
        flagged = factory_bound = 0
        for r in rows:
            t = T[r["id"]]; b = t["best"]
            lo = t["tcount_lo"]; hi = min(t["tcount_hi"], b["n"] - rank_of_checks(b["gates"], k))
            assert b["n"] == r["n"] and lo <= hi
            on = lo == hi and b["n"] == lo + 1
            flagged += bool(r["on_bound"]); factory_bound += on
            assert not (r["on_bound"] and not on)          # every layer flagged before stays on the bound
        res[str(k)] = dict(targets=len(rows), on_bound_in_pass_file=flagged, on_bound_with_factory_bound=factory_bound)
    json.dump(res, open(os.path.join(P, "full", "on_bound.json"), "w"), indent=1)
    print(res)


if __name__ == "__main__":
    main()
