"""The free diagonal gate D that a factory leaves on its checks and outputs, for every factory in outputs/staged_best.json.

A factory's kept rotations satisfy Eq. (4) (factory condition): every term of the phase polynomial that meets a check is
a free diagonal gate (Pauli Z at l=2; S/CZ-type Clifford at l>=3). Their product is D, so the circuit gives
D (U^dag|+>^k)|+>^(N-k), and D must be undone (or folded into the readout, measuring D X_q D^dag) before the X
measurements of the checks. This script expands each factory's rotations into monomials, coefficient
(-2)^(|Q|-1) sum_{Lambda contains Q} c_Lambda mod 2^l, keeps the terms that meet a check, and records per level:
    rows, not a factory (any non-free term: must be 0), D trivial, D nontrivial, D only Pauli Z on checks
    (a classical relabelling of the X outcomes), D with a CZ between a check and an output, D with a CZ between two
    checks, D with an S-type term on a check; and at l=3 the number of CZ terms of D (median, max, total).
Deterministic; reads outputs/staged_best.json only.
    python searches/clifford_correction.py   -> outputs/checks/clifford_correction.json
Known-answer check (asserted): [[8,3,2]] (CCZ012, l=3) has D trivial; every row passes the free-gate test.
"""
import os, json, itertools, collections, statistics
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def d_terms(row):
    L = row["l"]; M = 1 << L; k = row["k"]
    G = json.loads(row["gates"]) if isinstance(row["gates"], str) else row["gates"]
    coef = collections.defaultdict(int)
    for lam, c in G:
        for t in range(1, min(len(lam), L) + 1):
            for Q in itertools.combinations(sorted(lam), t):
                coef[Q] += c * (-2) ** (t - 1)
    return {Q: v % M for Q, v in coef.items() if v % M and any(q >= k for q in Q)}


def free(L, Q, v):
    t = len(Q)
    if L == 2:
        return t == 1 and v % 2 == 0
    return (t == 1 and v % (1 << (L - 2)) == 0) or (t == 2 and v % (1 << (L - 1)) == 0)


def main():
    rows = json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]
    out = {}
    for L in (2, 3, 4):
        R = [r for r in rows if r["l"] == L]; M = 1 << L; st = collections.Counter(); cz = []
        for r in R:
            k = r["k"]; D = d_terms(r)
            if any(not free(L, Q, v) for Q, v in D.items()):
                st["not_factory"] += 1; continue
            ncz = sum(1 for Q in D if len(Q) == 2); cz.append(ncz)
            if not D:
                st["D_trivial"] += 1; continue
            st["D_nontrivial"] += 1
            if all(len(Q) == 1 and v % (M // 2) == 0 for Q, v in D.items()):
                st["D_only_PauliZ_on_checks"] += 1
            if any(len(Q) == 2 and min(Q) < k for Q in D):
                st["D_CZ_check_output"] += 1
            if any(len(Q) == 2 and min(Q) >= k for Q in D):
                st["D_CZ_check_check"] += 1
            if any(len(Q) == 1 and v % (M // 2) for Q, v in D.items()):
                st["D_S_on_check"] += 1
        st["D_needs_Clifford_correction"] = st["D_nontrivial"] - st["D_only_PauliZ_on_checks"]
        out[f"l{L}"] = dict(rows=len(R), not_factory=st["not_factory"], **{k: v for k, v in st.items() if k != "not_factory"},
                           CZ_terms=dict(median=statistics.median(cz), max=max(cz), total=sum(cz)) if cz else None)
        assert st["not_factory"] == 0, (L, st)
    r832 = next(r for r in rows if r["l"] == 3 and r["spec"] == "CCZ012")
    assert not d_terms(r832), "[[8,3,2]] should have D trivial"
    os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
    json.dump(dict(source="outputs/staged_best.json", per_level=out),
              open(os.path.join(ROOT, "outputs", "checks", "clifford_correction.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
