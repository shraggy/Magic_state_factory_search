"""Checks behind the first-submission (arXiv:2606.28518v1) material restored in the revision.

    nice -n 19 python3 searches/lostv1_checks.py      -> outputs/checks/lostv1_checks.json   (about 2 min, one core)

1. Side-by-side table (v1 Table V): every l=2,3,4 entry is rebuilt from its catalogue row
   (outputs/large_unified_factories.csv) with merge_large_unified.regenerate_circuit and re-verified with
   unified_search.is_factory, both in monomial form and on the full 2^N phase table (brute_force=True), and with
   unified_search.classify_factory (n, k, d, Clifford-reduced output).
2. The C^{l-1}Z family at (s_total, s_O) = (1, 1), k = l: n = 2^N - 2^l on every catalogue row; the members with
   N <= 11 (up to [[2032,4,2]] at l=4) are rebuilt and re-verified as in 1.
3. The six-qubit parent of v1 Fig. 6(c): its cuts are verified with is_factory and classified.
4. Independent scripts in searches/lostv1/ (each its own code, no repository import), run and their output kept.
"""
import csv, json, os, subprocess, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
from merge_large_unified import regenerate_circuit

R = list(csv.DictReader(open(os.path.join(ROOT, "outputs", "large_unified_factories.csv"))))
out = {}


def verify(row):
    l, k, N = int(row["l"]), int(row["k"]), int(row["N"])
    g = regenerate_circuit(l, k, N, row["params"], row["search_stage"], row["raw"])
    f = dict(stage="x", k=k, r=N - k, gates=g, params="")
    info = us.classify_factory(f, l)
    return dict(l=l, n=int(row["n"]), k=k, N=N, stage=row["search_stage"], params=row["params"],
                is_factory=us.is_factory(f, l), brute_force=us.is_factory(f, l, brute_force=True) if N <= 12 else None,
                n_rebuilt=info["N"] if info else None, d=info["d"] if info else None,
                reduced=info["reduced"] if info else None, reduced_csv=row["reduced"])


def pick(l, n, k, stage="two-group", params=None):
    m = [r for r in R if (int(r["l"]), int(r["n"]), int(r["k"])) == (l, n, k) and r["search_stage"] == stage
         and (params is None or r["params"] == params)]
    return min(m, key=lambda r: int(r["N"])) if m else None


# 1. v1 Table V rows, by construction, at l = 2, 3, 4
V1 = {"RM all weights": [(2, 6, 1), (3, 14, 1), (4, 30, 1)],
      "RM odd weights": [(2, 7, 1), (3, 15, 1), (4, 31, 1)],
      "two-group deg 1, k=2": [(2, 6, 2), (3, 14, 2), (4, 30, 2)],
      "two-group deg 1, k=4": [(2, 8, 4), (3, 20, 4), (4, 44, 4)],
      "two-group deg 1, k=6": [(2, 10, 6), (3, 26, 6), (4, 58, 6)],
      "two-group deg 2 (CZ, CS, CT)": [(2, 4, 2), (3, 12, 2), (4, 28, 2)],
      "two-group deg 3 (CCZ, CCS)": [(3, 8, 3), (4, 24, 3)],
      "two-group deg 4 (CCCZ)": [(4, 16, 4)]}
PAR = {"two-group deg 1, k=2": "s_total=2,s_O=1", "two-group deg 1, k=4": "s_total=2,s_O=3",
       "two-group deg 1, k=6": "s_total=2,s_O=5"}
tab = {}
for name, mem in V1.items():
    tab[name] = []
    for l, n, k in mem:
        r = pick(l, n, k, params=PAR.get(name))
        tab[name].append(verify(r) if r else dict(l=l, n=n, k=k, missing=True))
out["v1_table_side_by_side"] = tab
ok = [e for v in tab.values() for e in v if e.get("is_factory") and e.get("brute_force") and e.get("n_rebuilt") == e["n"]]
out["v1_table_verified"] = f"{len(ok)} of {sum(len(v) for v in tab.values())}"

# One choice of steps, (s_total, s_O) = (2, k-1), even k: smallest-N catalogue row per (l, k); expected
# n = (2^{l-1}-1)k + 2^l on N = k + l qubits (H-code parameters at l=2, Bravyi--Haah at l=3, sqrt(T) analogue at l=4)
bh = []
for l in (2, 3, 4):
    for k in range(2, 13, 2):
        m = [r for r in R if r["l"] == str(l) and int(r["k"]) == k and r["search_stage"] == "two-group"
             and r["params"] == f"s_total=2,s_O={k - 1}"]
        if m:
            r = min(m, key=lambda r: int(r["N"]))
            bh.append(dict(l=l, k=k, N=int(r["N"]), n=int(r["n"]), formula=(2 ** (l - 1) - 1) * k + 2 ** l,
                           N_is_k_plus_l=int(r["N"]) == k + l, reduced=r["reduced"]))
out["two_group_s2_family"] = dict(rows=bh, formula_holds=f"{sum(b['n'] == b['formula'] and b['N_is_k_plus_l'] for b in bh)} of {len(bh)}")

# 2. C^{l-1}Z family, (1,1), k = l
fam = {}
for l, red in ((2, r"\mathrm{CZ}_{12}"), (3, r"\mathrm{CCZ}_{123}"), (4, r"\mathrm{CCCZ}_{1234}")):
    rows = sorted([r for r in R if r["l"] == str(l) and int(r["k"]) == l and r["search_stage"] == "two-group"
                   and r["params"] == "s_total=1,s_O=1" and r["reduced"] == red], key=lambda r: int(r["N"]))
    fam[l] = dict(members=[(int(r["N"]), int(r["n"])) for r in rows],
                  n_eq_2N_minus_2l=f"{sum(int(r['n']) == 2 ** int(r['N']) - 2 ** l for r in rows)} of {len(rows)}",
                  reverified=[verify(r) for r in rows if int(r["N"]) <= 11])
out["CkZ_family"] = fam

# 3. v1 Fig. 6(c) six-qubit parent (Quirk link 39 of the v1 PDF)
six = [(0, 1, 3, 4, 5), (0, 2, 3, 4, 5), (0, 1, 5), (0, 2, 5), (0, 1, 4), (0, 2, 4), (0, 1, 3), (0, 2, 3),
       (1, 2, 3, 4, 5), (1, 2, 5), (1, 2, 4), (1, 2, 3), (3, 4, 5)]
sp = dict(full_is_BI=us.is_factory(dict(k=0, r=6, gates=[(frozenset(S), 1) for S in six]), 3, brute_force=True))
for checks in ([0], [0, 1, 2]):
    outs = [q for q in range(6) if q not in checks]; m = {q: i for i, q in enumerate(outs + checks)}
    keep = [(frozenset(m[q] for q in S), 1) for S in six if set(S) & set(checks)]
    f = dict(k=len(outs), r=len(checks), gates=keep, stage="x", params="")
    info = us.classify_factory(f, 3)
    sp[f"cut_checks_{checks}"] = dict(is_factory=us.is_factory(f, 3, brute_force=True), n=info["N"], k=len(outs),
                                      d=info["d"], essential_dim=info["essential_dim"], reduced=info["reduced"])
out["v1_six_qubit_parent"] = sp

# 4. independent scripts
D = os.path.join(HERE, "lostv1")
for s in ("padded_842.py", "cube_output_phase.py", "small_n_l3.py", "even_weights.py", "relabel_invariance.py",
          "sym_distance.py", "six_qubit_parent.py", "clifford_orbit_ccz.py", "catalyst_split.py"):
    p = subprocess.run([sys.executable, s], cwd=D, capture_output=True, text=True)
    out[s] = p.stdout.strip().splitlines()[-25:]
json.dump(out, open(os.path.join(ROOT, "outputs", "checks", "lostv1_checks.json"), "w"), indent=1, default=str)
print(json.dumps({k: (v if k in ("v1_table_verified",) else "...") for k, v in out.items()}, indent=1))
print(out["v1_table_verified"], out["two_group_s2_family"]["formula_holds"],
      {l: (v["n_eq_2N_minus_2l"], len(v["reverified"]), all(e["is_factory"] and e["brute_force"] for e in v["reverified"])) for l, v in fam.items()})
