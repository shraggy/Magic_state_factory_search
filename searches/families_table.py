"""Table SI: the closed-form families that the two cheap stages of the cascade produce, each member looked up in the
larger-range factory catalogue (outputs/large_unified_factories.csv: two-group N<=18, gluing k<=9 at l=2,3; N<=16 at
l=4). A member is listed only if a verified factory with that [[n,k,2]], output class and stage is in the file.
    python searches/families_table.py     -> outputs/checks/families_table.json
"""
import os, re, csv, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
R = list(csv.DictReader(open(os.path.join(O, "large_unified_factories.csv"))))
def find(l, n, k, stage, red, d="2"):
    m = [r for r in R if r["l"] == str(l) and int(r["n"]) == n and int(r["k"]) == k and r["search_stage"] == stage
         and r["d"] == d and re.fullmatch(red, r["reduced"])]
    return min(m, key=lambda r: int(r["N"])) if m else None
T = r"\\mathrm\{T\}_\{\d+\}"; CS = r"\\mathrm\{CS\}_\{\d+\}"; CCZ = r"\\mathrm\{CCZ\}_\{\d+\}"
SK = r"\\mathrm\{(S|T|\\sqrt\{T\})\}_\{\d+\}"
fam = [
 ("two-group", r"Reed--Muller, all weights", [(l, 2 ** (l + 1) - 2, 1, r"\\mathrm\{(S|T|\\sqrt\{T\})\}_\{1\}", "2") for l in (2, 3, 4)]),
 ("two-group", r"Reed--Muller, odd weights", [(l, 2 ** (l + 1) - 1, 1, r".*", "3") for l in (2, 3, 4)]),
 ("two-group", r"hypercube $[[2^l,l,2]]$", [(l, 2 ** l, l, r"\\mathrm\{(CZ|CCZ|CCCZ)\}_\{\d+\}", "2") for l in (2, 3, 4)]),
 # malleable cuts of the all-weights parent, Eq. (chain), with two and three outputs (v1 Table V, degree 2 and 3)
 ("two-group", r"chain, two outputs $[[2^{l+1}-4,2,2]]$", [(l, 2 ** (l + 1) - 4, 2, r"\\mathrm\{(CZ|CS|CT)\}_\{12\}", "2") for l in (2, 3, 4)]),
 ("two-group", r"chain, three outputs $[[2^{l+1}-8,3,2]]$", [(l, 2 ** (l + 1) - 8, 3, r"\\mathrm\{(CCZ|CCS)\}_\{123\}", "2") for l in (3, 4)]),
 # one choice of steps, (s_total, s_O) = (2, k-1): H-code parameters (l=2), Bravyi--Haah (l=3), sqrt(T) analogue (l=4)
 ("two-group", r"$[[(2^{l-1}-1)k+2^l,k,2]]$, $k$ even", [(l, (2 ** (l - 1) - 1) * k + 2 ** l, k, f"({SK})+", "2") for l in (2, 3, 4) for k in (2, 4, 6)]),
 # hypercube cut with N-l checks: every rotation that touches a check, (s_total, s_O) = (1, 1), k = l
 ("two-group", r"$\\CkZ{l-1}$, $[[2^N-2^l,l,2]]_N$", [(l, 2 ** N - 2 ** l, l, r"\\mathrm\{(CZ|CCZ|CCCZ)\}_\{\d+\}", "2") for l in (2, 3, 4) for N in range(l + 1, l + 8)]),
 ("gluing", r"$m$ CCZ, $[[6m+2,3m,2]]$", [(3, 6 * m + 2, 3 * m, f"({CCZ})+", "2") for m in (2, 3)]),
 ("gluing", r"$m$ CS, $[[6m+6,2m,2]]$", [(3, 6 * m + 6, 2 * m, f"({CS})+", "2") for m in (2, 3)]),
]
out = []
for stage, name, members in fam:
    got = []
    for l, n, k, red, d in members:
        r = find(l, n, k, stage, red, d)
        if r: got.append(dict(l=l, n=n, k=k, N=int(r["N"]), d=int(d), reduced=r["reduced"], params=r["params"]))
    out.append(dict(stage=stage, family=name, members=got, asked=len(members)))
    print(f"{stage:9s} {name:32s} {len(got)}/{len(members)}:", [(g['l'], g['n'], g['k'], g['N']) for g in got])
os.makedirs(os.path.join(O, "checks"), exist_ok=True)
json.dump(out, open(os.path.join(O, "checks", "families_table.json"), "w"), indent=1)
