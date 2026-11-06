"""The cells first run on an arm64 laptop, rerun on x86-64 cluster nodes, compared class by class.

Cells (stored results, run on the laptop): l=2, k<=6 and l=4, k<=2 (outputs/staged_small), l=4, k=3
(outputs/staged_l4fix), and stages (i)-(ii) at l=3, k<=4 (outputs/staged_pc; stage (iii) there is the x86 chain).
Rerun on ccc4 with the same code and settings (LSF jobs 2246916-2246918, 2026-11-05, one core per task, x86_64,
z3 5.3.2), results copied to outputs/x86_rerun/. Commands (from the repository root on the cluster):
    python searches/staged.py --L 4 --k 3 --out outputs/x86_rerun/staged_l4fix --jobs 1
    python searches/staged.py --L 2 --k 1 2 3 4 5 6 --out outputs/x86_rerun/staged_small --jobs 1
    python searches/staged.py --L 4 --k 1 2 --out outputs/x86_rerun/staged_small_l4 --jobs 1
    python searches/staged_class.py --cell 3 k outputs/x86_rerun/staged_pc               (k = 1..4)
    python searches/staged_class.py --chunk 3 k I 2 outputs/x86_rerun/staged_pc          (I = 0, 1)
    python searches/staged_class.py --merge 3 k outputs/x86_rerun/staged_pc
Compares the final n per class (last pass present; ii_eg at l=3).
    python searches/x86_rerun_compare.py   -> outputs/checks/x86_rerun_compare.json
"""
import os, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs"); X = os.path.join(O, "x86_rerun")
PASSES = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg", "iii_r", "iii_f"]


def final(d, upto=None):
    ps = [p for p in PASSES if os.path.exists(os.path.join(d, f"{p}.json")) and (upto is None or PASSES.index(p) <= PASSES.index(upto))]
    return {o: int(b["n"]) for o, b in json.load(open(os.path.join(d, f"{ps[-1]}.json")))["best"].items()}


def main():
    cells = [(f"staged_small/l2_k{k}", f"staged_small/l2_k{k}", None) for k in range(1, 7)]
    cells += [(f"staged_small/l4_k{k}", f"staged_small_l4/l4_k{k}", None) for k in (1, 2)]
    cells += [("staged_l4fix/l4_k3", "staged_l4fix/l4_k3", None)]
    cells += [(f"staged_pc/l3_k{k}", f"staged_pc/l3_k{k}", "ii_eg") for k in range(1, 5)]
    rows, same, tot = [], 0, 0
    for a, b, upto in cells:
        A, B = final(os.path.join(O, a), upto), final(os.path.join(X, b), upto)
        s = sum(B.get(o) == n for o, n in A.items()); same += s; tot += len(A)
        rows.append(dict(cell=a, classes=len(A), same=s, differ=[o for o, n in A.items() if B.get(o) != n]))
    res = dict(same=same, total=tot, cells=rows, jobs=[2246916, 2246917, 2246918], host="ccc4, x86_64, z3 5.3.2")
    json.dump(res, open(os.path.join(O, "checks", "x86_rerun_compare.json"), "w"), indent=1)
    print(f"same n on x86-64 as on the laptop for {same} of {tot} classes")


if __name__ == "__main__":
    main()
