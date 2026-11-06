"""Per l=3 class: the largest N on which the all-coefficient free SAT query (z3) is unsatisfiable below the factory's n.
Unsatisfiable on N qubits means no factory with fewer inputs exists on N or fewer qubits (an idle check preserves a
factory). Reads the free_coef* records of the staged runs (outputs/staged_pc for k<=4, outputs/staged_A for k=5,6).
    python searches/free_sat_proofs.py   -> outputs/checks/free_sat_proofs.json
"""
import os, json, glob, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")


def main():
    best = {(r["k"], str(r["orbit"])): r for r in json.load(open(os.path.join(O, "staged_best.json")))["rows"] if r["l"] == 3}
    mx = collections.defaultdict(int)
    for k in range(1, 7):
        src = os.path.join(O, "staged_pc" if k <= 4 else "staged_A", f"l3_k{k}")
        for f in glob.glob(os.path.join(src, "free_coef*", "*.json")):
            for o, v in json.load(open(f))["per_class"].items():
                m = max([s["N"] for s in v.get("steps") or [] if s["last"] == "unsat"], default=0)
                mx[(k, o)] = max(mx[(k, o)], m)
    by_N = collections.Counter(mx.get(key, 0) for key in best)
    res = dict(classes=len(best), unsat_up_to_N={str(N): c for N, c in sorted(by_N.items(), reverse=True)},
               proved_on_own_N=sum(mx.get(key, 0) >= r["N"] for key, r in best.items()))
    json.dump(res, open(os.path.join(O, "checks", "free_sat_proofs.json"), "w"), indent=1)
    print(res)


if __name__ == "__main__":
    main()
