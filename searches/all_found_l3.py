"""Every factory recorded at l=3, k<=6 by stages (i)-(ii) of sweep A and stage (iii) of the final pipeline (not only the smallest per class), and the sizes (n, k) of the
Nezami--Haah classification (outputs/checks/nezami_haah_fig1.json) at which one is a T^{otimes k} factory; read by
searches/plot_literature.py (Fig. 4).

Sources: stages (i)-(ii) of sweep A (cell/*.json: best and rows, with exact class; cls/*.json, its non-(iii) results)
and stage (iii) of the final pipeline (chain records of the three deterministic runs outputs/solver_v2/staged_A_v4,
_v5, _v6, whose per-class best the paper reports; z3 synthesis excluded; --runs v4 gives the earlier v4-only set). The passes ask only for fewer inputs than the best so far, so larger factories are
recorded only as far as the passes met them on the way down.
    python searches/all_found_l3.py   -> outputs/checks/all_found_l3.json
"""
import os, json, glob, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
CELLS = [("staged_pc", k) for k in range(1, 5)] + [("staged_A", 5), ("staged_A", 6)]   # stages (i)-(ii) and their recorded finds
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--runs", default="v4,v5,v6", help="comma list of stage-(iii) runs, e.g. v4,v5,v6")
_p.add_argument("--out-dir", default=None, help="write all_found_l3.json here instead of outputs/checks")
_A = _p.parse_args()
CHAINS = [os.path.join(O, "solver_v2", f"staged_A_{v}") for v in _A.runs.split(",")]   # stage (iii) of the final pipeline


def collect():
    found, spec = collections.defaultdict(set), {}            # (n, k) -> classes
    for src, k in CELLS:
        d = os.path.join(O, src, f"l3_k{k}")
        spec.update({(k, o): s for o, s in json.load(open(os.path.join(d, "iii_r.json")))["specs"].items()})
        for f in glob.glob(os.path.join(d, "cell", "*.json")):
            x = json.load(open(f))
            for o, b in x["best"].items(): found[(b["n"], k)].add(o)
            for r in x.get("rows") or []: found[(r["n"], k)].add(str(r["orbit"]))
        for f in glob.glob(os.path.join(d, "cls", "*.json")):     # stage (iii) results of sweep A are not used
            x = json.load(open(f)); res = x["results"] if isinstance(x["results"], list) else eval(x["results"])
            for r in res:
                if not r["stage"].startswith("iii"): found[(r["n"], k)].add(str(r["orbit"]))
        for p in ("i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg"):        # one-process cells (k <= 4)
            f = os.path.join(d, f"{p}.json")
            if os.path.exists(f) and not glob.glob(os.path.join(d, "cls", "*.json")):
                for o, b in json.load(open(f))["best"].items(): found[(b["n"], k)].add(o)
        for f in [f for CHAIN in CHAINS for f in glob.glob(os.path.join(CHAIN, f"l3_k{k}", "chain", "*.json"))]:
            for o, v in json.load(open(f))["per_class"].items():
                for st in v["stages"]:
                    if st.get("n") is not None and "skipped" not in st and st["stage"] != "syn":
                        found[(st["n"], k)].add(o)
    return found, spec


def main():
    found, spec = collect()
    T = {key for key, s in spec.items() if all(t.startswith("T") for t in s.split(","))}
    # T^{otimes k} by class, not by the representative's spec: at k=6 the class list stores T^{otimes 6} under a
    # mixed-gate representative (orbit 2352:2177). Label T0,...,T_{k-1} with the cell's exact labeller.
    import sys as _s; _s.path[:0] = [os.path.join(ROOT, "searches"), os.path.join(ROOT, "classification")]
    import cascade_class as CC
    from cascade_class import index_of_poly
    from coset_certify import parse_target
    for k in range(1, 7):
        classes, ctx = CC.classes_of(3, k, os.path.join(O, "staged_A", "_cache")); C, off, lab, fiber = ctx
        P = {frozenset(t): c for t, c in parse_target(",".join(f"T{q}" for q in range(k)), 3).items()}
        if lab is not None:
            o = str(int(lab[index_of_poly(P, k, 3, C, off)]))
        else:
            lb = fiber.label(fiber.digits_of_poly(P)); o = lb if isinstance(lb, str) else f"{lb[0]}:{lb[1]}"
        T.add((k, o))
    tsize = {(n, k) for (n, k), cl in found.items() if any((k, o) in T for o in cl)}
    NH = [tuple(x) for x in json.load(open(os.path.join(O, "checks", "nezami_haah_fig1.json")))["points"]]
    hit = [p for p in NH if p in tsize]
    res = dict(sizes=len(found), factories=sum(len(v) for v in found.values()), t_sizes=sorted(tsize),
               nezami_haah=len(NH), nezami_haah_hit=len(hit), hit=sorted(hit), missed=sorted(set(NH) - set(hit)))
    json.dump(res, open(os.path.join(_A.out_dir or os.path.join(O, "checks"), "all_found_l3.json"), "w"), indent=0)
    print({k: res[k] for k in ("sizes", "factories", "nezami_haah", "nezami_haah_hit")})


if __name__ == "__main__":
    main()
