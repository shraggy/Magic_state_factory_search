"""Stages (iii)-(iv) of the unified search: exact target list, matched to the catalogue, as SAT tasks.

For each (l, k) the orbit table of enumerate_targets.py gives one target per Clifford class. Every catalogue factory
(stages i-ii) is mapped to its exact class, so each target carries the best input count the SAT-free stages reached
(or none). One 'cascade' task per genuine target: SAT scans circuit sizes upward from k+1 qubits and, at each size,
asks only for a factory with fewer inputs than the best known so far (from the catalogue or from smaller sizes).

    python searches/stage3_tasks.py --cells 2:1-6 3:1-5 4:1-3 --qmax 9 --out tasks_cascade.txt
writes the task file (one line per target, run_task.py format) and outputs/targets/catalogue_match.json
(exact class of every catalogue factory, and the exact class count per (l, k)).
"""
import os, sys, csv, json, time, argparse
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search  # noqa: E402,F401  (puts classification/ on sys.path)
import merge_large_unified as M  # noqa: E402
from classify import phase_polynomial  # noqa: E402
from enumerate_targets import enumerate_targets, index_of_poly  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CATALOGUES = [os.path.join(ROOT, "outputs", f) for f in ("all_factories.csv", "large_unified_factories.csv")]


def catalogue_rows(L, k):
    """Distance-2 catalogue rows at (l, k), paper convention (n = inputs, N = qubits), with their circuits."""
    out = []
    for path in CATALOGUES:
        for r in csv.DictReader(open(path)):
            if int(r["l"]) != L or int(r["k"]) != k or r["d"] != "2" or r["search_stage"] in ("padding", "SAT"):
                continue
            out.append(r)
    return out


def output_poly(r, L, k):
    g = M.regenerate_circuit(L, k, int(r["N"]), r["params"], r["search_stage"], r["raw"])
    return phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), L)


def parse_cells(cells):
    out = []
    for c in cells:
        L, ks = c.split(":")
        a, b = map(int, ks.split("-"))
        out += [(int(L), k) for k in range(a, b + 1)]
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", nargs="+", default=["2:1-6", "3:1-5", "4:1-3"])
    ap.add_argument("--qmax", type=int, default=9, help="largest circuit size (qubits) SAT is asked about")
    ap.add_argument("--time", type=int, default=86400, help="wall-clock limit per task, seconds")
    ap.add_argument("--out", default=os.path.join(ROOT, "tasks_cascade.txt"))
    a = ap.parse_args()
    lines, report = [], {}
    for L, k in parse_cells(a.cells):
        t0 = time.time()
        res, (C, off, lab) = enumerate_targets(L, k)
        t_enum = time.time() - t0
        best = {}                                         # orbit -> (inputs, qubits, stage)
        seen_labels = defaultdict(set)                    # orbit -> catalogue 'reduced' labels
        for r in catalogue_rows(L, k):
            orb = int(lab[index_of_poly(output_poly(r, L, k), k, L, C, off)])
            seen_labels[orb].add(r["reduced"])
            cand = (int(r["n"]), int(r["N"]), r["search_stage"])
            if orb not in best or cand < best[orb]:
                best[orb] = cand
        gen = [c for c in res["classes"] if c["genuine"]]
        reached = [c for c in gen if c["orbit"] in best]
        for c in gen:
            b = best.get(c["orbit"])
            label = f"l{L}_k{k}_o{c['orbit']}"
            lines.append(f"cascade {label} {L} {c['spec']} {k} {a.qmax} {b[0] if b else 1000000} {a.time}")
        report[f"l{L}_k{k}"] = dict(
            l=L, k=k, genuine_classes=len(gen), reached_by_stages=len(reached),
            catalogue_labels=sum(len(seen_labels[c["orbit"]]) for c in reached),
            merged=[sorted(v) for o, v in seen_labels.items() if len(v) > 1],
            targets=[dict(orbit=c["orbit"], spec=c["spec"], best_stage=best.get(c["orbit"])) for c in gen],
            enumerate_seconds=round(t_enum, 2), total_seconds=round(time.time() - t0, 2))
        print(f"l={L} k={k}: {len(gen)} genuine classes, {len(reached)} reached by stages i-ii "
              f"({report[f'l{L}_k{k}']['catalogue_labels']} catalogue labels), {time.time() - t0:.1f}s", flush=True)
    open(a.out, "w").write("\n".join(lines) + "\n")
    os.makedirs(os.path.join(ROOT, "outputs", "targets"), exist_ok=True)
    json.dump(report, open(os.path.join(ROOT, "outputs", "targets", "catalogue_match.json"), "w"), indent=1)
    print(f"{len(lines)} cascade tasks -> {a.out}")
