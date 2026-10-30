"""Factory-level benchmark in the metrics of Gong et al. (arXiv:2608.09727, App. D), for the Fig. 1 factories and the
prior construction for each target. Single-level distance-2 factories, injected T (l=3) or sqrt(T) (l=4) states.

Per factory (n inputs, N qubits, k outputs), at input error p (every input carries an independent Z fault):
  acceptance            P_acc = 1 - n p                       (first order; every single fault is detected at d=2)
  output error          c p^2 (any output wrong, accepted)   c = undetected weight-2 fault pairs, as in Fig. 1
  per-state error       e p^2, e = average over outputs of the pairs that corrupt that output (Gong et al. D.2)
  time rate             (outputs per run) / (inputs per run) * P_acc   (unit of time = one injection)
  volume rate           time rate / N                         (spatial footprint = circuit qubits, uncompressed)
The target is one k-qubit state, so a run yields one target state; the time rate is 1 / n per target state.
The prior construction distils the parts with independent factories, run in parallel: its inputs add, its
footprint is the sum of the parts' footprints, and its error is the sum of the parts' errors.

    python searches/benchmark_factories.py          # prints a table and writes outputs/benchmark_factories.json
"""
import csv, json, os, sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
csv.field_size_limit(sys.maxsize)
import unified_search as us  # noqa: E402
from error_coefficients import gates_of  # noqa: E402

ROOT = os.path.dirname(HERE)
P_VALUES = (1e-3, 1e-6)


def fault_stats(gates, k, L):
    """n, N, c (pairs corrupting any output) and per-output pair counts, for a d=2 factory."""
    out = set(range(k)); faults = []
    for S, c in gates:
        for _ in range(us.inputs(c, L)):
            faults.append((frozenset(S - out), frozenset(S & out)))
    assert all(ch for ch, _ in faults), "distance 1"
    N = max(max(S) for S, _ in gates) + 1
    any_bad, per_out = 0, [0] * k
    for (c1, o1), (c2, o2) in combinations(faults, 2):
        if c1 == c2 and o1 != o2:
            any_bad += 1
            for q in o1 ^ o2:
                per_out[q] += 1                         # Z on output q survives
    return dict(n=len(faults), N=N, c=any_bad, per_output=per_out)


def metrics(n, N, c, p):
    acc = 1 - n * p
    return dict(acceptance=acc, error=c * p * p, time_rate=acc / n, volume_rate=acc / (n * N))


if __name__ == "__main__":
    rows = list(csv.DictReader(open(os.path.join(ROOT, "outputs", "smallest_per_class.csv"))))
    fig1 = json.load(open(os.path.join(ROOT, "outputs", "error_coefficients.json")))
    print("Ours vs the best prior construction (single level, d=2). Error is c p^2 for the whole target state.")
    print(f"{'target':20s} {'n':>4s} {'N':>3s} {'c':>4s} | {'prior n':>7s} {'c':>4s} | "
          + " | ".join(f"p={p:g}: err ours/prior, time-rate ours/prior" for p in P_VALUES))
    out = []
    for r in fig1:
        res = dict(target=r["target"], l=r["l"], n=r["n"], c=r["c"], prior=r["prior"], prior_n=r["prior_n"],
                   prior_c=r["prior_c"])
        line = f"{r['target']:20s} {r['n']:4d}"
        # our factory's footprint N from the catalogue row with this n and target
        N = None
        for x in rows:
            if x["l"] == str(r["l"]) and int(x["n"]) == r["n"] and x["k"] == str(r["k"]):
                try:
                    g = gates_of(x)
                except Exception:
                    continue
                st = fault_stats(g, r["k"], r["l"])
                if st["c"] == r["c"]:
                    N = st["N"]; res["per_output"] = st["per_output"]; break
        res["N"] = N
        line += f" {N if N else '?':>3} {r['c']:4d} | {r['prior_n']:7d} {r['prior_c']:4d} |"
        for p in P_VALUES:
            m = metrics(r["n"], N or 1, r["c"], p)
            mp = metrics(r["prior_n"], 1, r["prior_c"], p)
            res[f"p={p:g}"] = dict(ours=m, prior_error=mp["error"], prior_time_rate=mp["time_rate"])
            line += f" {m['error']:.2e}/{mp['error']:.2e}, {m['time_rate']:.4f}/{mp['time_rate']:.4f} |"
        out.append(res)
        print(line, flush=True)
    json.dump(out, open(os.path.join(ROOT, "outputs", "benchmark_factories.json"), "w"), indent=1)
