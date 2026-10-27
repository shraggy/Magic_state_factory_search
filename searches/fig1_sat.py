"""SAT check of every Fig. 1 factory: is there a factory for the same target with fewer inputs on N or N+1 qubits?
UNSAT proves ours is the smallest on that many qubits; SAT returns a smaller one. Targets run in parallel; each result
is written to outputs/fig1_sat.json as it finishes. Run: python3 searches/fig1_sat.py [--time 3600] [--jobs 15]"""
import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, "outputs", "fig1_sat.json")
# (label, level, target spec on outputs 0..k-1, k, our input count n, our circuit size N)
FIG1 = [("CS12 T3", 3, "CS01,T2", 3, 20, 6), ("CS12 CS34 CS56", 3, "CS01,CS23,CS45", 6, 24, 8),
        ("CS12 CS34 T5", 3, "CS01,CS23,T4", 5, 26, 8), ("CCZ123 CS45 T6", 3, "CCZ012,CS34,T5", 6, 26, 9),
        ("CS12 T3 T4 T5", 3, "CS01,T2,T3,T4", 5, 30, 8), ("CS12 T3...T7", 3, "CS01,T2,T3,T4,T5,T6", 7, 36, 10),
        ("CT12 sT3", 4, "CT01,sT2", 3, 44, 7), ("CCS123 sT4", 4, "CCS012,sT3", 4, 48, 8),
        ("CT12 sT3 sT4", 4, "CT01,sT2,sT3", 4, 50, 8), ("CT12 CT34 sT5", 4, "CT01,CT23,sT4", 5, 68, 9),
        ("CCS123 CT45 sT6", 4, "CCS012,CT34,sT5", 6, 68, 10), ("CT12 CT34 sT5 sT6", 4, "CT01,CT23,sT4,sT5", 6, 70, 10)]


def job(label, L, spec, k, n_ours, N, tlimit):
    from sat_certify import certify
    from coset_certify import parse_target
    t = time.time()
    res, n_found, gates = certify(parse_target(spec, L), k, N, L, n_ours, tlimit)
    return dict(target=label, l=L, k=k, n=n_ours, N=N, result=res, smaller_n=n_found if res == "sat" else None,
                seconds=round(time.time() - t, 1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", type=float, default=3600.0, help="SAT time limit per (target, N), seconds")
    ap.add_argument("--jobs", type=int, default=15)
    ap.add_argument("--extra", type=int, default=1, help="also try N+1, ..., N+extra qubits")
    a = ap.parse_args()
    done = json.load(open(OUT)) if os.path.exists(OUT) else []
    have = {(d["target"], d["N"]) for d in done}
    todo = [(lab, L, sp, k, n, N + e) for lab, L, sp, k, n, N in FIG1 for e in range(a.extra + 1) if (lab, N + e) not in have]
    with ProcessPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(job, *t, a.time) for t in todo]
        for fu in as_completed(futs):
            d = fu.result(); done.append(d); json.dump(done, open(OUT, "w"), indent=1)
            print(f"{d['target']:18s} N={d['N']:2d}: {d['result']}" + (f" (n={d['smaller_n']})" if d['smaller_n'] else "") + f"  {d['seconds']}s", flush=True)
