"""Check that the cascade recovers every factory found by the earlier separate runs.

Reference factories (each mapped to its exact Clifford class with the enumerate_targets orbit table):
  - stages (i)-(ii): outputs/all_factories.csv and outputs/large_unified_factories.csv (d=2, genuine k)
  - SAT, 'after' runs on the cluster: ~/msf_scratch/sat_results/after_*.json with result 'sat'
  - entangled SAT (Distillation_factories/distance_2/entangled_sat), via its union_catalogue.csv rows
For every reference class, the cascade (outputs/cascade/l{L}_k{k}.json) must have a factory with n no larger.
Prints every miss; exits 1 if there is any.

    python searches/check_cascade.py --cells 3:3-5
"""
import os, sys, csv, json, glob, argparse
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search  # noqa: E402,F401
import merge_large_unified as M  # noqa: E402
from classify import phase_polynomial  # noqa: E402
from coset_certify import parse_target  # noqa: E402
from enumerate_targets import orbit_table, index_of_poly  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SAT_DIR = os.path.expanduser("~/msf_scratch/sat_results")
ENT = os.path.expanduser("~/Library/CloudStorage/OneDrive-IBM/Documents/GitHub/Distillation_factories/union_catalogue.csv")
NAME = {"T": "T", "CS": "CS", "CCZ": "CCZ"}
NMAX_REF = 99                                             # set by --nmax: compare only factories the sweep can reach
SKIPPED = 0


def ent_spec(s):
    """'T1 CS12 CCZ123' (1-based) -> 'T0,CS01,CCZ012'."""
    out = []
    for tok in s.split():
        name = tok.rstrip("0123456789")
        out.append(name + "".join(str(int(c) - 1) for c in tok[len(name):]))
    return ",".join(out)


def labeller(L, k):
    """Exact Clifford-class label for any output polynomial at (L, k), by the method the cascade uses there."""
    try:
        C, off, bits, lab = orbit_table(L, k, log=lambda *x: None)
        return lambda P: str(int(lab[index_of_poly(P, k, L, C, off)]))
    except ValueError:
        if L == 3:
            from enumerate_targets_fiber import enumerate_targets_fiber
            _, F = enumerate_targets_fiber(L, k, log=lambda *x: None)
        else:
            from enumerate_targets_2adic import enumerate_targets_2adic
            _, F = enumerate_targets_2adic(L, k, log=lambda *x: None)
        def f(P):
            lb = F.label(F.digits_of_poly(P))
            return lb if isinstance(lb, str) else f"{lb[0]}:{lb[1]}"
        return f


def references(L, k, lab_of):
    global SKIPPED
    ref = defaultdict(list)                               # orbit -> [(n, source)]
    for path in ("all_factories.csv", "large_unified_factories.csv"):
        for r in csv.DictReader(open(os.path.join(ROOT, "outputs", path))):
            if int(r["l"]) != L or int(r["k"]) != k or r["d"] != "2" or r["search_stage"] in ("padding", "SAT"):
                continue
            if r["essential_dim"] != r["k"]:
                continue                                  # padded copy of a smaller-k factory: not a k-output class
            if int(r["N"]) > NMAX_REF:
                SKIPPED += 1
                continue                                  # beyond the range of the sweep being checked (--nmax)
            g = M.regenerate_circuit(L, k, int(r["N"]), r["params"], r["search_stage"], r["raw"])
            P = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), L)
            ref[lab_of(P)].append((int(r["n"]), f"{r['search_stage']} ({path})"))
    for f in glob.glob(os.path.join(SAT_DIR, "after_*.json")):
        d = json.load(open(f))
        if d["l"] == L and d["k"] == k and d["result"] == "sat":
            ref[lab_of(parse_target(d["spec"], L))].append((d["smaller_n"], f"long SAT {d['label']} N={d['N']}"))
    if L == 3 and os.path.exists(ENT):
        for r in csv.DictReader(open(ENT)):                # columns by name: the file's layout has changed before
            if (r["level"] == "3" and int(r["k"]) == k and r["d"] == "2" and "entangled_sat" in r["source"]
                    and r.get("gate")):
                try:
                    orb = lab_of(parse_target(ent_spec(r["gate"]), 3))
                except (ValueError, AssertionError, KeyError):
                    continue
                ref[orb].append((int(r["N"]), f"entangled SAT {r['gate']}"))  # their N = our n (inputs)
    return ref


def true_size(spec, L, k):
    """A reference output that a CNOT frame reduces to fewer outputs is a padded copy of a smaller-k class:
    return (k', polynomial on k' outputs) for its Clifford-reduced form."""
    import unified_search as us
    P = {T: c for T, c in parse_target(spec, L).items() if c}
    ess, deg, red, exact = us.reduce_frame(P, k, L)
    qs = sorted({q for S in red for q in S}); idx = {q: i for i, q in enumerate(qs)}
    return len(qs), {frozenset(idx[q] for q in S): c for S, c in red.items()}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", nargs="+", default=["3:3-5"])
    ap.add_argument("--cascade-dir", default=os.path.join(ROOT, "outputs", "cascade"))
    ap.add_argument("--nmax", type=int, default=99, help="circuit-size range of the sweep (main 11, larger 18); "
                    "catalogue factories on more qubits are outside it and reported as skipped")
    a = ap.parse_args()
    NMAX_REF = a.nmax
    misses = 0
    for c in a.cells:
        L, ks = c.split(":"); L = int(L)
        lo, hi = map(int, ks.split("-"))
        for k in range(lo, hi + 1):
            lab_of = labeller(L, k)
            cas = json.load(open(os.path.join(a.cascade_dir, f"l{L}_k{k}.json")))["summary"]["best"]
            ref = references(L, k, lab_of)
            bad, padded = [], 0
            for orb, lst in ref.items():
                n_ref, src = min(lst)
                got = cas.get(orb)
                if got is None and src.startswith(("entangled SAT", "long SAT")):
                    # not a genuine k-output class: compare at its true size in that cell
                    spec = (ent_spec(src.split("SAT ", 1)[1]) if src.startswith("entangled")
                            else json.load(open(glob.glob(os.path.join(SAT_DIR, f"after_{src.split()[2]}_N*.json"))[0]))["spec"])
                    kk, Pk = true_size(spec, L, k)
                    if kk < k:
                        ck = json.load(open(os.path.join(a.cascade_dir, f"l{L}_k{kk}.json")))["summary"]["best"]
                        gk = ck.get(labeller(L, kk)(Pk))
                        if gk is not None and gk["n"] <= n_ref:
                            padded += 1
                            continue
                        bad.append((orb, n_ref, src + f" [padded, true k={kk}]", gk["n"] if gk else None))
                        continue
                if got is None or got["n"] > n_ref:
                    bad.append((orb, n_ref, src, got["n"] if got else None))
            print(f"l={L} k={k}: {len(ref)} reference classes, {len(ref) - len(bad)} recovered by the cascade"
                  + (f" ({padded} padded, recovered at their true k)" if padded else "")
                  + ("" if not bad else f", {len(bad)} MISSED"), flush=True)
            for orb, n_ref, src, got in bad:
                print(f"   class {orb}: reference n={n_ref} ({src}), cascade {got}", flush=True)
            misses += len(bad)
    print(f"catalogue factories outside the sweep's circuit range (N > {NMAX_REF}), not compared: {SKIPPED}")
    sys.exit(1 if misses else 0)
