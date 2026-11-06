"""The lambda(U) route at l=3, k <= 6: distil T (or CS, CCZ) states with the distance-2 factories of the catalogue, then
convert them into the target with Clifford gates, Pauli measurements and feedforward, no catalysts; against our factory.

lambda(U) = fewest |T> states from which stabilizer operations (Cliffords, Pauli measurements, feedforward) make |U>,
without catalysts. lambda(U) <= tau(U) (CNOT+T synthesis). What we use from Beverland, Campbell, Howard, Kliuchnikov,
arXiv:1904.01124 (full text read):
  (a) stabilizer nullity nu is monotone under stabilizer operations, measurements included, and additive (Sec. 2.1,
      Prop. 4.1), so lambda(U) >= nu(U). We compute nu exactly (Pauli spectrum, searches/full_clifford_classes.py).
  (b) stabilizer extent xi is monotone under stabilizer operations (Sec. 1.1, 2.2) and multiplicative on products of
      states of at most three qubits (Lemma 2.7(1)). Hence m T states make |psi> only if
      m >= log xi(psi) / log xi(T); Table 1 gives this ratio as 2.96818 for CS and 3.63356 for CCZ (W3: 3.63356,
      W4: 4.99907, W5: 5.93637). For a product of T, CS, CCZ on disjoint qubits, xi multiplies (each factor has <= 3
      qubits), so lambda >= ceil(a + 2.96818 b + 3.63356 c). The asterisk of Table 1 ("not known to hold for arbitrary
      catalysts") does not apply: we exclude catalysts, and the inputs T^m are products of 1-qubit states. For W4, W5 the
      extent itself is "computed to eight digits and reconstructed", not proved (App. A.4); those rows are flagged.
  (c) upper bounds without catalysts: tau(U) (exact T-count, staged_best 'tau'), minimized over the full-Clifford class
      (outputs/checks/full_clifford_classes.json): Clifford gates are free on this route. 4 T -> CCZ (Jones, [32] of
      Table 1) is tau of CS01,CS02, a member of the CCZ class, so it is covered by this minimum. For products of
      T, CS, CCZ on disjoint qubits also a + 3b + 4c (one T, three T for CS, four for CCZ, each without catalyst).
lambda is called known for a class when the best lower bound equals the best upper bound. Route cost at distance 2: the
fewest inputs for lambda T states from the cheapest distance-2 T^m factories of the catalogue
(outputs/large_unified_factories.csv, rows d=2 whose output is a product of T gates), batched; for products, also any
mix of directly distilled CS and CCZ states (pure CS^m, CCZ^m factories of the same catalogue) and T states for the
rest. To first order the error of the route is the sum over the distilled states, so it is a distance-2 route.
T supply, strongest distance-2 alternative we know besides the catalogue: Meier, Eastin, Knill, arXiv:1204.4221 (full
text read, Sec. II-III): the four-qubit-code routine takes ten |H> states to
two, accepts with a(p) = 1 - 10p + O(p^2), output error u(p) = 9p^2 + O(p^3) per output (every single error detected),
no catalyst. |H> = cos(pi/8)|0> + sin(pi/8)|1> is Clifford-equivalent to |T> = T|+> (checked numerically below), so it
is a T^2 factory with n = 10 for this route. Every T count below is the cheapest batching of the catalogue T^m factories
and MEK (10 inputs per 2 T states, surplus discarded). The CS and CCZ pieces use the cheaper of the catalogue pure
factories and T states (3 T per CS, 4 T per CCZ) from that supply.
Added 2026-11-06 (referee round 4), from Balada Gaggioli and Marecek, arXiv:2610.31115 (full text
read): (i) Prop. 7 and Table I, a nine-input distance-2 factory for one T
state (acceptance 1/128, error 27/2 p^2), and Prop. 5, no stabilizer protocol does it with eight, so T^1 costs 9 here
per run; per accepted output it costs more (that paper, Discussion: "per accepted output, the nine-input T factory
consumes more inputs than the ten-input protocol"), while the other factories used accept with 1 - n p + O(p^2). The
tallies are therefore also reported without it (tally_without_BGM, tally_full_clifford_without_BGM);
(ii) Prop. 10 (= Beverland et al. Prop. A.4): projecting one qubit of |CCZ> onto (|0> +- i|1>)/sqrt2 leaves CS|++> up
to a Clifford, so m CS states cost at most what m CCZ states cost (8 for one, [[8,3,2]]; Table I: n_min(CS,2) = 8).
Both are measurement routes without catalysts and are used for the T and CS pieces below (checked: cs_from_ccz_check);
the CS and CCZ pieces may also be drawn together from one batch of CCZ states (CCZ^m factories of the catalogue).
Full-Clifford groups are read from outputs/checks/full_clifford_classes.json; member labels are 'k:orbit:spec' and the
k=6 orbit labels contain colons ('6:0:1536:CCZ045,CCZ123'), so the orbit is everything between the first and last colon.
Invariants reported for the full-Clifford classes (all exact Clifford invariants, see invariants()): Pauli spectrum
(multiset of |<P>|), the number of Lagrangians inside the zero set, and the multiset of per-label Lagrangian counts.
Also recounted for every class: the nullity floor LB(nu) of outputs/nullity_l3.json (cheapest T/CS/CCZ supply of total
nullity >= nu; a lower bound for every such route, measurements included).

    python searches/lambda_route.py      -> outputs/checks/lambda_route.{json,md}
"""
import os, sys, json, csv, re, math, collections, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
OUT = os.path.join(ROOT, "outputs", "checks")
SB = {"paper_v4": "outputs/archive/v4_published/outputs/staged_best.json", "current": "outputs/staged_best.json"}
RATIO = {"T": 1.0, "CS": 2.96818, "CCZ": 3.63356}                 # Table 1, lower-bound column, log xi / log xi(T)
W_RATIO = {3: 3.63356, 4: 4.99907, 5: 5.93637}                     # Table 1 rows W3, W4, W5
MEK = {2: 10}                                                      # Meier-Eastin-Knill 10 -> 2 |H>, distance 2
BGM = {1: 9}                                                       # Balada Gaggioli-Marecek Prop. 7: 9 -> 1 T, distance 2


def cs_from_ccz_check():
    """<c| = (<0| - i<1|)/sqrt2 on qubit 3 of |CCZ> leaves amplitudes (1-i)(1,1,1,i) on |ab>: CS|++> (Prop. A.4)."""
    amp = {(a, b): 1 + (-1j) * (-1) ** (a * b) for a in (0, 1) for b in (0, 1)}
    v = [amp[(a, b)] / amp[(0, 0)] for a in (0, 1) for b in (0, 1)]
    assert max(abs(x - y) for x, y in zip(v, [1, 1, 1, 1j])) < 1e-12
    return True


def mek_check():
    """|H> = cos(pi/8)|0> + sin(pi/8)|1> equals, up to a global phase, a single-qubit Clifford applied to T|+>."""
    import numpy as np
    H = np.array([[1, 1], [1, -1]]) / np.sqrt(2); S = np.diag([1, 1j]); I = np.eye(2)
    Tp = np.array([1, np.exp(1j * np.pi / 4)]) / np.sqrt(2); h = np.array([np.cos(np.pi / 8), np.sin(np.pi / 8)])
    cl = [I]
    for _ in range(4):                                       # words in H, S of length <= 4 cover the 24 Cliffords mod phase
        cl = cl + [g @ c for c in cl for g in (H, S)]
    hit = [c for c in cl if abs(abs(np.vdot(h, c @ Tp)) - 1) < 1e-12]
    assert hit, "|H> not Clifford-equivalent to |T>"
    return True


def catalogue_pure():
    """{gate: {m: smallest n}} for distance-2 factories whose output is m disjoint copies of one gate (T, CS, CCZ)."""
    best = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(ROOT, "outputs", "large_unified_factories.csv"))):
        if r["l"] != "3" or r["d"] != "2":
            continue
        toks = re.findall(r"\\mathrm\{([A-Za-z]+)\}_\{([0-9,]+)\}", r["reduced"])
        if not toks or len({t[0] for t in toks}) != 1:
            continue
        g = toks[0][0]; size = {"T": 1, "CS": 2, "CCZ": 3}.get(g)
        qs = [q for t in toks for q in (t[1].split(",") if "," in t[1] else list(t[1]))]
        if size is None or len(qs) != len(set(qs)) or len(qs) != int(r["k"]) or len(toks) * size != int(r["k"]):
            continue
        m = len(toks); n = int(r["n"])
        best[g][m] = min(best[g].get(m, 10 ** 9), n)
    return {g: dict(sorted(v.items())) for g, v in best.items()}


def knap(m, fac):
    """Fewest inputs for m copies from batches of the factories fac {outputs: n} (surplus outputs discarded)."""
    if m == 0:
        return 0
    c = [0] + [10 ** 9] * m
    for i in range(1, m + 1):
        for j, n in fac.items():
            c[i] = min(c[i], c[max(0, i - j)] + n)
    return c[m]


def product_specs(kmax=6):
    """(a, b, c, spec) for T^a CS^b CCZ^c on disjoint qubits, a + 2b + 3c = k <= kmax, at least one gate."""
    out = []
    for c in range(0, kmax // 3 + 1):
        for b in range(0, (kmax - 3 * c) // 2 + 1):
            for a in range(0, kmax - 3 * c - 2 * b + 1):
                k = a + 2 * b + 3 * c
                if k == 0:
                    continue
                q = iter(range(k)); g = []
                g += ["CCZ" + "".join(str(next(q)) for _ in range(3)) for _ in range(c)]
                g += ["CS" + "".join(str(next(q)) for _ in range(2)) for _ in range(b)]
                g += ["T" + str(next(q)) for _ in range(a)]
                out.append((a, b, c, k, ",".join(g)))
    return out


def w_poly(n):
    """W_n of Beverland et al. Eq. (4): phase pi/4 * ((x_1 xor ... xor x_n) + sum x_i), as {support: coeff mod 8}."""
    from unified_search import poly_from_table
    tab = [((bin(x).count("1") & 1) + bin(x).count("1")) % 8 for x in range(1 << n)]
    return poly_from_table(tab, n, 8)


def extent_floor(ratio, fac):
    """Fewest inputs for any supply of distilled T, CS, CCZ states (pure factories of fac) whose extent ratio
    t + 2.96818 s + 3.63356 c reaches ratio = log xi(U) / log xi(T). xi is monotone under stabilizer operations and
    multiplicative on the supply (a product of <= 3-qubit states, Lemma 2.7(1)), so every route from such a supply,
    measurements included, no catalysts, needs at least this many inputs."""
    best = None
    for c in range(0, 4):
        for s_ in range(0, 5):
            t = max(0, math.ceil(ratio - RATIO["CS"] * s_ - RATIO["CCZ"] * c - 1e-9))
            if c > max(fac["CCZ"]) * 3 or s_ > max(fac["CS"]) * 3:
                continue
            v = knap(t, fac["T"]) + knap(s_, fac["CS"]) + knap(c, fac["CCZ"])
            best = v if best is None else min(best, v)
    return best


def main():
    import full_clifford_classes as FC
    fcc = json.load(open(os.path.join(OUT, "full_clifford_classes.json")))
    per = {(r["k"], r["orbit"]): r for r in fcc["per_class"]}
    group = {}
    for r in per.values():
        group[(r["k"], r["orbit"])] = [(r["k"], r["orbit"])]
    for g in fcc["groups"]:
        mem = [(int(m.split(":")[0]), ":".join(m.split(":")[1:-1])) for m in g]
        assert all(x in per for x in mem), mem
        for x in mem:
            group[x] = mem
    rc_ = fcc["recount"]["current"]                         # known answers of the full-Clifford count (current table)
    assert len({tuple(sorted(v)) for v in group.values()}) == fcc["full_clifford_classes"] == 622
    assert rc_["full_counts"] == {"same": 366, "fewer": 256} and rc_["diagonal_wins_that_survive"] == 260 \
        and rc_["diagonal_wins"] == 330, rc_
    rows = {nm: {(r["k"], str(r["orbit"])): r for r in json.load(open(os.path.join(ROOT, p)))["rows"] if r["l"] == 3}
            for nm, p in SB.items()}
    floor = {(r["k"], r["spec"]): r for r in json.load(open(os.path.join(ROOT, "outputs", "nullity_l3.json")))}
    fac = catalogue_pure()
    fac["T_catalogue"] = dict(fac["T"]); fac["T_MEK"] = dict(MEK)
    fac["T_BGM"] = dict(BGM); fac["CS_catalogue"] = dict(fac["CS"])
    fac["T"] = {m: min(fac["T"].get(m, 10 ** 9), MEK.get(m, 10 ** 9), BGM.get(m, 10 ** 9))
                for m in set(fac["T"]) | set(MEK) | set(BGM)}
    fac["CS"] = {m: min(fac["CS"].get(m, 10 ** 9), fac["CCZ"].get(m, 10 ** 9)) for m in set(fac["CS"]) | set(fac["CCZ"])}
    assert mek_check() and cs_from_ccz_check()
    assert fac["T"][1] == 9 and fac["T"][2] == 10 and fac["CS"][1] == 8 and fac["CCZ"][1] == 8
    Tcat = lambda m: knap(m, fac["T_catalogue"])
    nT = lambda m: knap(m, fac["T"])
    fT0 = {m: min(fac["T_catalogue"].get(m, 10 ** 9), MEK.get(m, 10 ** 9)) for m in set(fac["T_catalogue"]) | set(MEK)}
    nT0 = lambda m: knap(m, fT0)                             # T supply without the nine-input BGM factory
    # lower bounds
    lower = {x: (per[x]["nullity"], "nullity") for x in per}
    note = collections.defaultdict(list)
    prod_of = {}; ratio_of = {}
    for a, b, c, k, spec in product_specs():
        _, label = FC.labeller(k)
        o = label(FC.parse(spec))
        if (k, o) not in per:
            continue                                         # (not genuine on k qubits cannot happen for a product)
        rt = a * RATIO["T"] + b * RATIO["CS"] + c * RATIO["CCZ"]; lb = math.ceil(rt - 1e-9)
        for x in group[(k, o)]:
            prod_of[x] = (a, b, c, spec); ratio_of[x] = (rt, "product")
            if lb > lower[x][0]:
                lower[x] = (lb, f"extent, product T^{a} CS^{b} CCZ^{c} (Lemma 2.7(1), Table 1 ratios)")
    for n in (3, 4, 5):
        _, label = FC.labeller(n); o = label(w_poly(n))
        lb = math.ceil(W_RATIO[n] - 1e-9)
        for x in group.get((n, o), []):
            note[x].append(f"W{n}"); ratio_of[x] = (W_RATIO[n], f"W{n} (numerical extent)")
            if lb > lower[x][0]:
                lower[x] = (lb, f"extent of W{n} (Table 1; extent value numerical, App. A.4)")
    # upper bounds: min tau over the full-Clifford class
    upper = {x: min(rows["current"][y]["tau"] for y in group[x]) for x in per}
    upper_from = {x: "tau (min over full-Clifford class)" for x in per}
    for x, (a, b, c, _) in prod_of.items():                 # pieces: T 1, CS 3 (tau), CCZ 4 (Jones, measurement)
        if a + 3 * b + 4 * c < upper[x]:
            upper[x] = a + 3 * b + 4 * c; upper_from[x] = "product of pieces: T 1, CS 3, CCZ 4 (Jones)"
    known = {x for x in per if lower[x][0] == upper[x]}
    assert all(lower[x][0] <= upper[x] for x in per), "a lower bound exceeds an upper bound"
    # known answers
    ccz = next(x for x in per if per[x]["spec"] == "CCZ012"); cs = next(x for x in per if per[x]["spec"] == "CS01")
    t = next(x for x in per if per[x]["spec"] == "T0")
    ka = {"lambda(CCZ)=4": ccz in known and upper[ccz] == 4, "lambda(CS)=3": cs in known and upper[cs] == 3,
          "lambda(T)=1": t in known and upper[t] == 1}
    assert all(ka.values()), ka
    out_rows = []; tally = {nm: collections.Counter() for nm in SB}; tally0 = {nm: collections.Counter() for nm in SB}
    for x in sorted(known):
        lam = upper[x]
        routes = {"lambda T states (catalogue + MEK)": nT(lam)}
        rc0 = nT0(lam)
        routes_cat = {"lambda T states (catalogue only)": Tcat(lam)}
        if x in prod_of:
            a, b, c, _ = prod_of[x]
            best = None; bestc = None
            for bs, cs_ in itertools.product(range(b + 1), range(c + 1)):   # bs CS and cs_ CCZ distilled directly
                m = a + 3 * (b - bs) + 4 * (c - cs_)
                v = nT(m) + min(knap(bs, fac.get("CS", {})) + knap(cs_, fac.get("CCZ", {})),
                                knap(bs + cs_, fac.get("CCZ", {})))   # CS from a pooled CCZ batch by a Y measurement
                vc = Tcat(m) + knap(bs, fac["CS_catalogue"]) + knap(cs_, fac.get("CCZ", {}))
                v0 = nT0(m) + min(knap(bs, fac.get("CS", {})) + knap(cs_, fac.get("CCZ", {})),
                                  knap(bs + cs_, fac.get("CCZ", {})))
                rc0 = min(rc0, v0)
                best = v if best is None else min(best, v); bestc = vc if bestc is None else min(bestc, vc)
            routes["pieces: direct CS/CCZ + T states (catalogue + MEK) for the rest"] = best
            routes_cat["pieces, catalogue only"] = bestc
        rc = min(routes.values()); rc_cat = min(routes_cat.values())
        ef = extent_floor(ratio_of[x][0], fac) if x in ratio_of else None
        for nm in SB:
            r0 = rows[nm][x]
            tally0[nm]["factory fewer" if r0["n"] < rc0 else "same" if r0["n"] == rc0 else "route fewer"] += 1
        for nm in SB:
            r = rows[nm][x]
            res = "factory fewer" if r["n"] < rc else "same" if r["n"] == rc else "route fewer"
            tally[nm][res] += 1
            if nm == "paper_v4":
                out_rows.append(dict(k=x[0], orbit=x[1], spec=per[x]["spec"], nu=per[x]["nullity"], tau=r["tau"],
                                     lambda_=lam, lambda_lower_from=lower[x][1], lambda_upper_from=upper_from[x], W=note.get(x, []),
                                     product=prod_of.get(x, (None,) * 4)[3], routes=routes, route_n=rc,
                                     route_n_catalogue_only=rc_cat, route_n_without_BGM=rc0, mek_lowers_route=rc < rc_cat,
                                     group=[f"{y[0]}:{y[1]}:{per[y]['spec']}" for y in group[x]],
                                     extent_floor_any_TCSCCZ_supply=ef, extent_ratio=ratio_of.get(x),
                                     factory_below_extent_floor=(ef is not None and r["n"] < ef),
                                     factory_n_paper=r["n"], factory_n_current=rows["current"][x]["n"],
                                     result_paper=res,
                                     result_current="factory fewer" if rows["current"][x]["n"] < rc else
                                     "same" if rows["current"][x]["n"] == rc else "route fewer"))
    fl = {nm: collections.Counter() for nm in SB}
    for nm in SB:
        for x, r in rows[nm].items():
            lb = floor[(x[0], r["spec"])]["route_lower_bound"]
            fl[nm]["below floor" if r["n"] < lb else "at floor" if r["n"] == lb else "above floor (open)"] += 1
    gap = collections.Counter(upper[x] - lower[x][0] for x in per)
    inv = invariants(fcc, per, group)
    # route vs factory per full-Clifford class (lambda known for one member <=> for all: lambda is a Clifford invariant)
    fc_groups = {tuple(sorted(group[x])) for x in known}
    fc_tally = {nm: collections.Counter() for nm in SB}; fc_tally0 = {nm: collections.Counter() for nm in SB}
    fc_syn = {nm: {"with_BGM": collections.Counter(), "without_BGM": collections.Counter()} for nm in SB}
    for G in fc_groups:
        rr = [r for r in out_rows if (r["k"], r["orbit"]) in G]
        rc = min(r["route_n"] for r in rr); rc0 = min(r["route_n_without_BGM"] for r in rr)
        for nm in SB:
            n = min(rows[nm][y]["n"] for y in G)
            fc_tally[nm]["factory fewer" if n < rc else "same" if n == rc else "route fewer"] += 1
            fc_tally0[nm]["factory fewer" if n < rc0 else "same" if n == rc0 else "route fewer"] += 1
            sy = [rows[nm][y].get("syn_n") for y in G if rows[nm][y].get("syn_n") is not None]
            sy = min(sy) if sy else None                         # cheapest synthillation construction over the class
            for key, v in (("with_BGM", rc), ("without_BGM", rc0)):
                fc_syn[nm][key]["no construction" if sy is None else "route fewer than synthillation" if v < sy
                                else "same" if v == sy else "synthillation fewer"] += 1
    res = dict(doc=__doc__, known_answers=ka, mek_H_clifford_equivalent_T=True, invariants=inv,
               lambda_known_full_clifford_classes=len(fc_groups), of_full_clifford_classes=len({tuple(sorted(v)) for v in group.values()}),
               tally_full_clifford=({nm: dict(sorted(v.items())) for nm, v in fc_tally.items()}),
               tally_full_clifford_without_BGM=({nm: dict(sorted(v.items())) for nm, v in fc_tally0.items()}),
               route_vs_synthillation_full_clifford={nm: {k: dict(sorted(c.items())) for k, c in v.items()}
                                                     for nm, v in fc_syn.items()},
               tally_without_BGM={nm: dict(sorted(v.items())) for nm, v in tally0.items()},
               route_lowered_by_MEK=[r["spec"] for r in out_rows if r["mek_lowers_route"]], catalogue_pure_factories=fac,
               lambda_known=len(known), of=len(per), lambda_known_by_k=dict(sorted(collections.Counter(x[0] for x in known).items())),
               lambda_known_reason=dict(sorted(collections.Counter(lower[x][1].split(",")[0] for x in known).items())),
               lambda_gap_histogram=dict(sorted(gap.items())), rows=out_rows, tally={nm: dict(sorted(v.items())) for nm, v in tally.items()},
               nullity_floor={nm: dict(sorted(v.items())) for nm, v in fl.items()},
               per_class=[dict(k=x[0], orbit=x[1], spec=per[x]["spec"], lambda_lower=lower[x][0], lower_from=lower[x][1],
                               lambda_upper=upper[x], upper_from=upper_from[x]) for x in sorted(per)],
               sources=dict(staged_best=SB, full_clifford="outputs/checks/full_clifford_classes.json",
                            nullity_floor="outputs/nullity_l3.json", factories="outputs/large_unified_factories.csv",
                            literature="Beverland et al. arXiv:1904.01124, Sec. 2.1-2.2, Lemma 2.7, Prop. 4.1, Tables 1, 4; "
                                       "Meier, Eastin, Knill arXiv:1204.4221, Sec. II-III (10 -> 2, a(p), u(p)); "
                                       "Balada Gaggioli, Marecek arXiv:2610.31115, Props. 5, 7, 10, Table I"))
    json.dump(res, open(os.path.join(OUT, "lambda_route.json"), "w"), indent=1)
    write_md(res)
    print(json.dumps({k: v for k, v in res.items() if k not in ("doc", "rows", "per_class")}, indent=1))


def invariants(fcc, per, group):
    """Clifford invariants of each diagonal class and the number of distinct values among the full-Clifford classes.
    For a Clifford V: (i) |<V psi|V P V^dag|V psi>| = |<psi|P|psi>|, so the Pauli spectrum is invariant; (ii) V maps the
    zero set of psi onto that of V psi and Lagrangians onto Lagrangians (symplectic action), so the number of Lagrangians
    inside the zero set is invariant; (iii) for L in the zero set of psi and L' = V L V^dag, C_L' V C_L^dag maps <Z_i>
    to itself, so C_L' V psi and C_L psi lie in one diagonal class: the label of each Lagrangian is carried along, and
    the per-label Lagrangian counts (as a dict, hence also as a multiset of counts) are invariant. Checked on the data:
    every full-Clifford group has one value of each invariant."""
    G = {tuple(sorted(v)) for v in group.values()}
    f = {"spectrum": lambda r: (r["k"], json.dumps(r["spectrum"])),
         "spectrum + n_lagrangians": lambda r: (r["k"], json.dumps(r["spectrum"]), r["n_lagrangians"]),
         "spectrum + n_lagrangians + multiset of per-label counts":
             lambda r: (r["k"], json.dumps(r["spectrum"]), r["n_lagrangians"], tuple(sorted(r["labels"].values())))}
    out = {}
    for nm, fn in f.items():
        out[nm] = dict(groups_with_unequal_values=sum(len({fn(per[x]) for x in g}) > 1 for g in G),
                       distinct_values_over_full_classes=len({fn(per[g[0]]) for g in G}), full_clifford_classes=len(G))
        assert out[nm]["groups_with_unequal_values"] == 0, nm
    out["note"] = ("each count of distinct values is a lower bound on the number of full-Clifford classes; the exact "
                   "count (union of label sets, full_clifford_classes.json) is the number of groups")
    return out


def write_md(o):
    L = ["# The lambda(U) route (Hadamards, Pauli measurements, feedforward; no catalysts), l=3, k <= 6", "",
         "Script: `searches/lambda_route.py`. Data: `outputs/checks/lambda_route.json`. Needs "
         "`outputs/checks/full_clifford_classes.json` (searches/full_clifford_classes.py).", "",
         "## Method", "", "```", o["doc"].strip(), "```", "",
         "## Result", "",
         f"- Known answers reproduced: {o['known_answers']}.",
         f"- lambda is known (lower bound = upper bound) for **{o['lambda_known']} of {o['of']}** classes; by k: "
         f"{o['lambda_known_by_k']}; lower bound that closes it: {o['lambda_known_reason']}.",
         f"- Gap upper - lower over all {o['of']} classes (gap: classes): {o['lambda_gap_histogram']}.",
         f"- Distance-2 factories with pure outputs used for the route: {o['catalogue_pure_factories']}.", ""]
    nb = sum(r["factory_below_extent_floor"] for r in o["rows"])
    L += [f"- Extent floor over any supply of distilled T, CS, CCZ states (any conversion by stabilizer operations, no "
          f"catalysts; T supply includes MEK): our factory is below it for {nb} of {o['lambda_known']} of these classes"
          + (", so this monotone proves no factory optimal against the measurement route; it only bounds the route "
             "from below." if nb == 0 else "."),
          "- Rows where the route beats the factory (current table): "
          + (", ".join(f"{r['spec']} (route {r['route_n']}, factory {r['factory_n_current']}; full-Clifford group {r['group']})"
                       for r in o["rows"] if r["result_current"] == "route fewer") or "none") + ".",
          f"- Rows whose route cost the supplies outside the catalogue (MEK 10 -> 2 T, BGM 9 -> 1 T, CCZ + Y -> CS) lower: "
          f"{len(o['route_lowered_by_MEK'])} of {o['lambda_known']}: {o['route_lowered_by_MEK']}.",
          f"- Lower bound that closes lambda: {o['lambda_known_reason']} (W4, W5 extents are numerical, App. A.4 of "
          "Beverland et al.).",
          f"- Per full-Clifford class (lambda known for {o['lambda_known_full_clifford_classes']} of "
          f"{o['of_full_clifford_classes']}; route and factory each minimised over the members): {o['tally_full_clifford']}.",
          "", "## Clifford invariants of the full-Clifford classes", ""]
    for nm, v in o["invariants"].items():
        if nm != "note":
            L.append(f"- {nm}: {v['distinct_values_over_full_classes']} distinct values over {v['full_clifford_classes']} "
                     f"full-Clifford classes; groups with unequal values: {v['groups_with_unequal_values']}.")
    L += [f"- {o['invariants']['note']}.", ""]
    for nm, t in o["tally"].items():
        L.append(f"- Classes with lambda known, factory ({nm}) vs cheapest lambda route: {t} of {o['lambda_known']}.")
        L.append(f"- Same, without the nine-input T factory of arXiv:2610.31115 (acceptance 1/128): "
                 f"{o['tally_without_BGM'][nm]} of {o['lambda_known']}; per full-Clifford class: "
                 f"{o['tally_full_clifford_without_BGM'][nm]} of {o['lambda_known_full_clifford_classes']}.")
        L.append(f"- Route vs the cheapest synthillation construction over each full-Clifford class ({nm}): "
                 f"{o['route_vs_synthillation_full_clifford'][nm]} of {o['lambda_known_full_clifford_classes']}.")
    for nm, t in o["nullity_floor"].items():
        L.append(f"- All {o['of']} classes, factory ({nm}) vs nullity floor of every T/CS/CCZ route: {t}.")
    L += ["", "## Classes with lambda known", "",
          "| k | spec | nu | tau | lambda | lower bound from | lambda route n (catalogue only) | extent floor, any T/CS/CCZ supply | factory n (paper) | factory n (current) | factory vs lambda route (paper) | (current) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in o["rows"]:
        L.append(f"| {r['k']} | {r['spec']} | {r['nu']} | {r['tau']} | {r['lambda_']} | {r['lambda_lower_from']} | "
                 f"{r['route_n']} ({r['route_n_catalogue_only']}) | {r['extent_floor_any_TCSCCZ_supply']} | {r['factory_n_paper']} | {r['factory_n_current']} | {r['result_paper']} | {r['result_current']} |")
    open(os.path.join(OUT, "lambda_route.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
