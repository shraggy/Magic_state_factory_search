"""Tables SIV(a)-(c): representative smallest-n distance-2 factories at l=2, 3, 4, from the staged cascade
(outputs/staged_best.json, every factory re-verified there), with a Quirk link per row.

Rows: l=2 every class (k<=6); l=3 every class with k<=4, every class on the bound n = tau+1, for each gate-type mix at k=5,6, the largest saving over synthillation and the smallest n, and at k=7,8 one CCZ-layer factory per
distinct n of the cascade on a target list (outputs/k78_pipeline/full/l3_k{7,8}/ii_eg.json, iii_f.json); l=4 every class with k<=2 and the smallest and largest n at k=3 and k=4, over the classes that need a fourth-level
rotation: the 22 l=4 classes flagged third_level_target in staged_best.json (third-level targets, covered by the l=3
rows) are left out unless --with-third-level-targets is given.
Each row is re-verified here (unified_search.is_factory, brute force for N<=14, distance 2, n recomputed).
Stage: the first pass of the cascade whose best reached the final n (staged_best.json 'first').
ref.: an earlier paper that gives the factory (REFS, each read in the cited paper); 'syn' when no earlier paper gives
it but a synthillation construction has the same n (n = n_syn of outputs/synthillation_strongest_l3.json), or, at
k = 7, 8, when the factory is on the bound with one check, which by Corollary 3 (saturation) is the single-check form
of synthillation's case 11 (tau odd); 'new' otherwise (n < n_syn, dagger; or k = 8 off the bound). At l = 2, where every
output is a stabilizer state, a row without an earlier paper is marked '--'.
tau at l = 2 (S-count), where staged_best.json leaves it empty: the fewest rank-one terms y y^T summing to the symmetric
F2 matrix B of the output (S exponents on the diagonal, CZ pairs off it), rank(B) if diag(B) != 0 and rank(B) + 1
otherwise (Lempel 1975); asserted equal to the stored tau on every l = 2 class where it is stored.
    python searches/table_sv.py   -> outputs/checks/table_sv.json, Magic_State_Factories_tex/MSD_manuscript/table_sv.tex
                                     and .../circuits_tableVI.tex (one Quirk macro per row)
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
TEX = os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--sb", default="staged_best.json", help="staged-best file under outputs/")
_p.add_argument("--out-dir", default=None, help="write every output here instead of outputs/ (and the manuscript)")
_p.add_argument("--with-third-level-targets", action="store_true",
                help="keep the 22 l=4 classes whose target is a level-3 gate (default: left out)")
_A = _p.parse_args(); SB_FILE = _A.sb; OUTD = _A.out_dir; WITH_T3 = _A.with_third_level_targets

CODE = {"i_tg": "TG", "i_gs": "GS", "i_gh": "GS-c", "ii_dj": "GL-d", "ii_eg": "GL-e", "iii_r": "SAT-r", "iii_f": "SAT-f"}
# literature: (l, k, spec, n) -> bibtex keys; each was read in the cited paper
REFS = {(3, 1, "T0", 14): "webster2023transversal", (3, 2, "T0,T1", 14): "bravyi2012magic",   # (3,1): Webster et al. Example 4.1, Fig. 6
        (3, 2, "CS01", 12): "webster2023transversal,anqi2026qec",
        (3, 3, "CCZ012", 8): "bombin2007exact,jones2013composite,eastin2013distilling",
        (3, 3, "T0,T1,T2", 20): "nezami2022classification", (3, 4, "CS01,CS03,CS12", 18): "campbell2017unified",
        (3, 4, "T0,T1,T2,T3", 20): "bravyi2012magic", (3, 5, "T0,T1,T2,T3,T4", 26): "nezami2022classification",
        (3, 5, "CCZ014,CCZ023", 12): "campbell2017unified", (3, 6, "CCZ045,CCZ123", 14): "campbell2017unified",
        (2, 1, "S0", 6): "webster2023transversal",   # [[2^(l+1)-2,1,2]] at l=2: Webster et al. Example 4.1
        (2, 2, "CZ01", 4): "PhysRevA.54.4741", (2, 2, "S0,S1", 6): "jones2013multilevel",
        (2, 4, "CZ03,CZ12", 6): "PhysRevA.54.4741", (2, 6, "CZ05,CZ14,CZ23", 8): "PhysRevA.54.4741",   # Iceberg [[k+2,k,2]]
        (3, 4, "CS03,CS12", 18): "anqi2026qec", (3, 6, "CS05,CS14,CS23", 24): "anqi2026qec",            # (2m+2)CS -> mCS, m=2,3
        (4, 4, "CCCZ0123", 16): "11096910,webster2022xp",
        (4, 1, "sT0", 30): "webster2023transversal",          # Example 4.1 at level t = 4 (outputs/checks/webster_family_check.json)
        # l=3 factories with every coefficient doubled (each T rotation from two sqrt(T)); outputs/checks/doubled_l3_at_l4.json
        (4, 1, "sT0,sT0", 28): "x2:webster2023transversal", (4, 2, "CT01,CT01", 24): "x2:webster2023transversal,anqi2026qec",
        (4, 2, "sT0,sT0,sT1,sT1", 28): "x2:bravyi2012magic",
        (4, 3, "CCS012,CCS012", 16): "x2:bombin2007exact,jones2013composite,eastin2013distilling"}
TYPES = [("CS",), ("CCZ",), ("CCZ", "CS"), ("CCZ", "T"), ("CS", "T"), ("CCZ", "CS", "T")]
NAME = {2: {1: (r"\mathrm{S}",), 2: (r"\CZ",)}, 3: {1: (r"\mathrm{T}",), 2: (r"\CS",), 3: (r"\CCZ",)},
        4: {1: (r"\sqrt{\mathrm{T}}", r"\mathrm{T}"), 2: (r"\CT", r"\CS"), 3: (r"\mathrm{CCS}", r"\CCZ"),
            4: (r"\mathrm{CCCZ}",)}}


def terms(spec):
    t = collections.Counter()
    for tok in spec.split(","):
        nm = tok.rstrip("0123456789"); t[frozenset(int(c) for c in tok[len(nm):])] += 1
    return t


def gtypes(spec):
    return tuple(sorted({tok.rstrip("0123456789") for tok in spec.split(",")}))


def tex_output(spec, L, budget=24):
    groups = collections.defaultdict(list)
    for S, m in terms(spec).items():
        names = NAME[L][len(S)]
        if len(names) == 1:
            groups[(len(S), 0)].append(S)
        else:
            if m & 2: groups[(len(S), 1)].append(S)
            if m & 1: groups[(len(S), 0)].append(S)
    toks = [(NAME[L][s][i], ["".join(str(q + 1) for q in sorted(S)) for S in sorted(groups[(s, i)], key=sorted)])
            for (s, i) in sorted(groups, key=lambda x: (-x[0], x[1]))]
    w = lambda nm: 1.6 * len(nm.replace("\\mathrm", "").strip("{}\\"))
    lines, cur, width = [], [], 0.0
    for nm, subs in toks:
        part = []
        for sb in subs:
            add = (w(nm) if not part else 0) + len(sb) + 1
            if width + add > budget and (part or cur):
                if part: cur.append((nm, part))
                lines.append(cur); cur, part, width = [], [], 0.0; add = w(nm) + len(sb) + 1
            part.append(sb); width += add
        cur.append((nm, part))
    lines.append(cur)
    txt = ["$" + "".join(f"{nm}_{{{','.join(sb)}}}" for nm, sb in ln) + "$" for ln in lines if ln]
    return txt[0] if len(txt) == 1 else r"\begin{tabular}[t]{@{}l@{}}" + r" \\ ".join(txt) + r"\end{tabular}"


def verify(r):
    g = [(frozenset(S), c) for S, c in r["gates"]]; f = dict(k=r["k"], r=r["N"] - r["k"], gates=g)
    ok = us.is_factory(f, r["l"]) and (r["N"] > 14 or us.is_factory(f, r["l"], brute_force=True))
    return ok and us.distance(f, r["l"]) == 2 and sum(us.inputs(c, r["l"]) for _, c in g) == r["n"]


def select(SB):
    key = lambda r: (r["n"], r["N"], len(r["spec"].split(",")), r["spec"])
    out = collections.OrderedDict()
    def add(r, why): out.setdefault((r["l"], r["k"], r["orbit"]), (r, []))[1].append(why)
    for r in sorted([r for r in SB if r["l"] == 2], key=lambda r: (r["k"],) + key(r)): add(r, "l=2, every class")
    L3 = [r for r in SB if r["l"] == 3]
    for r in sorted([r for r in L3 if r["k"] <= 4], key=lambda r: (r["k"],) + key(r)): add(r, "k<=4")
    for r in L3:
        if r["n"] == r["tau"] + 1: add(r, "on the bound n = tau+1")
    for k in (5, 6):                                       # the largest saving over synthillation, one per gate-type mix
        for t in TYPES:
            T = [r for r in L3 if r["k"] == k and gtypes(r["spec"]) == t and r["n"] < r["syn_n"]]
            if T: add(min(T, key=lambda r: (r["n"] - r["syn_n"],) + key(r)), "largest saving over synthillation, type " + "+".join(t))
    for k in (5, 6):
        for t in TYPES:
            T = [r for r in L3 if r["k"] == k and gtypes(r["spec"]) == t]
            if T: add(min(T, key=key), "smallest n of type " + "+".join(t))
    L4 = [r for r in SB if r["l"] == 4 and (WITH_T3 or not r["third_level_target"])]
    for r in sorted([r for r in L4 if r["k"] <= 2], key=lambda r: (r["k"],) + key(r)): add(r, "l=4, k<=2")
    for k in (3, 4):
        C = [r for r in L4 if r["k"] == k]
        add(min(C, key=key), "smallest n at this k"); add(max(C, key=lambda r: (r["n"], -r["N"])), "largest n at this k")
    return [(r, why) for r, why in out.values()]


def proved_on_N():
    """(l, k, orbit) whose coefficient free-SAT query is UNSAT on at least the factory's own N (staged free_coef* files)."""
    out = set(); SB = json.load(open(os.path.join(O, SB_FILE)))
    best = {(r["l"], r["k"], r["orbit"]): r for r in SB["rows"]}
    for key, m in SB["cells"].items():
        L, k = int(key[1]), int(key.split("_k")[1])
        src = os.path.join(O, "staged_pc" if (L == 3 and k <= 4) else "staged_A", key) if L == 3 or key == "l4_k4" else os.path.join(ROOT, m["source"])
        for f in __import__("glob").glob(os.path.join(src, "free_coef*", "*.json")):   # the all-coefficient z3 proofs
            for o, v in json.load(open(f))["per_class"].items():
                mx = max([s["N"] for s in v.get("steps") or [] if s["last"] == "unsat"], default=0)
                if mx >= best[(L, k, o)]["N"]: out.add((L, k, o))
        f = os.path.join(ROOT, m["source"], "iii_f.json")              # one-process cells (l=2, l=4 k<=3): z3 free SAT records
        if L == 3: f = ""                                                # at l=3 iii_f is the XOR chain; proofs are the z3 files above
        if os.path.exists(f):
            for o, v in (json.load(open(f)).get("steps") or {}).items():
                st = v if isinstance(v, list) else v.get("steps")
                mx = max([x["N"] for x in st or [] if x.get("last") == "unsat"], default=0)
                if (L, k, o) in best and mx >= best[(L, k, o)]["N"]: out.add((L, k, o))
    return out


def tau_l2(spec, k):
    """S-count of an l=2 output: rank of its symmetric F2 matrix, +1 if the diagonal is zero (Lempel)."""
    rows = [0] * k
    for tok in spec.split(","):
        nm = tok.rstrip("0123456789"); q = [int(c) for c in tok[len(nm):]]
        if nm == "S": rows[q[0]] |= 1 << q[0]
        elif nm == "CZ": rows[q[0]] |= 1 << q[1]; rows[q[1]] |= 1 << q[0]
    diag = any(rows[i] >> i & 1 for i in range(k))
    B = []
    for v in rows:
        for b in B: v = min(v, v ^ b)
        if v: B.append(v)
    return len(B) if diag else len(B) + 1


def main():
    SB = json.load(open(os.path.join(O, SB_FILE)))["rows"]
    for r in SB:
        if r["l"] == 2:
            t = tau_l2(r["spec"], r["k"])
            assert r["tau"] is None or r["tau"] == t, (r["spec"], r["tau"], t)
            r["tau"] = t
    PROVED = proved_on_N()
    rows, macros, dropped = [], [], []
    for i, (r, why) in enumerate(select(SB)):
        if not verify(r):
            dropped.append(dict(l=r["l"], k=r["k"], orbit=r["orbit"])); continue
        name = "qSIII" + "".join(chr(97 + int(d)) for d in f"{i:03d}")
        g = [(frozenset(S), c) for S, c in r["gates"]]
        macros.append(r"\newcommand{\%s}{%s}" % (name, us.quirk_url(g, r["N"], r["l"], r["k"]).replace("%", r"\%").replace("#", r"\#")))
        rows.append(dict(r, macro=name, why=why, stage_code=CODE[r["first"]], ref=REFS.get((r["l"], r["k"], r["spec"], r["n"]), "new"),
                         proved=(r["l"], r["k"], r["orbit"]) in PROVED))
    for K in (7, 8):                     # k = 7, 8 CCZ layers beyond the class lists: per distinct n, the factory of the
        best = {}                        # cascade on the target list (outputs/k78_pipeline/full), tightest T-count first
        for p in ("ii_eg", "iii_f"):
            for b in json.load(open(os.path.join(O, "k78_pipeline", "full", f"l3_k{K}", f"{p}.json")))["best"].values():
                lo, hi = b["tau"]
                key = b["n"]
                if key not in best or (hi - lo, b["N"]) < (best[key]["tau"][1] - best[key]["tau"][0], best[key]["N"]):
                    best[key] = b
        for n_ in sorted(best)[:3]:                                   # the three smallest n: those on the bound first
            b = best[n_]; lo, hi = b["tau"]; exact = lo == hi
            g = [(frozenset(S), c) for S, c in b["gates"]]
            rr = dict(l=3, k=K, N=b["N"], n=b["n"], gates=b["gates"], spec=",".join("CCZ" + "".join(map(str, sorted(x))) for x in b["ccz"]))
            if not verify(rr):
                dropped.append(dict(l=3, k=K, spec="ccz layer")); continue
            name = "qSIIIk" + "".join(chr(97 + int(d)) for d in f"{len(macros):03d}")
            macros.append(r"\newcommand{\%s}{%s}" % (name, us.quirk_url(g, b["N"], 3, K).replace("%", r"\%").replace("#", r"\#")))
            st = b["stage"].split(":")[0]
            onb = exact and b["n"] == lo + 1 and b["N"] == K + 1 and lo % 2 == 1   # one check, on the bound: case 11 form
            rows.append(dict(rr, macro=name, why=["k=7, 8 CCZ layer"], stage_code={"B": CODE.get("i_gs", "GS"), "ii_eg": CODE.get("ii_eg", "GL-e")}.get(st, st),
                             ref="syn" if onb else "new", tau_text=str(lo) if exact else f"{lo}--{hi}", syn_n=None, block="k78",
                             tau=lo if exact else None))
    assert not dropped, dropped
    json.dump(dict(source=f"outputs/{SB_FILE}", rows=rows), open(os.path.join(OUTD or os.path.join(O, "checks"), "table_sv.json"), "w"), indent=0)
    open(os.path.join(OUTD or TEX, "circuits_tableVI.tex"), "w").write(
        "% Quirk links for Tables SIII(a-c), generated by searches/table_sv.py from outputs/staged_best.json\n" + "\n".join(macros) + "\n")
    open(os.path.join(OUTD or TEX, "table_sv.tex"), "w").write(tables(rows))
    print(f"{len(rows)} rows: " + ", ".join(f"l={L}: {sum(r['l'] == L for r in rows)}" for L in (2, 3, 4)))


def row(r, L, syn):
    bold = r["tau"] is not None and r["n"] == r["tau"] + 1
    nk = f"[[{r['n']},{r['k']},2]]_{{{r['N']}}}" + ("^\\circ" if r.get("proved") else "")
    tau = r.get("tau_text") or ("--" if r["tau"] is None else str(r["tau"]))
    if r["ref"] in ("new", "syn"):
        ref = r["ref"]
        if L == 2: ref = "--"
        elif ref == "new" and r.get("syn_n") and r["n"] == r["syn_n"]: ref = "syn"
    else:
        ref = (r"$\times2$\,\cite{%s}" % r["ref"][3:]) if r["ref"].startswith("x2:") else rf"\cite{{{r['ref']}}}"
    if syn and r.get("syn_n") and r["n"] < r["syn_n"]: ref += r"$^\dagger$"
    cols = [rf"$\mathbf{{{nk}}}$" if bold else f"${nk}$", rf"\quirk{{\{r['macro']}}}{{{tex_output(r['spec'], L)}}}",
            rf"\textbf{{{tau}}}" if bold else tau, r["stage_code"]]
    if syn: cols.append("--" if not r.get("syn_n") else str(r["syn_n"]))
    return " & ".join(cols + [ref]) + r" \\"


def body(L, R):
    syn = L == 3
    head = [r"$[[n,k,2]]_N$", "output", r"$\tau$", "stage"] + ([r"$n_\mathrm{syn}$"] if syn else []) + ["ref."]
    b = [r"\begin{tabular}[t]{l l c l " + ("c " if syn else "") + "l}", r"\toprule", " & ".join(head) + r" \\", r"\midrule"]
    prev = None
    for r in R:
        if prev is not None and r.get("block") != prev: b.append(r"\midrule")
        prev = r.get("block"); b.append(row(r, L, syn))
    return b + [r"\bottomrule", r"\end{tabular}"]


def tables(rows):
    """Tables SIV(a)-(c) as xtabular tables (package xtab), which break across columns and pages in two-column mode,
    so each table runs on from the last without white space."""
    SB4 = [r for r in json.load(open(os.path.join(O, SB_FILE)))["rows"] if r["l"] == 4]
    C4 = collections.Counter(r["k"] for r in SB4 if not r["third_level_target"]); T3N = sum(r["third_level_target"] for r in SB4)
    out = ["% generated by searches/table_sv.py from outputs/staged_best.json; do not edit by hand",
           r"\stepcounter{table}\edef\tabVInum{S\Roman{table}}\addtocounter{table}{-1}"]
    cap = {2: r"Smallest distance-2 factory found for each output class at $l=2$ (generators $\mathrm{S}$, $\CZ$), $k\le6$. "
              r"Output: the class, qubits labelled from $1$, linked to its circuit; $\tau$: its $Z(\pi/2^l)$-count, the least "
              r"number of $Z(\pi/2^l)$ rotations that prepare it (Corollary~\ref{cor:bound}; ``--'': not computed); stage: the "
              r"first stage reaching this $n$ (App.~\ref{app:cascade}): TG two-group, GS group-symmetric, GS-c with a group on "
              r"the checks, GL-d gluing on disjoint blocks, GL-e gluing on entangled support, SAT-r restricted SAT, SAT-f free "
              r"SAT; ref.: an earlier paper that gives this factory, or, where none does, ``syn'' if a synthillation construction~\cite{campbell2017unified} has the same $n$ and ``new'' if none does (at $l=2$, where every output is a stabilizer state, ``--''). $^\circ$: free SAT over all coefficients is unsatisfiable below this $n$ on $N$ qubits, so $n$ is smallest on $N$ or fewer qubits. Bold: on the bound $n=\tau+1$, so smallest at any $N$. Every class is in the "
              r"repository~\cite{singh_shraggymsd_simulations_2026}.",
           3: r"As Table~\ref{tab:smallest-l2}, at $l=3$ (generators $\mathrm{T}$, $\CS$, $\CCZ$): every class with $k\le4$, and "
              r"at $k=5,6$ the classes on the bound and, for each combination of gate types in the listed representative, the smallest $n$ and the largest saving over "
              r"synthillation, of the $766$ classes with $k\le6$. $[[14,1,2]]_4$ is also the $[[14,2,2]]$ of Ref.~\cite{bravyi2012magic} with one output kept, and its parent, Fig.~\ref{fig:parent}(a), is the $15$-rotation identity of Fig.~16(a) of Ref.~\cite{litinski2019magic} up to the signs of its odd-weight rotations, a Clifford change. $n_\mathrm{syn}$: fewest inputs of a synthillation "
              r"construction~\cite{campbell2017unified} (App.~\ref{app:compare}); $^\dagger$: $n<n_\mathrm{syn}$. Classes are taken up to CNOT, Pauli $X$ and diagonal Clifford gates on the outputs, so outputs related only through Hadamard gates, such as $\CS_{12,13}$ and $\CCZ_{123}$, have separate rows. Last block: "
              r"$\CCZ$ layers at $k=7,8$, beyond the class lists, from the cascade on a target list (End Matter, ``Beyond six "
              r"outputs''; Fig.~\ref{fig:new-factories}): the three smallest $n$ at each $k$, with $\tau$ exact at $k=7$ and "
              r"for $[[18,8,2]]_9$, and between proved lower and upper bounds for the other two rows at $k=8$; a factory on the bound with one check has the single-check form of case~11 of synthillation (Corollary~\ref{cor:saturation}), marked ``syn''.",
           4: r"As Table~\ref{tab:smallest-l2}, at $l=4$ (generators $\sqrt{\mathrm{T}}$, $\CT$, $\mathrm{CCS}$, $\mathrm{CCCZ}$; "
              r"their squares are $\mathrm{T}$, $\CS$, $\CCZ$): every class with $k\le2$, and the smallest and largest $n$ "
              + (r"among the $%s$ classes at $k=3$ and the $%s$ at $k=4$, of the $%s$ classes with $k\le4$ that need a fourth-level rotation. The other $%s$ classes of the $l=4$ list are third-level targets, covered by Table~\ref{tab:smallest-l3} with $T$ inputs (App.~\ref{app:target-classes})." % tuple(f"{x:,}".replace(",", "{,}") for x in (C4[3], C4[4], sum(C4.values()), T3N)) if not WITH_T3 else
               r"among the classes at $k=3$ and $k=4$, the $%s$ third-level targets included (App.~\ref{app:target-classes}); $\times2$: the $l=3$ factory of Table~\ref{tab:smallest-l3} for the same output, from the cited paper, with every coefficient doubled, so that each $T$ rotation is made from two $\sqrt{T}$ inputs." % T3N)}
    # Fixed column blocks (no automatic breaking): (a) and the first part of (b) share the left column of the first page,
    # (b) continues in the right column and in the left column of the next page, and (c) follows it there.
    def block(L, R, label, contd, env="table"):
        syn = L == 3
        head = [r"$[[n,k,2]]_N$", "output", r"$\tau$", "stage"] + ([r"$n_\mathrm{syn}$"] if syn else []) + ["ref."]
        spec = "l l c l " + ("c " if syn else "") + "l"
        b = [rf"\begin{{{env}}}[H]", rf"\renewcommand{{\thetable}}{{\tabVInum({'abc'[L - 2]}){' (cont.)' if contd else ''}}}"]
        if L > 2 or contd: b.append(r"\addtocounter{table}{-1}")
        ctext = (r"Table~\ref{tab:smallest-l%d}, continued." % L) if contd else cap[L]
        b += [r"\centering", r"\caption{\shnew{" + ctext + "}}"]
        if not contd: b.append(rf"\label{{tab:smallest-l{L}}}")
        b += [r"\scriptsize\setlength{\tabcolsep}{2.4pt}", rf"\begin{{tabular}}{{{spec}}}", r"\toprule",
              " & ".join(head) + r" \\", r"\midrule"]
        prev = None
        for r in R:
            if prev is not None and r.get("block") != prev: b.append(r"\midrule")
            prev = r.get("block"); b.append(row(r, L, syn))
        return b + [r"\bottomrule", r"\end{tabular}", rf"\end{{{env}}}", ""]
    R3 = [r for r in rows if r["l"] == 3]
    lines = lambda r: 1 + tex_output(r["spec"], 3).count(r"\\")
    cut, acc = [], 0                                  # split (b) into two blocks of about equal height
    total = sum(lines(r) for r in R3)
    for i, r in enumerate(R3):
        acc += lines(r)
        if acc >= total / 2: cut = i + 1; break
    # in-place ([H]) blocks, each column ended by hand: (a)+(b) | (b cont.) on one page, then (c) at the top of the next
    out += [r"\clearpage"] + block(2, [r for r in rows if r["l"] == 2], 2, False) + block(3, R3[:cut], 3, False)
    out += [r"\vfill\newpage"] + block(3, R3[cut:], 3, True)
    out += [r"\vfill\newpage"] + block(4, [r for r in rows if r["l"] == 4], 4, False)
    return "\n".join(out)


if __name__ == "__main__":
    main()
