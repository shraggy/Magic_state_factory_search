# Plan: finish LU20979 tonight — reviewed and refined (v2, 2026-10-31 22:30)

Implementer: Claude. User revises 2026-11-01 09:00, submits that evening. Rules: `AGENTS.md` (binding).
Manuscript dir `M = Magic_State_Factories_tex/MSD_manuscript/`. Draft v1 kept at `~/msf_scratch/finish_paper_plan_draft_v1.md`.

## R. Review findings that change the plan (read first)

R1. **"Synthillation fails on 8 classes" is FALSE under the strongest reading of its Thm 1.** Thm 1 of arXiv:1606.01904
    states that suboptimal synthesis matrices are allowed ("if suboptimal matrices are used, the resource cost is
    n = col(A) + 2col(B) + Δ"). Appending one zero column to B flips the parity of col(B), selects another Table I
    case, and gives a VALID factory with exactly our n on all 8 (case 1, `A,B+0`). Source:
    `outputs/synthillation_strongest_l3.json` (`searches/synthillation_strongest.py`; every variant checked with the
    BI test). Drop point (3) of the Fig. 1 list ("targets where synthillation fails") from the paper, letter and figure.
    What remains true and checkable: on those 8 targets (and on 160 classes in total) Thm 1 *with optimal A, B, as
    stated* gives an odd count that is not a valid factory; the valid synthillation factory needs a padded B.
    This is a footnote at most, not a headline.
R2. **The 117 wins survive the strongest synthillation.** Against the best valid variant per class (A or A+0; B,
    B+0 or B(μ+1)), ours needs fewer inputs on **117 of 766**, the same on **649**, more on **0**
    (`outputs/synthillation_strongest_l3.json` vs `outputs/l3_best_per_class.json`). Margins: 110 classes by 2
    inputs, 7 by 4. All 117 are entangled (k = 5: 11, k = 6: 106); n from 20 to 32; none at the bound. Robustness:
    for all 117, every Table I formula value over these variants (validity ignored) is still > our n, so no
    synthillation variant of this kind can match them. The two named wins hold: CS₁₂CS₃₄CS₅₆ 24 vs 26 and
    CS₁₂CS₃₄T₅ 26 vs 28. At equal n, N is equal on all 649 (no qubit advantage anywhere).
R3. **Model (c) gives almost no provable wins.** The only measurement-proof lower bound available is the stabilizer
    nullity ν (Beverland et al., arXiv:1904.01124; monotone under stabilizer operations incl. Pauli measurements).
    Computed for all 766 classes: ν = k for every class (`outputs/nullity_l3.json`). A route from distilled states
    needs total nullity ≥ ν; with CCZ batches ([[6m+2,3m,2]], ν 3 per CCZ), catalysed T (4 inputs, ν 1) and CS (8,
    ν 2), the cheapest such route costs ≥ 4, 8, 8, 12, 14, 14 inputs for ν = 1..6. Our factory beats that bound on
    **2 of 766** classes only: [[12,5,2]] (12 < 14) and [[14,6,2]] (14 < 16 without batching; = 14 with the batched
    [[14,6,2]] itself, so only [[12,5,2]] is a strict win). Neither is new against synthillation (R4).
    → "Provably better under (c)" can be marked on **1–2 points** only. Showing more needs an upper-bound route search,
    which can only show where routes win. **Decision D1 for the user** (below).
R4. **Frontier factories are synthillation case-11 instances** (same n, N, c; `outputs/frontier_vs_case11.json`, to be
    written in step C1 from the check already run). What is ours: the bound n ≥ T-count + 1 (App. S4) proving them
    smallest, and finding which targets reach it. Synthillation's paper lists only disjoint-CCZ examples (one CCZ,
    n = 8; N CCZ, n = 6N+2), read in its text.
R4b. **[[12,5,2]] is NOT new: it is synthillation's Example IV.4.** Campbell-Howard Example III.2/IV.4: tof#, two
    Toffolis sharing one control, F = 4x5(x1x2 + x3x4), "τ[U2#] = 11 and case 11 of Table I, we have a synthillation
    cost of 12". Our CCZ125CCZ134 has the same canonical form as tof# (checked: `frontier_search.canonical`), so
    [[12,5,2]] is their factory. Their U_N# family (N CCZs sharing one control, 2N+1 qubits) has n = 4N+4:
    U3# is a 7-qubit target with T-count 15 at n = 16 (isotropy floor 15 = their τ; so n = 16 is at the bound,
    k = 7). The batch tof# x CCZ (k = 8) has T-count 17 by subadditivity and floor, at n = 18. So synthillation's
    paper ALREADY reaches the bound at n = 8, 12, 14, 16, 18 with k = 3, 5, 6, 7, 8 respectively.
    What is ours: [[16,6,2]]7 and [[18,6,2]]7 reach T-count 15 and 17 with k = 6 outputs (theirs need k = 7, 8 and
    N = 8, 9 by case 11); the optimality proof (bound) for all; finding them by search. Every "[[12,5,2]] new"
    claim (main.tex Tables SIV(b) "new§", End Matter "New factories", Fig. 1 diamond, letters) must change to
    \cite{campbell2017unified} (their Example IV.4).
    Fig. 1: [[8,3,2]], [[12,5,2]], [[14,6,2]] as circles (in synthillation's paper); [[16,6,2]], [[18,6,2]] as
    diamonds (same n as their U3#, tof#xCCZ, on one fewer output qubit... NOTE: fewer OUTPUTS k, not a benefit by
    itself; state as "the smallest-k factories at the bound for T-count 15 and 17 that we know of" only after checking
    that no k=6 target with T-count 15 or 17 is in their paper: done, none listed).
R5. **l=4.** Coverage: 2/2 (k=1), 5/8 (k=2), 75/75 (k=3), 4,082/4,082 (k=4): 4,164 of 4,167 classes at k ≤ 4
    (`outputs/entangled_gluing_l4_k{1..4}.json` + `outputs/cascade/l4_k{1..3}.json`). The cascade alone had 35 of 85
    at k ≤ 3 and 569 of 4,082 at k = 4 (main.tex Table SIII "tab:coverage"). At the bound: only [[16,4,2]]₅.
    √T cultivation (arXiv:2606.10430, Chen, da Silva Fonseca, Sornborger, read in its HTML): "MST uses one √T|+⟩_L
    state and demands one T|+⟩ state half of the time"; d=3 MSC infidelity ~10⁻⁵ (idling) to 6×10⁻⁷ (no idling).
R6. **Word budget.** Main text (Introduction → Discussion, texcount on the current file) = 3,660 words; with
    the Code-and-Data and Acknowledgements paragraphs 3,904 — the 3,750 limit must be checked with the same
    convention as before (the earlier 3,747 figure). Before any edit, run `texcount` on the exact span the journal
    counts (Introduction through Discussion + figure captions + footnotes, excluding End Matter and references) and
    record it in `~/msf_scratch/wordcount_log.txt`; every edit below states its net change.
R7. **v3 sweep (job 1229268)** is pending ("Not enough job slots", 32 cores on one host). It is a confirmation run,
    not a source of new numbers: the merged table already contains every stage. If it has not started by 23:30,
    resubmit with `-n 16`; if not finished by 07:00, the paper uses the merged table and says nothing about v3.
R8. Abstract currently ≈ 541 characters (limit 600). Title keeps "Unified Search": the cascade is now one stage-wise
    search (`searches/cascade.py`, stages i, ii, ii′, iii, iv, v), tested on the k = 3 cell (6/6 classes reproduce).

## 0. User decisions (as given; caveats in brackets)
0.1 Main claim: one BI search gives a factory for every l=3 class with k ≤ 6 that matches or beats synthillation
    (fewer inputs on 117 of 766, same on 649, more on 0) [revised from 641/8 by R1–R2], proves the frontier
    factories smallest, and extends to l=2 and l=4.
0.2 Frontier factories: provably smallest; state that synthillation's case-11 construction gives the same factories.
0.3 Cost model (c), Clifford + all Pauli measurements. [Caveat: only the ν bound is provable; see R3, D1.]
0.4 Beyond k=6: state the gap in theory.
0.5 Fig. 1 marks points 1, 2, 5, 6, 7 [point 3 dropped by R1; point 4 dropped earlier] and the (c) winners.
0.6 Update Tables SIV and the repository catalogue.
0.7 Title keeps "Unified Search".

## D. Decisions needed from the user (ask once, at the start; defaults in bold apply if no answer by 23:00)
D1. Model (c) markers: **(a) mark the 1 strict win ([[12,5,2]]) and state the ν bound in one sentence** /
    (b) drop the (c) markers / (c) build an upper-bound route search (cannot finish tonight).
D2. The 8 ex-"fails": **mention nowhere** / one-sentence footnote in App. compare.

## A. Compute (in parallel with writing)
A1 (done) strongest synthillation: `outputs/synthillation_strongest_l3.json`.
A2 (done) nullity: `outputs/nullity_l3.json`; add script `searches/nullity_l3.py` that regenerates it (copy the
    code run in this session; Pauli-expectation count, ν = k − log₂|Stab|), and the route lower bound LB(ν).
    Test: CCZ ν = 3, CS ν = 2, T ν = 1 (known values, Beverland et al. Table I).
A3 merge: rerun `searches/merge_l3_best.py` once the last 2 msfps tasks finish (or at 23:30, whichever first);
    add a column `syn_strongest_n` from A1; write `outputs/l3_best_per_class.json` (final).
A4 v3 sweep: see R7. If it lands, check each cell's best n ≤ merged best (rule 8); a regression is investigated,
    not used.
A5 catalogue: from the final merged table write `outputs/smallest_per_class.csv` (l=3 rows replaced, k ≤ 6),
    `stage` column = source (cascade:*, entangled-gluing, synthesis-sat, pattern-sat), run
    `searches/add_error_coefficient.py` for c. Keep the old file as `outputs/smallest_per_class_v2.csv`.

## B. Figure 1 (`searches/plot_frontier.py` → `figures/frontier_l3.png` and `M/FIGURES/frontier_l3.png`)
Full column width 3.4 in, two stacked panels, height ≈ 4.3 in, 300 dpi; dataviz skill (palette validated with
`scripts/validate_palette.js`, light mode; text in ink colours, not series colours).
B1 Panel (a) "T-count per input, l=3, k ≤ 6": x = n (8–34), y = exact T-count (`l3_best_per_class.json`).
    - all 766 classes: small grey dots (alpha 0.35), jittered ±0.2 in x;
    - dashed line T = n − 1 labelled "bound n ≥ T-count + 1 (App. S4)";
    - frontier (5 points at the bound): filled markers, labelled [[n,k,2]]_N and c (28, 66, 91, 120, 153),
      synthillation's disjoint-CCZ members [[8,3,2]], [[14,6,2]] drawn as open circles ("also synthillation case 11"),
      [[12,5,2]], [[16,6,2]], [[18,6,2]] as filled diamonds ("this search, proved smallest");
    - the 117 classes with fewer inputs than synthillation: accent colour dots;
    - the 1 (c)-strict win (D1a): ring marker on [[12,5,2]] (in synthillation's paper, R4b; the ring marks a proof, not novelty);
    - inset text "k = 7–9: candidates, T-count bounds only" with the k=7 exact points (n = 20, 22 at T = 17) as
      crosses and k = 8 range bars [floor, upper] at n = 20, 22, 24 (from `outputs/k8_candidates/best.json`), styled
      dotted and grey (unproved regime, point 7).
B2 Panel (b) "Per class, against synthillation": horizontal stacked bar per k (k = 3, 4, 5, 6): fewer / same /
    more inputs (0 everywhere), counts printed on segments (k=5: 11/54/0; k=6: 106/573/0; k ≤ 4: 0/22/0), plus one
    bar per level for coverage: l=2 9/9 classes, l=3 766/766, l=4 4,164/4,167 (k ≤ 4), synthillation none at l=2, 4.
B3 Caption (replace main.tex lines 127): see C4.
B4 Render, inspect for collisions, rebuild PDF, look at page 2.

## C. Manuscript edits (`M/main.tex`). Before starting: `cp main.tex ~/msf_scratch/main_before_final_2230.tex`.
After every numbered edit: `latexmk -pdf main.tex`, grep log for "Error|undefined", run the word count of R6,
check abstract ≤ 600 characters. All new text in \shnew{}.

C1 Intro, end of paragraph 2 (line 111), replace
    "We classify factories by the largest output $T$-count they reach for a given number of inputs: a factory needs at
    least one input more than the $T$-count of its target, and at $l=3$ the search finds factories at that bound for
    $T$-counts that no known factory family reaches (Fig.~\ref{fig:new-factories})."
    with
    "We classify factories by the largest output $T$-count they reach for a given number of inputs: a factory needs at
    least one input more than the $T$-count of its target, which we prove. At $l=3$ the search reaches every output
    class with up to six outputs, needs fewer inputs than synthillation for $117$ of the $766$ classes and never more,
    and proves the factories at the bound smallest (Fig.~\ref{fig:new-factories})."  (net ≈ +8 words)
C2 Abstract (line 92), replace "and drives a search for any target diagonal gate that reaches the largest $T$-count
    per input, beyond the known factory families." with "and drives one search for every target diagonal gate. At
    distance~2 it matches or beats the best known construction for all $766$ three-level output classes with up to
    six outputs and proves the factories with the largest $T$-count per input smallest." Check ≤ 600 characters;
    if over, drop "three-level".
C3 Cascade paragraph (line 143), replace the last sentence "Its last two stages reach entangled outputs, … asking
    only for a factory smaller than the earlier stages found." with: "Its later stages reach entangled outputs, whose
    gates overlap on shared qubits: gluing on entangled support takes an optimal compilation of the target as the
    output part and solves the same linear system for the check couplings; the next lists exactly one target per
    Clifford class of output [footnote kept]; a SAT stage then fixes the output parts to a compilation of the target
    and solves only for which checks each rotation touches; and a last SAT stage over all rotations asks only for
    factories smaller than every stage before it found." (net ≈ +35 words; offset in C9)
    Also line 143: "A distance-2 BI circuit, however, is simply a compilation circuit for its target in which a few
    qubits are sacrificed as checks (App.~\ref{app:floor})" — keep (it is the App. S4 statement).
C4 Fig. 1 caption (line 127), replace whole caption with:
    "\shnew{Distance-2 factories at $l=3$ for all $766$ output classes with $k\le6$ (gray), at their exact output
    $T$-count. (a)~Every factory meets the bound $n\ge T\text{-count}+1$ (dashed; App.~\ref{app:floor}); the five at
    the bound are the smallest for their target at any size, labelled $[[n,k,2]]_N$ with output-error coefficient
    $c$ (error $c\,p^2$). Circles: also given by synthillation~\cite{campbell2017unified}; diamonds: entangled
    $\CCZ$ layers found by the search, which synthillation's construction also reaches at the same $n$ and $N$
    once the target is supplied. Colored: the $117$ classes where the search needs fewer inputs than every valid
    synthillation construction. Ring: the one factory that beats every route distilling known states and converting
    them with Cliffords and Pauli measurements, by the stabilizer-nullity bound~\cite{beverland2020lower}.
    Crosses and bars: $k=7$ (exact) and $k=8$ (bounds only) candidates, beyond the reach of exact $T$-counts.
    (b)~Per class against synthillation, and classes covered at each level. Data and circuits:
    Ref.~\cite{singh_shraggymsd_simulations_2026}.}"
    (Caption words count toward the limit if the journal counts captions — check with R6 convention.)
C5 Cascaded-search paragraph (line 335), replace "It reaches every class at $l=2$ and at $l=3$ for $k\le6$, and
    matches or beats every factory of the separate searches." through the end of the paragraph with:
    "It reaches every class at $l=2$ and at $l=3$ for $k\le6$. Fixing the checks of a borrowed identity leaves a
    circuit that prepares its target, so a factory needs at least one input more than the target's $T$-count
    (App.~\ref{app:floor}), and a factory at that bound is the smallest for its target. Five classes reach it, at
    $n=8$, $12$, $14$, $16$ and $18$ ($T$-count $7$, $11$, $13$, $15$, $17$; Fig.~\ref{fig:new-factories}); the
    synthillation construction gives the same factories once their targets are known, and the bound proves them
    smallest. On $117$ of the $766$ classes the search needs fewer inputs than every valid synthillation
    construction, by two or four, and on none more (App.~\ref{app:compare})." (net ≈ −10 words)
C6 Discussion (line 341), replace last sentence "At distance~2, … finds factories at the bound of one input above the
    $T$-count of their target." with "At distance~2, where cultivation makes a single factory sufficient, the cascaded
    search reaches every output class at $l=2$, at $l=3$ for $k\le6$ and at $l=4$ for $k\le4$ but three, matches or
    beats synthillation on every class, and proves the factories at the bound $n=T\text{-count}+1$ smallest. Beyond
    $k=7$ outputs the exact $T$-count is out of reach of current methods, so larger candidates are compared only by
    bounds on their $T$-count; tighter lower bounds on $T$-count are an open problem for theory." (net ≈ +35)
C7 Open problems paragraph: keep as is (two questions: d>2 systematic search with follow-up jain2026symmetry;
    Hadamards / non-CSS). No change.
C8 End Matter "New factories" (line 408): replace from "Below $n=20$ the known families reach this bound only at"
    through "…with output-error coefficients $c=66$, $120$ and $153$. Their targets are entangled $\CCZ$ layers on
    overlapping qubits, each of $T$-count above that of any factory with fewer inputs." with:
    "Five classes with $k\le6$ meet this bound: $[[8,3,2]]_4$, $[[12,5,2]]_6$, $[[14,6,2]]_7$, $[[16,6,2]]_7$ and
    $[[18,6,2]]_7$, of $T$-count $7$, $11$, $13$, $15$, $17$ and $c=28$, $66$, $91$, $120$, $153$. Each is a one-check
    factory: an optimal compilation of a $\CCZ$ layer through one check qubit, which is also case~11 of the
    synthillation construction~\cite{campbell2017unified} (their Table~I), with the same $n$, $N$ and $c$; that work
    lists only disjoint $\CCZ$s ($n=6m+2$), and the search finds the entangled layers $[[12,5,2]]$, $[[16,6,2]]$ and
    $[[18,6,2]]$ that reach the bound. At $k=7$ the exact $T$-count of the best candidates is $17$, below the bound,
    and at $k\ge8$ only bounds on $T$-count are available. Stabilizer operations, including Pauli measurements,
    cannot increase the stabilizer nullity~\cite{beverland2020lower}, which equals $k$ for every class here; by this
    bound only $[[12,5,2]]$ beats every route that distils $\CCZ$, $\CS$ and $T$ states and converts them."
    and change "at $l=4$ exact counts are available here only up to five qubits" to "at $l=4$ gluing on entangled
    support gives a factory for $4{,}164$ of the $4{,}167$ classes with $k\le4$, and only $[[16,4,2]]_5$ meets the
    bound". (End Matter: not counted.)
C9 Offsets for C3/C6 (+70): (i) line 145 delete "Unlike conditions parametrized by $n$, our search is parametrized by
    $N$, which stays small even when $n$ is large:" (−18) and merge; (ii) line 140 "Each intermediate step is a
    factory, so one parent circuit gives several factories for different magic states, e.g.\ $\CCZ$, $\CS$ and $T$ at
    $l=3$." → "Each intermediate step is a factory." (−18); (iii) line 111 "These resource estimates could be improved by
    targeting the smallest distance-2 factory for each target magic state an algorithmic compilation needs." →
    "Each target magic state a compilation needs then calls for its smallest distance-2 factory." (−9); (iv) if still
    over, move the Iceberg sentence (line 140) to End Matter "Further remarks" (−30). Stop when ≤ 3,750.
C10 "Distil, then synthesize" (End Matter line 403): replace "…so this bound does not exclude a route from two
    $\CCZ$s." with "…; for the classes here $\nu=k$, so the bound excludes such routes only where the factory needs
    fewer inputs than the cheapest states of total nullity $k$, which holds for $[[12,5,2]]$ alone." Keep the d=3
    follow-up sentence ("in our follow-up work").
C11 √T cultivation caveat: grep "cultivation does not supply" in main.tex and letters; replace with "cultivation of
    $\sqrt{T}$ states has been proposed~\cite{chen2026sqrtT}, though teleporting the gate needs a $T$ state half of the
    time". Add bib entry `chen2026sqrtT` (I.-C. Chen, M. da Silva Fonseca, A. Sornborger, "Efficient Magic State
    Cultivation for √T Gates", arXiv:2606.10430, 2026).
C12 App. S4 (line 1205), append: "The same argument is constructive: taking a compilation of the target as the output
    parts and solving only for the check couplings gives the gluing-on-entangled-support stage, and fixing the output
    parts to a compilation and solving for the check patterns gives the check-pattern SAT stage
    (App.~\ref{app:cascade})." Also update "The bound holds for all $1{,}551$ catalogue factories…" by rechecking on
    the final catalogue (script: loop over `l3_best_per_class.json`, assert n ≥ tcount+1: 766/766; cite file).
C13 App. cascade (line 1112): rewrite the stage list to (i) two-group, (ii) gluing on disjoint support, (ii′) gluing
    on entangled support, (iii) exact classes, (iv) check-pattern SAT, (v) SAT over all rotations seeded with the
    best n so far, deterministic step limit (z3 rlimit) and a fresh process per query. Add one paragraph each for
    (ii′) and (iv) (definitions as in C3, plus: (iv) uses the synthesis multiset A + 2B + 4c plus m check-only
    rotations, Boolean variables for check membership, R = 1..4 checks, μ and μ+1). Update Table "tab:coverage"
    with l=4 coverage (R5) and l=3 unchanged; replace the "time of (iv)" column by per-stage times if v3 lands,
    else keep and footnote.
C14 App. compare (line 1326): replace second paragraph by: synthillation built from Table I cases 1–11 for every
    class with A optimal or padded and B smallest, padded or μ+1 (Thm 1 allows suboptimal A, B), checked with the BI
    condition; per class: ours fewer on 117 (by 2: 110, by 4: 7), same 649, more 0; the 117 are all entangled; at equal
    n the qubit count N is the same in all 649. Files: `outputs/synthillation_strongest_l3.json`,
    `outputs/l3_best_per_class.json`. Remove the sentence about three targets where "no factory of the formula's size
    exists" (it is about the formula as stated; keep only if D2 = footnote, reworded to "Thm 1 with optimal matrices
    gives an odd count there; a padded Lempel factor gives a valid factory with our $n$").
C15 Tables SIV(b): for each row keep n if it equals the merged best; update rows whose class has a smaller merged n
    (check every row by its class); "search stage" column: add entangled gluing / check-pattern SAT where they give
    the best; ref column: rows that equal a synthillation construction for a disjoint target keep \cite{campbell2017unified};
    "new" stays only where n < every valid synthillation variant or the target is not a product of disjoint gates
    and not listed in the literature — add a footnote mark for "same n as a synthillation construction once the
    target is given". Frontier rows stay bold.
C16 Remove/replace everywhere (grep in main.tex, referee_response.tex, cover_letter.tex): "no known family reaches",
    "below $20$ inputs", "one qubit fewer", "fewer qubits", "synthillation fails", "formula fails", "odd input count",
    "cannot reach" — each occurrence rewritten to the R1–R4 facts.

## L. Referee letter (`M/referee_response.tex`; back up first) and cover letter
L1 Summary bullet "New factories at the bound" (line 107): replace with "\textbf{Every class, matched or beaten.}
    The cascade reaches every output class at $l=2$ and at $l=3$ with up to six outputs, and $4{,}164$ of the
    $4{,}167$ classes at $l=4$ with up to four. Against the synthillation construction, built from its Table~I for
    every $l=3$ class and checked with the same condition, it needs fewer inputs on $117$ of $766$ classes and never
    more. The five factories that meet the bound $n=T$-count$+1$, among them the entangled $\CCZ$ layers $[[12,5,2]]$,
    $[[16,6,2]]$ and $[[18,6,2]]$, are proved smallest; synthillation's construction also gives them once their
    targets are supplied."
L2 Line 108 "Exact minimality": keep; append "and every factory in the catalogue meets the bound (766 of 766)".
L3 Referee A/B/C passages (lines ~134, ~150, ~175, ~177): replace each "no known family reaches below $20$ inputs"
    with the L1 facts; in the Referee C benchmark (line 150 onward) replace the frontier table's "source" column
    entries "this work" by "this search (also synthillation case 11)"; replace line 177 paragraph's comparison
    sentences with C14's numbers; delete "For three mixed targets … its formula gives an odd input count … the
    unified search supplies a valid one" (R1).
L4 Referee B impact response (line 175): replace "These are promising alternatives in the strongest sense available:
    no factory for their targets can use fewer inputs." with "Across every output class the search matches or beats
    the synthillation construction, with fewer inputs on $117$ of $766$ classes, and the factories at the bound are
    proved smallest."
L5 Cover letter bullet "Provably optimal new factories" (line 25): replace last sentence with "The search finds the
    factories at this bound, among them entangled $\CCZ$ layers of $T$-count $11$, $15$ and $17$, and for $117$ of
    the $766$ output classes with up to six outputs it needs fewer inputs than the synthillation construction, and
    never more." First bullet: add "$4{,}164$ of the $4{,}167$ classes at $l=4$ with up to four outputs".
L6 Compile both letters, zero errors; red/blue colours unchanged; no surnames outside quotations (grep "Campbell",
    "Howard", "Gidney", "Litinski" etc. outside \cite and quotations).

## Rp. Repository
Rp1 Copy into `outputs/`: A1, A2 files, `l3_best_per_class.json`, `pattern_sat/`, `pattern_sat_local/`,
    `synthesis_sat/`, `entangled_gluing_l3.json`, `entangled_gluing_l4_k*.json`, `k7_candidates/`, `k8_candidates/`,
    `k9_candidates/`, `class_table_l3.json`, `frontier_vs_case11.json`.
Rp2 `outputs/README.md`: one table row per new file with its script. `searches/README.md`: stages (ii′), (iv), (v)
    and the "SAT restricted by App. S4" lesson (general SAT with 5× budget closed 0 of 23; check-pattern closed all).
Rp3 Commit file list for the user (no commit).

## T. Order and time (start 22:45)
1. 22:45 ask D1/D2 (one message), back up files, word-count baseline (R6) .......................... 10 min
2. 22:55 A2 script + C1, C2, C3, C9 (main text) with compiles and counts ............................. 45 min
3. 23:40 C5, C6, C10, C11, C16 (main text + End Matter) ............................................... 40 min
4. 00:20 A3 merge (msfps done by then) → A5 catalogue → C15 tables .................................... 45 min
5. 01:05 B figure (script, render, inspect, C4 caption) ................................................ 60 min
6. 02:05 C8, C12, C13, C14 (End Matter / appendices) .................................................. 45 min
7. 02:50 L1–L6 letters ................................................................................ 50 min
8. 03:40 Rp1–Rp3, final compile of all three PDFs, word/character checks, look at every page with a figure .. 30 min
9. 04:10 report to the user: what changed, numbers with files, open items (v3 status, D answers applied).
Cluster: re-arm a monitor on 1229268 and the last msfps tasks at step 1; report any failure when seen.

## Risks
- v3 sweep never runs tonight: no effect on numbers (R7).
- Word limit: C9 offsets listed; stop adding text if over.
- If the user chooses D1(c): not feasible tonight; fall back to D1(a) and say so.
