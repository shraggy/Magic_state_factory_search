# Context for the plan review: everything learnt, tried and decided (2026-10-31)

Read with `plans/finish_paper_plan.md`. Every number below has its source file; "unchecked" marks anything not
backed by a rerunnable check. Repository root: `Magic_state_factory_search/`. Rules: `AGENTS.md` (read it).

## A. The paper and its situation
- PRL manuscript LU20979, "Borrowed Identities: Malleability and Unified Search of Distillation Circuits Across
  the Clifford Hierarchy", second-round revision. Files: `Magic_State_Factories_tex/MSD_manuscript/main.tex`,
  `referee_response.tex`, `cover_letter.tex`, `distill.bib`. Main text <= 3,750 words (End Matter excluded),
  abstract <= 600 characters. Main text currently ~3,747 words. User revises 2026-11-01 09:00, submits that evening.
- Referees (first round): excited about the theory (borrowed-identity condition), not about the search or
  malleability. Referee B: the framework did not find promising alternatives to known schemes ([[11,1,2]] worse
  than [[10,2,2]]); limited at d>2; speed not relevant. Referee C asked for a systematic benchmark.
- Style rules: new text in \shnew{}; plain compact English; letter comments red, responses blue; no author
  surnames outside quotations; thank referees; follow-up works (arXiv:2609.30860, 2610.06535, jain2026symmetry)
  are called "follow-up" (this work came first); no failed attempts described in the paper; no symbol in the main
  text without definition (no tau, n_S in main text except the synthillation footnote); notation [[n,k,d]]_N,
  \mathcal{C}_l for hierarchy levels, \mathbb{C} for circuits.

## B. Theory facts (proved or checked)
- Borrowed identity (BI): a distance-2 diagonal factory is a circuit of Z(pi/2^l) parity rotations on outputs +
  checks that is the identity on the checks (up to Clifford) and deposits the target on the outputs.
- App. S4 (`app:floor`, main.tex ~line 1205): fixing the checks to a basis value leaves a circuit preparing the
  target, so n >= tau(U) + n_S; a CNOT-shift lemma gives n_S >= 1, so n >= T-count + 1 for every factory.
  This is a BI theorem (uses diagonal rotations and |+> checks). Synthillation states only a cost, no lower bound.
- One-check condition at level l (verified 4,000/4,000 at each l): overlaps of |A| outputs vanish mod
  2^(e(|A|+1)-|A|). l=3: self-orthogonal; output is a CCZ layer. Optimal compilation + one check (+ check-only
  rotation if tau odd) gives a factory at n = tau+1 for every pure-CCZ target (109 random layers verified).
- Isotropy floor (Mazumder, arXiv:2610.01024): T-count >= 2(k-d)+1 for pure cubic; m disjoint CCZ have T-count
  exactly 6m+1. Stabilizer nullity (Beverland et al. arXiv:1904.01124).
- Exact T-count: Amy-Mosca RM decoder up to k=6 (`classification/tcount.py`); k=7 via `searches/tcount7.py`
  `exists_below2` (two information sets; agrees with exact decoder 12/12 at k=6; seconds at k=7). k=8: coset
  dimension 163 on 255 coordinates, infeasible (timed out at 280 s, combinatorial).

## C. Synthillation (Campbell-Howard arXiv:1606.01904) - read from the text
- Thm 1: n = tau[U] + 2 mu[U] + Delta, 0 <= Delta <= 11; constructive via Table I (11 cases of matrix G built from
  A (synthesis), B (Lempel factor), c (linear terms)). CCZ circuits: case 9 (col(A) even, G=[A;1], Delta=0) or
  case 11 (col(A) odd, G=[A 0;1 1], Delta=1): ONE check row. Explicit examples: one CCZ n=8 (Eq. 114); N CCZ n=6N+2.
  No entangled CCZ layers listed; no T-count 11/15/17 examples (searched text). Claims "reduced resource cost";
  no optimality claim for n.
- Our implementation: `searches/synthillation_compare.py` (cases 1-8 match Table I; cases 9-11 initially used a
  wrong 2-row template - FIXED by rebuilding case 9/11 directly). `searches/class_table_l3.py` builds it for every
  class (mu and mu+1) and checks with our BI test: valid on 758 of 766; 8 fail (CCZ layers with one T).
- With true case 11 the frontier factories are IDENTICAL to synthillation's (same n, N, c): CCZ (8,4,28),
  [[12,5,2]] (12,6,66), [[14,6,2]] (14,7,91), [[16,6,2]] (16,7,120), [[18,6,2]] (18,7,153). Our earlier claim
  "one qubit fewer" was a script artefact and is withdrawn.

## D. The search (cascade) and today's changes
- `searches/cascade.py` stages now: (i) two-group, (ii) gluing (disjoint symmetric blocks; = synthillation family
  for product targets), (ii') gluing on entangled support (NEW: output part = optimal network of the target, check
  couplings from the same 2-adic solve), (iii) exact target classes, (iv) check-pattern SAT (NEW, l=3: output
  parts fixed to a synthesis multiset A + 2B + 4c + m check-only rotations; SAT over Boolean check patterns only),
  (v) quick SAT over all rotations, seeded with best n so far, deterministic z3 rlimit. Old cascade backed up at
  `~/msf_scratch/cascade_before_stages.py`. Test: k=3 cell reproduces stored best n on 6/6 classes.
- Clifford freedom in gluing generalised to every level (cliff_exp); l=3 regression tests pass (CCZ n=8,
  T1T2CCZ125CCZ134 n=22).
- Main sweep v3: job 1229268 on ccc4, 18 cells (l=2,3,4 x k=1..6), 32 cores/cell, rlimit 6e8 -> 
  `outputs/cascade_main_v3/` (running at time of writing).

## E. Results (with files)
- Per-class l=3, k<=6, 766 classes (`outputs/l3_best_per_class.json`, `searches/merge_l3_best.py`): vs
  synthillation fewer 117, equal 641, more 0, synthillation fails 8. By k: k<=4 0/21/0/1; k=5 11/52/0/2; k=6
  106/568/0/5. Best-factory source: entangled gluing 336, synthesis-SAT 195, cascade 147, check-pattern 88.
  Every factory verified with is_factory, distance, and output deposit (882 SAT factories rechecked 882/882).
- History of the comparison (do NOT put in paper): initial cascade lost 614 of 758 to synthillation (short SAT
  budget, symmetric-only gluing); entangled gluing -> 282 losses; synthesis-restricted SAT (`synthesis_sat.py`,
  job 1217503) -> 79; check-pattern SAT (`pattern_sat.py`, msfps arrays) -> 0. Larger rlimit general SAT (job
  1218993) closed none of 23 tried: the restriction from App. S4 is what works.
- Last class: T5T6CS14CS23CCZ135CCZ156CCZ245, n=34 N=9 (= synthillation), `outputs/pattern_sat_local/657.json`;
  UNSAT for n=29..33 holds only for that synthesis multiset with 3 checks.
- Frontier at l=3 (`outputs/frontier_l3.json`): n = 8, 12, 14, 16, 18 at T-count 7, 11, 13, 15, 17; at bound.
  Class table: exactly 5 classes at the bound among 766.
- k=7 (`outputs/k7_candidates/`): candidates [[20,7,2]], [[22,7,2]] have exact T-count 17 (not 19/21): no frontier
  point. k=8/9 (`outputs/k8_candidates/best.json`, `k9_candidates/best.json`): 4,463 / ~4,500 candidates; proved
  floor <= 17 (k=8), <= 19 (k=9) never reaches n-1 at new n; heuristic upper bound = n-1 for many (e.g. 182 of 240
  at n=20, k=8) - unproved. This is the "gap in theory".
- l=4 (`outputs/entangled_gluing_l4_k{1..4}.json`): entangled gluing builds a factory for 75/75 classes at k=3 (47
  new vs cascade), 4,082/4,082 at k=4 (cascade had 0). Coverage k<=4: 4,164 of 4,167 (cascade alone had 35 of 85 at
  k<=3); 3 uncovered classes at k=2 (sT0 CT01, sT0 sT1 CT01, sT0 sT0 sT1). Only [[16,4,2]]5 (CCCZ, sqrt(T)-count 15)
  at the bound; median gap ~24. Old l=4 greedy counts (job 1206194) matched exact 16/16 at k=5 but are upper bounds.
  Cultivation of sqrt(T): arXiv:2606.10430 (T correction about half the time) - caveat.
- l=2: Iceberg [[k+2,k,2]] closed form (one-check two-group), certified Z(pi/4)-count k+1 at k=2,4,6; cite Steane
  PhysRevA.54.4741 and Gidney et al. yoked surface codes.
- Distill-then-synthesize (unitary, `outputs/min_synthesis_l3_k{3,4,5}_*.json`, k=6 missing: orbit table >
  MAX_BITS): 84 classes k=3-5. (a) catalysed T (4/T): direct fewer 30, equal 17, route fewer 37. (b) MEK (5/T):
  48/8/28. CS/CCZ-only targets: direct 13, equal 2, route 2 of 17. Route wins are targets with T parts. The user's
  conjecture "CCZ-type factories capture all lower-T-count targets" does NOT hold as stated; user chose model (c)
  (Clifford + all measurements) as the correct competitor; not yet computed.

## F. Things that are false or retracted (must not appear)
- "No known family reaches these T-counts below 20 inputs" (synthillation formula does).
- "One qubit fewer than synthillation" for frontier factories (false; same N).
- Any "x% fewer" cost claims, App. S6 cost comparison (removed), d>2 malleability conjecture (removed).
- "The cascade matches or beats every factory" without the merged table and counts.
- l=4 "counts" from greedy as exact.

## G. Compute rules and state
- Laptop: 16 cores, leave 5 to system; workers <= 16 - 5 - load. Cluster ccc4 (LSF): workers = cores requested,
  one host, deterministic z3 rlimit, fresh process per query; pending-job threshold forces chunked submission.
- Running now: 1229268 (v3 sweep), last 3 msfps tasks. Monitor task b2hoe23t4.

## H. User's open requests (keep)
- Finish paper, letter, cover letter, Fig. 1 tonight; Fig. 1 must carry points 1,2,3,5,6,7 and mark factories
  provably better and new under model (c).
- Model (c) competitor computed for the per-class comparison.
- Repository: copy all results into outputs/, update READMEs, prepare commit list (user commits).
