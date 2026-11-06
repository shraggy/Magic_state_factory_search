# Rewrite plan for LU20979 (draft for the user's approval; nothing has been edited yet)

Date: 2026-11-01. Backup of the current state: `~/msf_scratch/refine_0905/`.

## A. Your inputs since 08:30, and how the plan answers each

| # | Your input | How the plan handles it |
|---|---|---|
| 1 | Merge iv and v? Do they differ in runtime? | One "SAT" step with two encodings: restricted (App. S4, fast) then unrestricted (catches what lies outside a synthesis). Data: restricted/gluing strictly better on 619 classes, unrestricted on 16, all 16 of which beat synthillation. Section D, stage 4. |
| 2 | How much does two-group help? Merge paragraphs (c) two-group and (d) cascade | Two-group gives the best factory for 9 of 766 classes at k <= 6, but 1,242 of 1,291 factories in the large sweep (reach, closed forms). (c) and (d) become one paragraph built on one condition. Section D. |
| 3 | Two-group equations belong in the appendix; [[12,2,2]] is from the symmetric case; only [[56,6,2]] is from two-group | Lemma 4, Theorem 6 and Eq. (N-asym) move to App. S3 with their labels. [[12,2,2]] is credited to the symmetric/sequential construction, [[56,6,2]] to two-group. |
| 4 | Describe all stages with the interesting [[n,k,d]] each reproduces, total vs best vs new counts | One table in the main text (Table I): per stage, what is fixed, what is solved for, classes reached, best found, wins over synthillation, and examples. Numbers from `outputs/stage_contributions_l3.json`. |
| 5 | Keep the new condition's notation close to the previous conditions | Same symbols as Theorem 3: a sum of rotation counts over supersets, in $2^{l-|T|+1}\mathbb{Z}$, for $1\le|T|\le l$. Theorem 3 is its $S_N$ case and Theorem 6 its $S_k\times S_{N-k}$ case. Exhaustively checked. |
| 6 | End Matter looks like 2.5 pages | Cut to 2 pages or fewer (PRL limit, counted in pages). Table E1 and malleability remarks go to the Supplement. Measured on the PDF. |
| 7 | Refine thoroughly: whole main.tex, letter, cover letter; one message; free reframing; letter answers first, then points to where moved content now lives; cover letter on the same premise; do this every time | Sections B–G. Saved as a standing rule in memory. |
| 8 | Choose the title to fit the premise | Section B. |
| 9 | Full freedom (you have a backup); a crisp PRL on the borrowed-identity idea | The whole main text is restructured around the BI condition, not around the search's history. |
| 10 | Campbell–Howard PRL 118, 060501 as inspiration (style only, no content or equations); two single-column figures, several equations and matrices, about 4.5 pages | Same shape: one framing idea, the main condition early, constructions inline, two single-column figures, compact equations. None of their equations or text. |

## B. Single message and title

**Message.** A distillation circuit need only return one state, $\ket{+}^{\otimes N}$, to itself. This borrowed identity is a single divisibility condition, the same at every level of the Clifford hierarchy. Each known construction and each search stage is this condition with different parts held fixed. Fixing the checks proves that every factory needs at least one input more than its target's $T$-count.

**Title (proposed).** "Borrowed Identities: One Condition for Magic-State Distillation Across the Clifford Hierarchy".
Alternatives:
- "Borrowed Identities for Magic-State Distillation";
- keep the current "…Malleability and Unified Search…" (the letter already explains that title change).

**Abstract (draft, ≈590 characters, to be measured).**
"A magic-state distillation circuit need only act as the identity on a single input state, a *borrowed identity*. We show this is one divisibility condition, unchanged across the Clifford hierarchy. It contains the Reed–Muller, hypercube and synthillation factories as special cases, gives one parent circuit several factories (malleability), and proves that a factory needs one input more than its target's $T$-count. One search over this condition reaches every distance-2 factory class we can list, and matches or beats synthillation on all 766 three-qubit-level output classes with up to six outputs."
(Wording to be tightened to ≤ 600 characters. "three-qubit-level" will be replaced by a correct short phrase, e.g. "$T$-type".)

## C. Main text (target ≈ 3,650 words; two single-column figures)

1. **Introduction** (≈ 650 words; existing paragraphs compressed)
   - Cultivation makes a single distance-2 factory enough; the figure of merit is $n$ per target (keep the $np$ argument).
   - Heisenberg picture versus the BI (Litinski observation).
   - One paragraph stating the results: the condition, closed forms, malleability, the bound, the search, and the synthillation comparison. No history of the search.
2. **Framework** (≈ 450 words)
   - Gates (Eq. 1) and Definition 1, as now.
   - **Theorem 1 (BI condition), new main-text statement:**
     $$\sum_{S\supseteq T} c_S \in 2^{\,l-|T|+1}\mathbb{Z}\quad\text{for all } T,\ 1\le|T|\le l.$$
     Verified by `searches/bi_condition_check.py`: 2,121,136 of 2,121,136 circuits agree with brute force.
   - Factory extraction in one sentence, and the target $U$ as the removed rotations.
   - **Corollary (bound):** fix the checks, and the CNOT-shift argument gives $n\ge\tau(U)+1$. Five-line proof, moved from App. S4 into the main text because it is a headline.
3. **Symmetric solutions** (≈ 500 words)
   - Lemma 2 and Theorem 3 kept, labelled as the $S_N$ case of Theorem 1.
   - Reed–Muller, hypercube and Iceberg as one paragraph (Iceberg comes back from the End Matter).
4. **Malleability** (≈ 450 words)
   - Proposition, four-qubit parent: [[8,3,2]] → [[12,2,2]] → [[14,1,2]], plus [[10,2,2]] and [[11,1,2]]; one sentence on the five-qubit and $l=4$ chains.
5. **One search, five ways to fix the condition** (≈ 700 words + Table I)
   - Table I, from `outputs/stage_contributions_l3.json`:

     | stage | held fixed | solved for | reached | best (first) | wins over synthillation | examples |
     |---|---|---|---|---|---|---|
     | two-group | $S_k\times S_{N-k}$ | signs | 13 | 9 | 0 | [[8,3,2]], [[56,6,2]], [[20,4,2]] |
     | gluing, disjoint | blocks | check couplings | 18 | 15 | 2 | [[14,6,2]], [[18,5,2]] |
     | gluing, entangled | compilation | check couplings | 766 | 425 | 76 | every class |
     | SAT | compilation (restricted) / nothing (unrestricted) | check pattern / all rotations | 675 / 742 | 301 / 16 | 23 / 16 | [[16,6,2]], [[18,6,2]] |

   - Totals: 766 classes, 117 wins, 5 at the bound. Also: large sweep 1,291 factories at $l=3$ ($N\le18$, $k\le12$), and $l=4$ 4,167 of 4,167 classes at $k\le4$.
   - Fig. 1 (single column, existing).
6. **Results against prior work** (≈ 400 words)
   - Synthillation built from its Table I, strongest variants; 117 / 649 / 0.
   - The frontier, stated honestly (case 11; Example IV.4 = [[12,5,2]]).
   - The nullity sentence.
7. **Discussion** (≈ 350 words)
   - Message restated; the $T$-count theory gap beyond $k=7$; the two open problems (d > 2 with follow-up jain2026symmetry; Hadamards/non-CSS).

**Second single-column figure.** Fig. 2 = the four-qubit parent circuit, shrunk to one column (existing `FIGURES/qk_four.pdf`, panel a). It shows malleability visually. Alternative: keep only Fig. 1 if the page budget is tight.

## D. What moves where

- **To App. S3 (two-group):** setup paragraph, Lemma 4 (`lem:dicke-asym`), Theorem 6 (`thm:asymmetric-master`), Eqs. `phase-asym`, `asymmetric-master`, `N-asym`, with labels unchanged so all 30 cross-references resolve.
- **To the Supplement:** Table E1 (stage timings), the malleability remarks, the Clifford-padded discussion.
- **Into the main text:** the bound proof from App. S4 (App. S4 keeps the 766/766 check and the constructive remark).
- **Cascade appendix:** replace the printed monomial condition, which disagrees with the code on 54 of 240 test circuits, with Theorem 1 as stated in the main text. The stages are described as five ways to fix Theorem 1.
- **Removed:** duplicated stage descriptions (main-text paragraphs (c) and (d), the End Matter "Unified distance-2 search").

## E. End Matter (≤ 2 pages, measured on the PDF)

1. **Scope and completeness** (≈ 250 words)
2. **New factories and the frontier** (≈ 350 words, from the current "New factories")
3. **Distil, then synthesize; cultivation at l=4** (≈ 150 words: nullity, √T caveat)
4. **Fig. E1** (literature recovered), kept if it fits; else moved to the Supplement and cited.

## F. Referee letter (rewrite in place; red comments and blue responses kept; no surnames)

- **Summary of main changes:** rewritten to the new premise and title; stages as in Table I; coverage numbers current.
- **For each comment about content that has moved:**
  - Lemma 4, Theorem 5/6, Eq. (7), (0,0) sector — Referee A comment 2, Referee B comment 5, Referee C comments 3 and 4;
  - stage descriptions — Referee A comment 4;
  - Table E1 — Referee A comment 6, Referee B speed paragraph.

  The pattern for each is: (i) answer the question as now; (ii) "In the resubmitted draft, [statement] is in App. S3 (Lemma 4, Theorem 6); the main text now states the general condition, Theorem 1, of which it is the $S_k\times S_{N-k}$ case."
- **Referee A comment 4:** stage list replaced by the Table I description; remove "58 classes in 45 minutes".
- **Referee B minor comment 1:** describe the new Figs. 1 and 2.
- **Page and section references** ("p. 3", "p. 4, after Theorem 5", ...) updated after the final compile.

## G. Cover letter

- Same premise and title; the "What is new" bullets become:
  1. one condition;
  2. the bound;
  3. one search and its result against synthillation;
  4. corrections.
- Restructuring paragraph: say that the main text now leads with the condition and the specific constructions are its special cases.
- New paragraph "Why this is a Letter" (section J): the qualities that make a PRL Letter, stated for this paper.

## H. Checks after each step

- Compile all three: zero errors, no undefined refs or citations.
- Main text ≤ 3,750 words (texcount, abstract through Discussion), abstract ≤ 600 characters, End Matter ≤ 2 pages (PDF page span).
- Every number quoted with a file in `outputs/` (rule 2). No new equation unless checked (Theorem 1: `bi_condition_check.py`; bound: App. S4 proof).
- Backups before each document.


## J. What the reference PRL shows a PRL Letter does (Campbell and Howard, PRL 118, 060501; read in full, style only)

Read from its text (`/tmp/chprl/p.txt`, arXiv:1606.01906); nothing of its content or equations is reused.
1. **One named paradigm and one reframing.** The abstract names the standard approach in a phrase
   ("distill-then-synthesize"), then states a single unifying alternative. → Ours: the Heisenberg-picture
   condition versus the borrowed identity, a one-state condition.
2. **Main result in the abstract and on page 1.** The headline cost theorem appears on p. 2, before any
   construction; the rest of the Letter proves it. → Ours: Theorem 1 (BI condition) and its corollary (the bound)
   stated early in the Framework, everything after as its special cases.
3. **One concrete, practically relevant example, carried by a figure.** Their Fig. 2 shows one gate (tof#) with a
   named saving. → Ours: Fig. 1 with every class against the bound and against synthillation, and one named
   example in the text (e.g. CS12CS34CS56 at 24 vs 26; the [[16,6,2]] and [[18,6,2]] frontier).
4. **Proof sketched inline for one case; the rest deferred.** They prove the simplest case in the Letter and say
   the others "follow the same methodology" (companion paper). → Ours: prove the bound inline (five lines); the
   symmetric and two-group cases as special cases; proofs in the Supplement.
5. **Self-contained definitions, short.** Notation is introduced once and reused (their D3, tau, mu). → Ours:
   theta = pi/2^l, D_l, n_T, tau(U); nothing else in the main text.
6. **Practical payoff stated plainly, no overclaiming.** They state savings as a ratio with its conditions. →
   Ours: "never more inputs than synthillation; fewer on 117 of 766 classes, by 2 or 4", with the conditions.
7. **Length.** About 3,400 words of text, three figures, 5 pages with references. → Ours: main text at the
   3,750 limit, two single-column figures.

**Use in the cover letter ("PRL type of additions").** A paragraph "Why this is a Letter", built on these points
without naming or comparing to that paper:
- one conceptual change (the borrowed identity) with consequences across a field (every Clifford level, every output type);
- a headline theorem a reader can use (n >= T-count + 1, proved from one line);
- a concrete, checkable payoff for current architectures (cultivation plus one d=2 factory; every class matched or beaten);
- the specialist detail (search stages, catalogues, proofs) in End Matter and Supplement, the main text readable by
  the broad quantum-computing readership.

## I. Questions for you before I start

1. Title: the proposed one, an alternative, or keep the current one?
2. A second figure (four-qubit parent circuit, single column): yes or no?
3. Move the bound proof into the main text (≈ 120 words): yes or no?
4. Keep Fig. E1 (literature recovered) in the End Matter, or move it to the Supplement?
