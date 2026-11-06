# Supplement restructure: assessment and plan (no edits made yet; Fig. 1 header removed)

The Supplement is currently 7,644 words over pages 7–24. The main text was restructured around Theorem 2 (BI condition) and
Corollary 3 (bound), but the Supplement still follows the *old* paper's order and has three problems, as you said:
1. **Order is wrong.** Malleable/sequential circuits come before the search in the main text, but S3 (sequential) sits *after*
   S2 (two-group). Two-group is now a *stage of the cascade*, yet it is its own section, apart from S4 (cascade).
2. **Redundant.** The same material appears more than once:
   - Stage descriptions appear in S4's intro, S4.1 (each stage again), Table SII (an old three-stage table), and S2's opening.
   - Output classification by degree (S2, "Output classification") is obsolete: every output is now classified exactly (S5),
     and its "degree-d layer" text duplicates S5's intro.
   - Table SII describes the *old* three-stage search ("all three stages", "SAT (optional)", 591 factories), contradicting
     Table I of the main text. It is also *ridiculously large*: a single-column table squeezed by \resizebox.
   - Fig. S4(b) "Factories by search stage" repeats Table SII graphically for the old stages.
   - "Distance saturation" (S1.2) and "Search complexity" (S2) repeat main-text or End-Matter remarks.
   - Two example verifications (S1.4) re-derive [[8,3,2]] already worked in the main text (Eq. 3) and in S3.
3. **Misplaced.** S1.3 (proof of the two-group lemma) belongs with the two-group condition; Fig. S3 (literature recovered)
   is a results figure but sits in the method section.

## Proposed structure (each section answers one main-text pointer)

**S1. The borrowed-identity condition** (proofs behind Theorem 2 and Corollary 3)
- S1.1 Proof of Theorem 2 and the factory condition (moved from S4: the parity expansion, the Clifford threshold e_t, the
  numerical checks). *Labels kept:* `eq:monomial-condition`, `app:cascade` → now points here, see label notes.
- S1.2 Relation to the Heisenberg picture (current S1.1, Proposition 6, unchanged).
- S1.3 The bound on every catalogue factory (current S4.4, 95 words).

**S2. Symmetric and sequential circuits** (Theorem 4, malleability) — current S1.2 + S3, merged and *before* two-group
- S2.1 Proof of Lemma 7 (Dicke phase), current S1.2 minus "Distance saturation" (moved to one sentence in S2.2).
- S2.2 Sequential construction and Theorem 11 (current S3.1), Figs. S1–S2 (parent circuits, kept).
- S2.3 Worked examples (current S3.2), which already verify [[8,3,2]]…[[15,1,3]] step by step; the duplicate [[8,3,2]]
  check in S1.4 is deleted. Catalytic conversions and gate lists kept (Referee B asked for them).
- S2.4 Sequential search results (current S3 "Results", 1 paragraph) and Fig. S4(a) only.

**S3. The cascaded search** — current S2 (two-group) merged into current S4, in stage order
- S3.1 Overview: one paragraph, the stage list once, pointing to Table I; the obsolete Table SII is **removed** and its
  still-used numbers (stage times, larger-range factory counts) go into one compact two-column-width table, Table SII, see below.
- S3.2 Stage (i), two-group: the two-group condition with Lemma 8 and Theorem 9 *and* the proof of Lemma 8 (current S1.3),
  search procedure, sign assignment, runtime; Table SI (two-group catalogue, kept, already table*). The obsolete
  "Output classification" paragraph is replaced by one sentence pointing to S4.
- S3.3 Stage (ii), gluing; lower-level rotations (current S4.1 text).
- S3.4 Stage (iii), exact target classes; Table SIII (current, kept).
- S3.5 Stage (ii′), gluing on entangled support.
- S3.6 Stages (iv)–(v), SAT (restricted and unrestricted, current S4.1 text, merged into one subsection).
- S3.7 Coverage and completeness; Table SIV (coverage, kept; current S4.3).
- S3.8 Clifford padding (current S4.1 paragraph, moved here as the last step).
- Fig. S5 (stage scaling) kept in S3.1.

**S4. Factories found** (current S5)
- S4.1 Smallest factory per output type + Tables SV(a–c) (kept); "Exact minimality" paragraph kept.
- S4.2 Comparison with prior constructions (kept).
- Fig. S3 (literature recovered) moves here, as it is a results figure.

## Table SII (the oversized one)

- **Remove** its "main range" block: obsolete, describes the old three-stage search.
- **Keep** what the text cites: larger-range factories per level, largest n, stage times. Make it a compact table* across
  both columns of the one-column Supplement, five rows × four columns:

  | | l=2 | l=3 | l=4 | time |
  |---|---|---|---|---|
  | two-group factories (N ≤ 18) | 1,457 | 1,242 | 510 | 0.4–4 s (N ≤ 11); minutes (N ≤ 18) |
  | gluing factories | 27 | 49 | 66 | 1–300 s |
  | largest n | 262,142 | 262,142 | 65,534 | |
  | classes reached, k ≤ 6 (k ≤ 4 at l=4) | 9 | 766 | 4,167 | |

  The counts are from `outputs/large_unified_factories.csv` (by stage, checked) and the coverage files.
  The timing numbers come from the existing table / `outputs/stage_timing` and stay as they are.
- Since the Supplement is single-column (\onecolumngrid), the \resizebox to \columnwidth there is what made it look
  huge: the new table uses natural width.

## Removed as redundant (≈1,100 words)

- S1.4 "Example verifications" ([[8,3,2]] duplicate; [[20,4,2]] two-group example kept, moved into S3.2 as the worked
  example of Theorem 9, since Referee C's table checks use it).
- S2 "Output classification" by degree (obsolete; one sentence remains).
- S4 intro stage list (duplicated by S3.1).
- Fig. S4(b) "Factories by search stage" (old stages; Table I and Table SII carry this). Fig. S4(a) (sequential) kept.
- "Distance saturation" paragraph → one sentence in S2.2.

## Labels and cross-references

Every \label stays on the same object; only section order changes, so \ref's resolve. The letter cites App. S-numbers
by \ref through xr-hyper, so they update automatically on recompile. Hard-coded names in the letter ("App. S3",
"Appendix F" in referee quotes) will be checked after compile.

## Checks after the restructure

- Compile main and letter: zero errors, no undefined refs.
- Every table and figure fits its width; render each Supplement page and look at it.
- Supplement word count before → after.
- The letter's pointers (App. numbers) re-read against the new numbering.

## Questions

1. Delete Fig. S4(b) (old stage scatter)? I recommend yes; Table I and the new Table SII cover it.
2. Keep the [[20,4,2]] worked example (moved to the two-group stage)? I recommend yes, it is the one fully worked
   two-group check.
