# k=7, 8 in the pipeline (beyond the class lists, CCZ-layer targets)

Saved by the main agent from the k78 agent's final hand-off (2026-11-03); the agent could not write report files.
Status: complete (all cluster jobs ended, results collected in all_isd/ and full/); the FINAL section below replaces
the pending parts. Per listed layer: 142 of 142 at k=7 and 634 of 4,043 at k=8 are on the bound n = tau+1 (582 decided without the factory, 52 more by its own bound tau <= n - r; full/on_bound.json)
(full/l3_k{7,8}/ii_eg.json, rows with on_bound true).

## Why k=7, 8 were outside the pipeline
Two code reasons, neither in the theory:
1. **No class list.** `staged.ctx` (searches/staged.py:36-40) uses `cascade_class.classes_of`, which exists only for
   k <= 6; `label_fix` (staged.py:47-53) credits a factory only to a listed class, so nothing is kept at k >= 7.
2. **No exact synthesis.** `_eglue` (cascade.py:52-81) calls `synthillation_compare.opt_A`, which enumerates all 2^dim
   words of the Reed-Muller coset (synthillation_compare.py:43): dim 64 at k=7, 163 at k=8.

## The pass
`searches/k78_eg.py`, no existing file edited (so ii_eg at k <= 6 and every sweep A seed are unchanged):
- the same `entangled_gluing.glue`, r = 1..4, both Clifford forms;
- targets from a list (CCZ layers), one per class of a GL(k,2)-invariant fingerprint (`searches/k78_orbits.py`);
- synthesis `frontier_k7.bound` (fixed seed, restarts), optionally lowered by information-set decoding (`--isd`,
  fixed seeds, every word checked to lie in the coset);
- a one-check candidate (all synthesis columns through one check);
- `--exact7`: exact T-count at k=7 by `tcount7.exists_below2`;
- every factory checked: `is_factory` (monomial and brute force on all 2^N states), distance 2, `deposits`.
`searches/k78_pass.py` merges a run into a staged.py-schema pass file (previous pass = sweep B stage (i) + ii_dj,
credit to the earlier stage on ties).

## Results (files under outputs/k78_pipeline/)
**Reproduction** (targets_paper62.json; paper62/, paper62_isd/; re-verified independently by the main agent):
57 of 57 on-bound factories at the paper's n, N (1 at k=7, n=16, N=8; 56 at k=8, n=18, N=9), all distance 2; 5 of 5
Table SIII(b) rows reproduced, 4 smaller:

| row (paper n) | pipeline n, N | T-count | how proved | file |
|---|---|---|---|---|
| k=7, 20 | 18, 8 | 17 exact | `exists_below2`; weight-17 synthesis from ISD | exact7_paper62_isd.json, paper62_isd/ |
| k=7, 22 | 18, 8 | 17 exact | same | same |
| k=8, 20 | 18, 9 | 17 exact | isotropy floor 17 = ISD synthesis 17 | tau_upper_paper62_isd_row_k8_n20.json |
| k=8, 22 | 20, 9 | 17-19 | undecided | tau_upper_paper62_isd_row_k8_n22.json |
| k=8, 24 | 20, 9 | 17-19 | undecided | tau_upper_paper62_isd_row_k8_n24.json |

Three more factories on the bound: two at k=7, [[18,7,2]]_8 (the two rows are one class up to CNOTs), one at k=8,
[[18,8,2]]_9. T-count upper bounds from the factories' syntheses, checked as polynomials and in the coset
(searches/k78_tau_upper.py, tau_upper_paper62.json: 62 of 62 valid). Laptop vs cluster: identical n, N, bound and
gates on the 58 of 62 paper targets whose cluster tasks finished (paper62_isd_cluster/); 4 not compared.

**Errors made and fixed:** (1) glue's greedy kernel reduction (`symfree_search.reduce_support`) missed the one-check
form for front_k8_44 (n=20 instead of 18; first run kept in paper62_glue_only/); fixed by the one-check candidate.
(2) ISD returned no word when the coset base already met the target weight; fixed. (3) Job 1861551 failed on all 400
tasks: the cluster had the laptop's arm64 `searches/csrc/capped_minweight` ("Exec format error"); an x86_64 build from
the same .c is in `searches/csrc_x86_64/` (cluster), picked on non-Darwin hosts; resubmitted as job 107577.

**Seeds.** Sweep B stage (i): 317 CCZ-layer outputs at k=7, 202 at k=8 (staged_B_ccz.json, searches/k78_staged_B_ccz.py);
B's labels are frame-dependent (317 outputs, 33 labels, 3 fingerprint classes), so matching is by fingerprint. B's
best per class (staged_B_ccz_best.json, re-verified): k=7 n = 16, 16, 22; k=8 n = 18, 18, 20, 20, 20, 22; N = k+1.

Stage (ii) disjoint gluing, new sweep B pass (searches/k789_ii_dj.py, `cascade.stage_gluing_multi`, rmax=4;
full/l3_k{7,8,9}/ii_dj.json; k=7,8 re-verified by the main agent):

| k | partitions solved and verified | CCZ layers | larger-sweep gluing rows reproduced | CPU-s |
|---|---|---|---|---|
| 7 | 8 of 8 | 0 | 8 of 8 | 93 |
| 8 | 10 of 10 | 0 | 10 of 10 | 356 |
| 9 | 12 of 12 (11 distinct outputs) | 1 (three disjoint CCZ blocks) | 12 of 12 | 1354 |

None of the larger sweep's 231 k=7,8 rows is a CCZ layer (searches/k78_large_ccz.py, large_ccz.json; 231 of 231
rebuilt and verified), so for the k=7,8 CCZ classes the running-best seed is stage (i) alone.

**Targets** (searches/k78_targets.py, targets_all.json): 4,185 distinct genuine CCZ layers (142 at k=7, 4,043 at
k=8) from the paper, the candidate files, sweep B and random orbit samples; 6 fingerprint classes at k=7, 20 at k=8
(lower bounds on the number of classes: the fingerprint is not proved injective).

**Known-answer checks:** fingerprints reproduce the pure-CCZ class counts of class_table_l3 at k=3..6 (1, 0, 1, 3);
`exists_below2` agrees with the exact decoder on 38 of 38 random layers (searches/k78_tcount7_selftest.py,
tcount7_selftest.json); ISD reaches the exact T-count on 30 of 30 random layers at k=5,6 with 0 invalid words
(searches/k78_isd_selftest.py, isd_selftest.json); XOR driver: front_k7_0 at N=8 gives n=16 = tau+1.

**Cost:** cluster CPU per target median 71 s (k=7), 88 s (k=8), up to ~350 s at k=8 when tau is undecided.
Reproducing the 62: ~1.5 CPU-h. All 4,185 layers: ~100-200 CPU-h, ~1-1.5 h on 400 tasks. T-count side: k=7 exact
(~7 s per layer); k=8 no exact decision, bounds are isotropy floor and ISD/compilation, on the bound only where they
meet. SAT at k=7,8: feasible at the solver agent's sizes, not measured.

## Remaining (superseded by FINAL)
1. `bjobs -a 107577`; resubmit EXIT/RUNLIMIT tasks with the same tasks_k78_all.txt line (chunks already present are skipped).
2. Copy all_isd/ to outputs/k78_pipeline/all_isd/.
3. `python searches/k78_pass.py outputs/k78_pipeline/all_isd outputs/k78_pipeline/full` -> full/l3_k{7,8}/ii_eg.json.
4. `python searches/k78_exact7.py outputs/k78_pipeline/all_isd` -> exact7_all_isd.json.
5. XOR pass: copy searches/k78_xor.py, k78_xor_pass.py to the cluster; one target per fingerprint class not on the
   bound; N = k+1, k+2; `--wmax 5 --confl 2000000`, then `--wmax 6`; cluster only; merge with
   `python searches/k78_xor_pass.py outputs/k78_pipeline/full <w5dir> <w6dir>` -> iii_f.json (~26 classes x 2 N x 2
   variants, well under 200 CPU-h).
6. Per-class table (pending the merge); the count of stage (i)+(ii-disjoint) targets at k=7-9 and their XOR cost.

## Proposed paper framing (wording for the user and main agent to decide)
"Beyond six outputs no class list exists, so we run stage (ii) on a target list: the CCZ layers on k=7,8 outputs, one
per class of a GL(k,2)-invariant fingerprint. Their output rotations are fixed to the best compilation found, and the
T-count is bounded below by the isotropy bound, or decided exactly at k=7. Where the two bounds meet, the pass reaches
n = tau(U)+1, the form Corollary 3 requires: [N_7] classes at k=7 and [N_8] at k=8, including [[16,7,2]]_8,
[[18,7,2]]_8, [[18,8,2]]_9." N_7 and N_8 wait for the merge; until then the reproducible statement is the 57 factories
plus the 3 new on-bound ones.

## FINAL (2026-11-03; replaces the pending parts above)
All jobs ended: 107577 395 DONE + 5 RUNLIMIT, the 5 rerun as 2006234 (5/5 DONE, table unchanged; they held 52 targets on
the slow hosts cccxc550/552); XOR job 2006780 36/36 DONE. Files: full/l3_k{7,8}/ii_eg.json, iii_f.json. Class = GL(k,2)
fingerprint (k78_orbits.py). Every class best re-verified twice: by the agent (k78_verify_pass.py, brute force) and by the
main agent (is_factory, distance 2, n recomputed): 6/6 at k=7, 20/20 at k=8, 0 invalid.

**k=7: N_7 = 6 of 6 classes on the bound n = tau+1.** k7f0, k7f3 n=16 (tau 15, isotropy floor), from sweep B stage (i);
k7f1 n=18 (tau 17), k7f2 n=20 (tau 19; lowers B's 22), k7f4, k7f5 n=20 (tau 19), from ii_eg, tau exact by
tcount7.exists_below2. Over the target list 142/142 k=7 layers have exact tau and an on-bound factory (exact7_all_isd.json).

**k=8: N_8 = 2 of 20 on the bound** (k8f0, k8f1, n=18, tau 17, from B stage (i)). 20/20 classes reached, 14 not by B.
The other 18: n=20 (f2, f10, f19 from B; f3, f5, f6 from ii_eg), n=22 (f4 from B; f7-f9, f11-f15 from ii_eg), n=24
(f17, f18), n=26 (f16); tau between the floor 17 and n-1, undecided. Note: the files record the upper bound as n (not
n-1) for f4, f6, f11, f12; the factory itself gives tau <= n-1, so those stored bounds are loose, not wrong.

**XOR iii_f on the 18:** 0 lowered. wmax 5: N=9 UNSAT below the cap for 15/18, N=10 UNSAT for 2 (f2, f10); rest UNKNOWN
at 2e6 conflicts. 2.75 CPU-h.

**Cost:** gluing pass k=8 1,696,000 CPU-s (~471 CPU-h), k=7 8,061 CPU-s; ~4 h wall on 400 tasks plus reruns (one
representative per class would be ~1 CPU-h). ii_dj k=7/8/9: 93 / 356 / 1,354 CPU-s. Stage (i)+(ii-disjoint) targets at
k=7-9: ~296 reduced labels (65+8, 85+10, 117+11; frame-dependent upper bounds, B stage (i) partial at k=8,9); XOR on them
estimated ~90 CPU-h, not run.

**Paper rows (Table SIII(b)):** k=7 rows n=20,22 -> 18 (both on the bound); k=8 row n=20 -> 18 (on the bound); k=8 rows
n=22, 24 -> 20 (tau 17-19).

**Proposed framing (superseded: the paper states the per-layer counts above):** "All 6 classes at k=7 and 2 of 20 at k=8 reach n = tau(U)+1, including
[[16,7,2]]_8, [[18,7,2]]_8, [[20,7,2]]_8 and [[18,8,2]]_9; for the other 18 classes at k=8, n is the compiled T-count
plus one and stage (iii) finds nothing smaller."
