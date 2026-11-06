# Outputs

Result files of the search and of the checks behind the paper. `n` is the number of input magic states, `N` the number
of circuit qubits and `k` the number of outputs, as in `[[n,k,2]]`. Every factory here has distance at least 2, is built
from `Z(pi/2^l)` parity rotations only, and was verified against the borrowed-identity condition. Every number quoted
below was read from the file named on the same line.

## Sweep A: the staged cascade on every cell with a class list

| path | from | contents |
|---|---|---|
| `staged_A/l{L}_k{k}/` (`l=3, k=5,6`; `l=4, k=4`) | `searches/staged_class.py` | the large cells, one task per class: `cell/` (sweeps), `cls/` (per-class passes, every factory found with its stage, `n`, `N`, gates, CPU seconds), `free_<form>/` (free SAT per form), and the merged pass files `i_tg.json ... iii_f.json` |
| `staged_pc/l3_k{1..4}/` | `searches/staged_class.py` | the same for `l=3, k<=4` |
| `staged_final_b3/l3_k{1..6}/` | `searches/build_staged_best3.py` | the `l=3` cells as the paper uses them: stages (i)-(ii) of sweep A, stage (iii) (`iii_r`, `iii_f`) per class the best of the three deterministic runs `solver_v2/staged_A_v4`, `_v5`, `_v6`; `attribution.json` (per class: run, layer, form), `forms.json` (per run and layer); `summary.json` (credit per form, its dependence on the tie-break between runs) |
| `staged_final/l3_k{1..6}/` | `searches/build_staged_final.py` | the earlier `l=3` cells with stage (iii) from `staged_A_v4` alone (the v4-only pipeline; its derived outputs are in `archive/v4_published/`) |
| `solver_v2/` | `searches/xor_*.py`, `chain_iii*.py` | stage-(iii) runs: `staged_A_v4`, `staged_A_v5`, `staged_A_v6` (the three deterministic runs whose per-class best is the final pipeline), earlier history-dependent runs `staged_A_v2`, `staged_A_v3`, coverage runs `cov/`, the synthillation recheck `optB_strongest*`; `REPORT.md` |
| `outside_pipeline/` | `searches/outside_pipeline.py` | 4 smaller verified factories the pipeline does not reproduce (1 from an earlier run, 35 with the v4-only pipeline; 3 from the arm64 Test A), with how each was obtained (README) |
| `staged_l4fix/l4_k3/`, `staged_small/l{2,4}_k*/` | `searches/staged.py` | `l=4, k=3`; `l=2, k<=6`; `l=4, k<=2` (one process per cell) |
| `staged_*/_cache/` | `searches/cascade_class.py` | pickled exact class tables (gitignored; rebuilt on first use, up to ~1 h for `l=3, k=6` and `l=4, k=4`) |
| `archive/v4_published/` | copies | the outputs of the v4-only pipeline that the per-class best of v4, v5, v6 replaced (README there) |
| `staged_best.json` | `searches/staged_best.py` | **the smallest factory found for each of the 4,942 classes** (`l=3` from `staged_final_b3`) (`l=2`: 9, `l=3`: 766, `l=4`: 4,167), re-verified, with gates, exact `Z(pi/2^l)`-count `tau`, synthillation `n` at `l=3`, `c`, final stage and the first stage reaching its `n` |
| `staged_A/compare_l3.json` | `searches/staged_compare.py` | per class against every earlier run: `l=3, k=5`: 57 equal, 7 better, 1 worse of 65; `k=6`: 605 equal, 66 better, 8 worse of 679 |
| `staged_A_rerun/` | `searches/rerun_A_tasks.py` | rerun of the classes a limit may have touched; not merged (the verified smaller factories are in `outside_pipeline/`) |
| `earlier_catalogues/` | earlier scripts (now in `trash/`) | the `l=3` catalogue and CSVs before sweep A, kept because they record the classes where an earlier SAT run is still smaller |

Each pass file has `settings`, `seconds` (CPU seconds summed over tasks when `seconds_are_cpu_sum`), `classes`, `reached`,
`new`, `lowered`, `inputs_saved`, and `best` (per class: `n`, `N`, `stage`, `gates`).

## Catalogues (from `staged_best.json`)

| file | from | contents |
|---|---|---|
| `smallest_per_class.csv.zip` | `searches/export_csv.py` | one row per class (4,942), the earlier column layout plus a Quirk link (`clifford_class = L{l}-k{k}-{class}`; padding columns empty) |
| `all_factories.csv`, `all_factories_with_quirk.csv.zip` | same | the same rows without and with the Quirk column |
| `l3_catalogue_k6.csv` | same | the 766 `l=3` rows: `k, spec, tcount, n, N, c, source, syn_strongest_n, gates, quirk, verified` |

Every Quirk link contains exactly its factory: 4,942 of 4,942 (`checks/quirk_check.json`). All 4,942 circuits have
`N <= 14`, so every row has a link.

## Per-class tables at l=3

| file | from | contents |
|---|---|---|
| `class_table_l3.json` | `searches/class_table_l3.py` | per class: exact T-count, `mu`, staged `n`, synthillation with optimal matrices |
| `synthillation_strongest_l3.json` | `searches/synthillation_strongest.py` | per class, every synthillation variant on its Table I case, checked, and the smallest valid one; valid for all 766 |
| `stage_contributions_l3.json` | `searches/stage_contributions.py` | Table I, and totals: fewer than synthillation on 330 of 766 (243 by two, 87 by four), same on 436, more on 0; 5 at the bound (295, 52, 471 with stage (iii) from v4 alone) |
| `checks/text_counts.json` | `searches/text_counts.py` | prose counts no other script writes: entangled / product classes among the fewer-than-synthillation ones, `N=8,9`, gap >= 17, equal `n` with equal `N`, per-form lowerings of stage (iii); `old` = archived v4-only build, `b3` = current |
| `nullity_l3.json` | `searches/nullity_l3.py` | stabilizer nullity per class (`nu = k` for 766 of 766) and the bound on conversion routes (beaten on 1) |
| `entangled_gluing_l3.json`, `entangled_gluing_l4_k*.json` | `searches/entangled_gluing.py`, `entangled_gluing_l4.py` | gluing on entangled support run on its own; the `l=4` files supply `tau` to `staged.py` |
| `pattern_sat/`, `pattern_sat_local/`, `synthesis_sat/`, `cascade/`, `cascade_main_v3/`, `cascade_v4/` | earlier runs of the same stages | the reference runs that `staged_compare.py` compares sweep A against |

## Class lists and the regimes beyond them

| path | from | contents |
|---|---|---|
| `targets/l{L}_k{k}.json`, `target_classes_l{2,3,4}.json` | `searches/enumerate_targets*.py`, `searches/target_classes.py` | one representative per class; class counts (`l=3`: 1, 2, 6, 13, 65, 679; `l=4`: 2, 8, 75, 4,082; `l=2`: 1, 2, 1, 2, 1, 2) |
| `k7_candidates/`, `k8_candidates/` | `searches/k7_candidates.py` (`best.json`, `exact2_20_22.json` inline) | CCZ-layer candidates at `k=7, 8`, among the targets of `k78_pipeline/targets_all.json` |
| `large_unified/`, `large_unified_*` | `searches/large_unified.py`, `merge_large_unified.py` | the wider range: the cheap forms of stages (i) and (ii) (two-group, disjoint gluing), to `N=18` (`l=4`: 16) |
| `staged_B/` | `searches/b_stage_i.py` | stage (i) at `l=3`, `k=7..9`, one file per combination (`k=8,9` partial) |
| `k78_pipeline/` | `searches/k789_ii_dj.py`, `k78_*.py` | stage (ii) disjoint gluing at `k=7..9`; the cascade on a target list of `CCZ` layers at `k=7, 8` (`full/`, `all_isd/`, `REPORT.md`) |
| `beyond_class_lists_with_quirk.csv.zip` | `searches/export_beyond_csv.py` | every factory of the regimes beyond the class lists, with Quirk links (7,832 rows) |
| `literature_refs.json` | `searches/unified_search.py` | literature references keyed `l,k,n,N` (read by `unified_search.py`) |

## Checks (`checks/`)

| file | from | contents |
|---|---|---|
| `check_A/` | `searches/check_A.py` | after-A check of sweep A's own files (`staged_A`: `l=3`, `k=5,6`; `l=4`, `k=4`): 4,826 of 4,826 pass, 4,819 also by brute force; the final `l=3` stage (iii) is checked by `searches/check_A2.py` in `solver_v2/staged_A_v4/checks/`: 766 of 766, 761 by brute force |
| `quirk_check.json` | `searches/quirk_check.py` | 4,942 of 4,942 Quirk links contain exactly their factory |
| `stage_table.json`, `.tex` | `searches/stage_contributions.py` | Table SII |
| `table_sv.json`, `.tex` | `searches/table_sv.py` | Tables SIV(a-c), 70 rows, each re-verified |
| `families_table.json` | `searches/families_table.py` | Table SI: 9 families, 48 members |
| `lostv1_checks.json` | `searches/lostv1_checks.py` (+ `searches/lostv1/*.py`) | first-submission material restored in the SM: v1 Table V entries re-verified (21 of 21), $(2,k-1)$ family 16 of 16, C^{l-1}Z family 16/15/12 of 16/15/12, six-qubit parent cuts, padded [[8,4,2]], [[8,3,2]] output phase, small-N brute force, even weights, relabelling, symmetric distance (244 of 244), CCZ Clifford orbit, catalyst split |
| `v1_two_group_sweep.json` | `searches/v1_two_group_sweep.py` | v1 two-group sweep statistics quoted in the response letter: 7,203 tuples, 107,940 sign patterns, 21,920 valid, 689 groups, 535 non-trivial |
| `catalytic_check.json` | `searches/catalytic_check.py` | `[[10,2,2]]` and `[[11,1,2]]` by statevector simulation; `[[10,2,2]]` output error `9p^2` per output, `13p^2` joint |
| `syn_c_compare.json` | `searches/syn_c_compare.py` | `n`, `c` of the factory and of every valid synthillation variant, per `l=3` class |
| `clifford_correction.json` | `searches/clifford_correction.py` | the free diagonal gate `D` each factory leaves on checks and outputs, per level |
| `route_compare_l3.json` | `searches/route_compare_l3.py` | distil-then-synthesize routes (`4 tau` catalysed; per-gate teleportation) against the factory |
| `k78_class_count.json` | `searches/k78_class_count.py` | Burnside count of `CCZ`-layer classes at `k<=8`; fingerprints of the `k=7,8` target list |
| `class_invariants_l3.json` | `searches/class_invariants_l3.py` | per `l=3` class: single gate type, product on disjoint qubits |
| `closed_forms.json`, `two_group_check.json`, `two_group_distance.json`, `sequential_check.json`, `bi_condition_check.json`, `catalogue_bound_check.json`, `syndrome_rank_check.json`, `best_symmetry_l3.json`, `orbit_group_l3.json` | the scripts of the same names | the condition and distance checks of the paper (see the root README table) |

Large circuits: the `large_unified_l{2,3,4}_circuits.txt.zip` files hold gate lists (`q,q,...:c` means `Z(c*pi/2^l)` on
those qubits); two-group circuits are fixed by `(k, N, s_total, s_O)` and rebuilt exactly by
`merge_large_unified.regenerate_circuit(l, k, N, params)`.
