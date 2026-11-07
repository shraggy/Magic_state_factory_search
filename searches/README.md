# Searches

The cascaded search for distance-2 borrowed-identity factories, the enumeration of exact output classes, and the
checks and figure scripts behind the paper. Run every script from the repository root. See the root
[README](../README.md) for the map from paper claims to scripts and outputs.

## The Clifford hierarchy, briefly

Every script takes a level `l`. A factory consumes noisy copies of the level-`l` magic state for
`Z(theta)`, `theta = pi/2^l`, and Clifford operations, and deposits a diagonal level-`l` state on `k` outputs.

| `l` | weight-1 gate | weight-2 gate | weight-3 gate | weight-4 gate |
|---|---|---|---|---|
| 2 | `S`  | `CZ`  | —     | —      |
| 3 | `T`  | `CS`  | `CCZ` | —      |
| 4 | `√T` | `CT`  | `CCS` | `CCCZ` |

`n` is the number of input magic states, `N` the number of circuit qubits, `k` the number of outputs.

## The cascade

`staged.py` runs the cascade as passes, each reading the best factory per class from the pass before, asking only for
fewer inputs, and writing its own file; `staged_class.py` splits the same passes over cluster tasks (one per class).
The settings are `staged.DEFAULT`; the root README gives the commands of sweep A.

| pass | stage | files | what it fixes and solves |
|---|---|---|---|
| `i_tg` | (i) two-group | `Two_group.py`, via `unified_search.stage_two_group` | outputs and checks each symmetric; binomial sums |
| `i_gs` | (i) group-symmetric | `orbit_group.py` | outputs under `G`, checks symmetric; swept over a group menu, or aimed at one class with `G = Aut(U)` and `G = {e}` (CP-SAT) |
| `i_gh` | (i) check groups | `orbit_group2.py` | as `i_gs` with a group `H` on the checks |
| `ii_dj` | (ii) gluing, disjoint | `symfree_search.py`, via `cascade.stage_gluing_multi` | output rotations fixed to symmetric blocks on disjoint outputs; 2-adic linear solve |
| `ii_eg` | (ii) gluing, entangled | `entangled_gluing.py` (`l=3`), `entangled_gluing_l4.py` (`l=4`) | output rotations fixed to an optimal network of the target; same solve |
| `iii_r` | (iii) restricted SAT | `xor_sat2.py` driven by `chain_iii2.py`; parts from `synthesis_sat.parts` with `det_parts.py` | each rotation's output part drawn from a synthesis of the target; parity equations as XOR constraints (CryptoMiniSat), supports unlimited then `<= 6` |
| `iii_f` | (iii) free SAT | `xor_sat2.py` driven by `chain_iii2.py`; minimality proofs `quick_sat.py` on `sat_certify.py` | all rotations on `N = k+1 .. k+3` qubits, supports `<= 5`, `<= 6`, unlimited (XOR); proofs from z3 over all coefficients |

After the run: `merge_A_v3.py` per stage-(iii) run and `build_staged_best3.py` (stage (iii) per class the best of the runs v4, v5, v6 into `outputs/staged_final_b3`; `build_staged_final.py` gives the earlier v4-only `outputs/staged_final`), `staged_best.py`
(best per class over every cell, re-verified), `check_A.py`/`check_A2.py` (after-run checks), `outside_pipeline.py`
(earlier smaller factories, kept apart), `b_stage_i.py` (stage (i) beyond the class lists).

**SAT-variant comparison** (the table in the root README, "Stage (iii) with different SAT layers"). Each run is one
task per class (`k`, index `I`, total `M`; usage `chain_iii2.py k I M SEEDDIR OUT CHAIN`, with `OUT` below the sweep-A
output directory used as `SEEDDIR`), seeded with the best `n` after `ii_eg`, then merged into per-form pass files:

```
python searches/chain_iii2.py k I M OUT outputs/solver_v2/staged_A_v4 xs:None:2e6,xs:6:2e6,x:5:2e6,x:6:2e6,x:None:2e6
python searches/chain_iii2.py k I M OUT outputs/solver_v2/staged_A_v5 xsb:8:None:2e6,xsb:8:6:2e6,x:5:2e6,x:6:2e6,x:None:2e6
python searches/chain_iii2.py k I M OUT outputs/solver_v2/staged_A_v6 xs:None:2e6,xs:6:2e6,x:5:2e6,x:6:2e6,x:None:2e6,xsb:8:None:2e6,x2:5:2e6,x2:6:2e6,x2:None:2e6
python searches/merge_A_v3.py outputs/solver_v2/staged_A_v4    # likewise _v5, _v6
python searches/sat_variants.py          # -> outputs/checks/sat_variants.json (per run: layers, CPU-h per form, fewer, inputs fewer; best of three; without SAT)
python searches/build_staged_best3.py    # -> outputs/staged_final_b3 (per class the best of v4, v5, v6; summary.json: credit and tie-break counts)
python searches/staged_best.py           # -> outputs/staged_best.json
python searches/stage_contributions.py   # -> outputs/stage_contributions_l3.json, Table I and Table SII (outputs/checks/stage_table.{json,tex}; also
                                         #    Magic_State_Factories_tex/MSD_manuscript/stage_table.tex when that directory exists)
python searches/stage_contributions.py --out-dir DIR   # -> the same three files in DIR only, nothing under outputs/ or the manuscript
python searches/text_counts.py           # -> outputs/checks/text_counts.json (prose counts)
```

Older drivers still used as libraries: `cascade.py` (stage functions), `cascade_class.py` (exact class lists and
labellers), `unified_search.py` (`is_factory`, `distance`, `inputs`, `stage_two_group`, `quirk_url`). `large_unified.py`
with `merge_large_unified.py` runs the cheap forms of stages (i) and (ii) over the wider range (`outputs/large_unified_*`);
`merge_large_unified.regenerate_circuit` rebuilds any of its circuits from its parameters.

Exact output classes:

- `target_classes.py`: number of genuine `k`-output classes per `(l, k)` by Burnside's lemma, checked by brute force
  where feasible; writes `outputs/target_classes_l{L}.json`.
- `enumerate_targets.py` (direct orbit table), `enumerate_targets_fiber.py` (`l=3, k=6`, fibered over the cubic part),
  `enumerate_targets_2adic.py` (`l=4, k=4`, fibered over the even elements): one representative per class and an exact
  labeller; write `outputs/targets/l{L}_k{k}.json`.

Commands for every regime are in the root README ("One search, three regimes").

## Running SAT stages reproducibly

1. **Deterministic limits only, fresh process per query.** z3 queries have a step limit (`rlimit`), CryptoMiniSat a
   conflict limit; every query runs in its own process and every synthesis in a fresh z3 context (`det_parts.py`), so
   an answer does not depend on the queries before it. No query has a wall-clock limit (`sat_certify.certify(time_limit=None)`). A query that reaches
   its limit proves nothing, and the class keeps its best factory.
2. **One worker per core you actually get.** On a cluster request `W` cores on one host (`bsub -n W -R
   "span[hosts=1]"`) and pass `--jobs W`; per-class tasks use one core each.
3. **Test on a known cell first**, e.g. `python searches/staged.py --L 3 --k 3`, which must give `[[8,3,2]]` for `CCZ`.
4. **Gap to the bound.** Every factory has `n >= T-count + 1` (Corollary 3); a gap of 0 proves it smallest at any `N`,
   and restricted SAT and the seeded free-SAT form skip such classes (the unseeded coefficient form runs on every class).

## Checks and analysis

| file | what it checks or computes | output |
|---|---|---|
| `bi_condition_check.py` | borrowed-identity condition (Theorem 2) against brute-force phase tables: exhaustive at small `(l, N)`, random beyond | printed |
| `factory_condition_check.py` | the factory form, Eq. (4), against brute force, on factories and non-factories at `l=2,3,4` | printed |
| `closed_form_checks.py` | binomial identities, Theorem 4, and the closed-form families | `outputs/checks/closed_forms.json` |
| `distance3_pair_check.py` | the pair condition for `d>=3` against `unified_search.distance` | printed |
| `parent_cuts_check.py` | every cut of the parent circuits of Fig. 1, by simulation | printed |
| `families_table.py` | Table SI, every family member looked up in the larger catalogue | `outputs/checks/families_table.json` |
| `lostv1_checks.py`, `lostv1/` | first-submission material restored in the SM, each claim re-verified (`is_factory`, brute force, independent scripts) | `outputs/checks/lostv1_checks.json` |
| `v1_two_group_sweep.py` | the first-submission two-group sweep (21,920 valid configurations) | `outputs/checks/v1_two_group_sweep.json` |
| `stage_contributions.py` | Table I: per-stage reach at `l=3` | `outputs/stage_contributions_l3.json` |
| `class_table_l3.py` | per `l=3` class: exact T-count, cascade `n`, synthillation with optimal matrices | `outputs/class_table_l3.json` |
| `synthillation_compare.py` | builds the synthillation template (Campbell-Howard Thm 1, Table I) for a target; library for the next two | — |
| `synthillation_strongest.py` | strongest synthillation per class, padded matrices included, each checked with the BI condition | `outputs/synthillation_strongest_l3.json` |
| `error_coefficients.py` | output-error coefficient `c` (library, used by `staged_best.py` and `check_A.py`) | — |
| `nullity_l3.py` | stabilizer nullity of every `l=3` class and the bound it gives on conversion routes | `outputs/nullity_l3.json` |
| `route_compare_l3.py` | achievable distil-then-synthesize routes (catalysed `4 tau`; per-gate distillation and teleportation; CS/CCZ batched through `[[6m+2,3m,2]]`) against the factory | `outputs/checks/route_compare_l3.json` |
| `route_kappa_l3.py` | exact leading-order `kappa` and first-order acceptance of the three routes against the factory (2 worker processes) | `outputs/checks/route_kappa_l3.json` |
| `route_class_match_l3.py` | factory output and route label in the same exact class, 766 classes | `outputs/checks/route_class_match_l3.json` |
| `route_table_l3.py` | counts of the factory against the three routes (inputs, `kappa`, rejection) overall, on the 330 and by `k`, for the route table (Table SIV of the Letter, Table II of the PRX Quantum version); the top-level sections take each route on the diagonal representative only (659/47/60, 314/11/5, 602 of 645: a weaker variant, not printed in the paper); asserts agreement with the route files; section `full_clifford` (printed in the paper: 608/62/96, 291/23/16, 567 of 645) takes each route as the cheapest over the full-Clifford group of the class and adds the lambda route as a fourth, with the known-answer assert CS01,CS02 -> 8 inputs, `kappa` 28 | `outputs/checks/route_table_l3.json` |
| `syn_vs_routes.py` | synthillation alone against the same full-Clifford routes: per diagonal class (591/53/122, 85.5% of the saving), synthillation minimised over the full-Clifford class (613, 91.4%), and on equal footing over the 622 full-Clifford classes with all three minimised (521 and 529, 88.3%) | `outputs/checks/syn_vs_routes.json` |
| `syn_c_compare.py` | output-error coefficient `c` of the factory against every valid synthillation variant | `outputs/checks/syn_c_compare.json` |
| `synthillation_bound.py` | the counts Theorem 1 of synthillation allows, case 2 of its Table I, margins of the cascade | `outputs/synthillation_bound_l3.json` |
| `clifford_correction.py` | the free diagonal gate `D` (Clifford correction) each factory leaves before readout | `outputs/checks/clifford_correction.json` |
| `class_invariants_l3.py` | gate types and products counted on the class, not the label | `outputs/checks/class_invariants_l3.json` |
| `k78_class_count.py` | exact number of `CCZ`-layer classes at `k<=8` (Burnside over `GL(k,2)`) and completeness of the `k=7,8` list | `outputs/checks/k78_class_count.json` |
| `coset_certify.py` | exact smallest `n` on `N` qubits at `l=3` (minimum weight of a punctured Reed-Muller coset) | printed |
| `tcount7.py` | exact T-count by capped coset decoding (`exists_below2`); designed for up to 8 qubits, but no `k=8` exact T-count has been obtained (used at `k=7`) | printed |
| `isotropy_floor.py` | proved lower bound on the T-count of a CCZ layer (Mazumder, isotropy floor) | — |
| `frontier_k7.py`, `k7_candidates.py` | CCZ-layer candidates at `k>=7`: heuristic upper bound and isotropy floor; their layers are among the targets of `k78_targets.py` | `outputs/k7_candidates/`, `k8_candidates/` |
| `b_stage_i.py`, `k789_ii_dj.py` | beyond the class lists at `l=3`, `k=7..9`: stage (i) sweeps, stage (ii) disjoint gluing | `outputs/staged_B/`, `outputs/k78_pipeline/full/l3_k{7,8,9}/ii_dj.json` |
| `k78_targets.py`, `k78_staged_B_ccz.py`, `k78_orbits.py`, `k78_eg.py`, `k78_pass.py`, `k78_exact7.py`, `k78_tau_upper.py`, `k78_verify_pass.py`, `k78_xor.py`, `k78_xor_pass.py`, `k78_large_ccz.py`; self-tests `k78_tcount7_selftest.py`, `k78_isd_selftest.py` | the cascade on a target list of `CCZ` layers at `k=7, 8`: targets, entangled gluing (ii), exact `k=7` T-count, fast free SAT (iii), merge and re-verification | `outputs/k78_pipeline/` (`REPORT.md`) |
| `export_beyond_csv.py` | catalogue of the regimes beyond the class lists, with Quirk links | `outputs/beyond_class_lists_with_quirk.csv.zip` |

## Figures

| file | figure | output |
|---|---|---|
| `plot_parent.py` | Fig. 1, parent circuits and their cuts | `figures/parent_l3.{pdf,png}` and the manuscript `FIGURES/parent_l3.pdf` |
| `plot_frontier.py` | Fig. 2, every `l=3` class against the bound and synthillation | `figures/frontier_l3.png` and the manuscript copy |
| `plot_literature.py` | Fig. 4, factories found and the literature factories they recover | `figures/literature_recovered.png` and the manuscript copy |

## Notebooks

- `sequential.ipynb`: sequential borrowed identities (the construction behind the malleable chains), the worked
  examples of App. S2.2; the claims are checked by `sequential_check.py`.

Scripts and results that nothing in the search (any of its three regimes) or in the paper uses are in `../trash/` (gitignored).

**Caveat when extending any search:** a term with coefficient nonzero mod `2^l` can still be Clifford. Check
genuineness with `classification/classify.py`, not a bare "coefficient != 0" test.
