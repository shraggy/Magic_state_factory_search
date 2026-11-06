# Archive: outputs of the v4-only pipeline (replaced 2026-11-05)

These are the files the paper used while stage (iii) at `l=3`, `k<=6` came from one deterministic run,
`outputs/solver_v2/staged_A_v4` (SAT layers `xs:None, xs:6, x:5, x:6, x:None`), assembled into `outputs/staged_final`
by `searches/build_staged_final.py`. Copied unchanged (timestamps kept) before they were overwritten.

**What replaced them.** Stage (iii) is now, per class, the best of the three deterministic runs v4, v5 and v6
(`outputs/solver_v2/staged_A_v{4,5,6}`), built into `outputs/staged_final_b3` by `searches/build_staged_best3.py`
(per-class credit in `attribution.json`, totals in `summary.json`). `searches/staged_best.py` reads `staged_final_b3` by
default, and every script below was rerun with its default arguments into its normal output path. The earlier build is
reproduced by `python searches/staged_best.py --final-dir staged_final --out staged_best_v4.json`
(checked 2026-11-05: 4,942 of 4,942 rows equal to the archived file; stages (i)-(ii) and all cells other than `l=3` are the same in both builds).

| archived file | written by | what changed in the replacement |
|---|---|---|
| `outputs/staged_best.json` | `searches/staged_best.py` | `l=3` stage (iii) from v4, v5, v6 |
| `outputs/stage_contributions_l3.json` | `searches/stage_contributions.py` | Table I and totals (below) |
| `outputs/checks/stage_table.{json,tex}`, `manuscript/stage_table.tex` | `searches/stage_contributions.py` | nothing (Table SII quotes v4's seconds; the reached counts are the same) |
| `outputs/checks/table_sv.{json,tex}`, `manuscript/table_sv.tex`, `manuscript/circuits_tableVI.tex` | `searches/table_sv.py` | nothing |
| `outputs/checks/gong_compare.json` | `searches/gong_compare.py` | nothing |
| `outputs/outside_pipeline/smaller_than_pipeline.csv`, `README.md` | `searches/outside_pipeline.py` | 35 classes -> 1 |
| `outputs/smallest_per_class.csv.zip`, `all_factories.csv`, `all_factories_with_quirk.csv.zip`, `l3_catalogue_k6.csv` | `searches/export_csv.py` | the `l=3` rows |
| `outputs/checks/literature_points.json`, `figures/literature_recovered.png`, `manuscript/FIGURES/literature_recovered.png` | `searches/plot_literature.py` | `l3_sizes.tie_only` 21 -> 20 (not quoted) |
| `outputs/checks/all_found_l3.json` | `searches/all_found_l3.py` (now `--runs v4,v5,v6` by default) | 3,114 -> 3,188 factories; same 228 sizes, same 19 of 68 Nezami--Haah sizes |
| `outputs/checks/syndrome_rank_check.json` | `searches/syndrome_rank_check.py` | summary unchanged (1,490 checked, 0 violations, 729 equalities); one `l=3` row has a different factory |
| `outputs/nullity_l3.json` | `searches/nullity_l3.py` | 69 of 766 rows (their `n`); quoted counts unchanged (766 of 766 with nullity `k`; route bound above the factory for 1, equal for 2, below for 763) |
| `outputs/checks/free_sat_proofs.json`, `catalogue_bound_check.json`, `quirk_check.json`, `sat_variants.json` | the scripts of the same names | nothing (rerun, byte-identical) |
| `figures/frontier_l3.png`, `manuscript/FIGURES/frontier_l3.png` | `searches/plot_frontier.py` | manuscript copy replaced by the same script run on the new data: panel (b) `k=6` bar 273 -> 308, points of panel (a); `figures/frontier_l3.png` not yet regenerated |
| `searches/*.py` | | the scripts as they were before their defaults were switched |

**Headline numbers, old -> new** (`outputs/stage_contributions_l3.json` here and in `outputs/`; inputs fewer in total
from `outputs/checks/sat_variants.json`, `v4` and `best_of_three`): fewer inputs than synthillation on 295 -> 330 of 766
classes (by two: 243 -> 243; by four: 52 -> 87; `k=6`: 273 -> 308, `k=5`: 22 -> 22), the same on 471 -> 436, more on
0 -> 0; inputs fewer in total 694 -> 834; on the bound 5 -> 5; gap to the bound of 17 or more 355 -> 320.
Table I "first" column: (i) 177 -> 176, (ii) 204 -> 188, restricted SAT 293 -> 321, free SAT 92 -> 81.
