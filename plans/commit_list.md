# Files for the user to commit (2026-11-01, night before submission). Claude did not commit.

The manuscript folder `Magic_State_Factories_tex/` is gitignored in this repo; its files are listed separately.

## Code (new)
searches/entangled_gluing.py, searches/entangled_gluing_l4.py      stage (ii'): gluing on entangled support
searches/pattern_sat.py, searches/synthesis_sat.py                   stage (iv): SAT restricted by App. S4
searches/synthillation_strongest.py, searches/class_table_l3.py      synthillation built from Table I, all classes
searches/merge_l3_best.py, searches/build_l3_catalogue.py           best per class, verified catalogue with circuits
searches/nullity_l3.py                                               stabilizer-nullity bound (model c)
searches/tcount7.py, searches/isotropy_floor.py, searches/k7_candidates.py, searches/frontier_k7.py,
searches/frontier_search.py, searches/classify_frontier.py           T-count frontier, k >= 7 candidates
searches/plot_frontier.py                                            Fig. 1
searches/sat_rerun.py, searches/l4_candidates.py, searches/add_error_coefficient.py,
searches/measured_parts.py, searches/synthesis_batched.py, searches/synthesis_l4.py, searches/synthesis_routes.py

## Code (modified)
searches/cascade.py (stages ii', iv, v; --rlimit), searches/quick_sat.py (rlimit), searches/sat_certify.py
(rlimit, output_parts), searches/min_synthesis.py, searches/plot_new_factories.py, searches/README.md

## Data (new)
outputs/l3_catalogue_k6.csv, outputs/l3_best_per_class.json, outputs/class_table_l3.json,
outputs/synthillation_strongest_l3.json, outputs/nullity_l3.json, outputs/entangled_gluing_l3.json,
outputs/entangled_gluing_l4_k{1..4}.json, outputs/pattern_sat/ (766 files), outputs/pattern_sat_local/,
outputs/synthesis_sat/, outputs/k7_candidates/, outputs/k8_candidates/, outputs/k9_candidates/,
outputs/l4_candidates/, outputs/cascade/*.json (new cells), outputs/frontier_l3.json, outputs/frontier_l4.json,
outputs/frontier_search_*.json, outputs/l3_k6_best_merged.json, outputs/min_synthesis_l3_k*_c2t*.json,
outputs/cost_c2t_l3.json, outputs/synthesis_*_l3.json, figures/frontier_l3.png

## Data (modified)
outputs/README.md, outputs/*.csv (c column), outputs/fig1_l3_rows.json, figures/new_factories.png, .gitignore

## Plans (optional to commit)
plans/finish_paper_plan.md, plans/context_for_review.md, plans/commit_list.md

## Manuscript folder (gitignored here)
Magic_State_Factories_tex/MSD_manuscript/{main.tex, referee_response.tex, cover_letter.tex, distill.bib,
circuits_tableVI.tex, FIGURES/frontier_l3.png} and the three PDFs. Backups of the pre-edit versions:
~/msf_scratch/*_before_final_2321.*

## Not to commit
.vscode/
