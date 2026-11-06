# Figures

| file | shows | produced by |
|---|---|---|
| `parent_l3.pdf`, `parent_l3.png` | Fig. 1: parent borrowed identities as support matrices, cut into malleable factories | `searches/plot_parent.py` (cuts checked by `searches/parent_cuts_check.py`) |
| `frontier_l3.png` | Fig. 2: every `l=3` class at (exact T-count, `n`) against the bound `n = T-count + 1` and synthillation | `searches/plot_frontier.py` |
| `literature_recovered.png` | Fig. 4 (End Matter): factories found, with the literature factories they recover | `searches/plot_literature.py` |
| `circuits/qk_*.pdf` | parent circuits drawn from Quirk links, each `Z(pi/2^l)` parity rotation one box (not in the current paper) | `figures/circuits/make_all.py` (`litinski.py`), then `pdflatex` each `qk_*.tex` |

All scripts run from the repo root, e.g. `python3 searches/plot_frontier.py`. Each writes to `figures/` and, if the
manuscript folder is present, to `Magic_State_Factories_tex/MSD_manuscript/FIGURES/`. For the circuits:
`cd figures/circuits && python3 make_all.py && for f in qk_*.tex; do pdflatex $f; done`.
