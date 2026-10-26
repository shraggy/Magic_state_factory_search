# Figures

| file | shows | produced by |
|---|---|---|
| `new_factories.png` | Fig. 1: each new factory against the best prior construction for its target | `searches/plot_new_factories.py` |
| `literature_recovered.png` | Fig. E1: every factory found (gray), with the literature factories it recovers | `searches/plot_literature.py` |
| `stages.png` | Fig. S1: factories by search stage (color = level `l`, shape = stage) | `searches/plot_stages.py` |
| `stage_scaling.png` | Fig. S2: wall-clock time of each search stage against `N`, one line per `k` | `searches/plot_stage_scaling.py` |
| `sequential.png` | Fig. S3: sequential search results, `l=2,3,4` | `searches/plot_supp_figures.py` (or `searches/sequential.ipynb`) |
| `circuits/qk_*.pdf` | Figs. S4-S5: annotated parent circuits, drawn from the paper's Quirk links (each `Z(pi/2^l)` parity rotation one box, angle on top) | `figures/circuits/make_all.py` (`litinski.py`), then `pdflatex` each `qk_*.tex` |

All scripts run from the repo root, e.g. `python3 searches/plot_stages.py`. For the circuits:
`cd figures/circuits && python3 make_all.py && for f in qk_*.tex; do pdflatex $f; done`.
