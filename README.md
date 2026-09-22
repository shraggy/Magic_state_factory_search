# Magic State Factory Search

Code and data for:

> S. Singh, C. Gidney, C. Jones, **"Borrowed Identities: Malleable Distillation Factories and a
> Unified Numerical Search,"** [arXiv:2606.28518](https://arxiv.org/abs/2606.28518) (2026).

The paper introduces a *borrowed-identity* condition for magic-state distillation factories — a
circuit only has to act as the identity on **one** input state, not on an entire codespace — which
holds at every level of the Clifford hierarchy and turns factory-finding into a search problem.
This repository implements that search, classifies what it finds, and catalogues the results. See
the paper for the theory; this README covers what's here and what it found.

## Repository structure

| folder | contents |
|---|---|
| [`searches/`](searches/) | The three search methods (two-group, symmetry-free, sequential/malleable) plus the closed-form family's distance verifier. **[README](searches/README.md)** |
| [`classification/`](classification/) | Turns a factory's gate list into `(degree, essential_dim, t_count, decomposition)`, and reconstructs a factory's explicit circuit or a real Quirk link from a catalogue row. **[README](classification/README.md)** |
| [`outputs/`](outputs/) | The factory catalogues (one CSV per level) and raw sweep logs. **[README](outputs/README.md)** |
| [`figures/`](figures/) | Generated plots and Quirk circuit diagrams. **[README](figures/README.md)** |

[`pipeline.ipynb`](pipeline.ipynb), at the root, is a runnable walk-through: search → classify →
reconstruct a circuit (matrix or Quirk link) → independently re-verify it, or look up a factory
already in the catalogue.

## Main findings

The Clifford-hierarchy level `l` fixes which magic state a factory consumes — `θ=π/2^l`, with
`l=3` the level `T` lives at (see the hierarchy table this repo's other READMEs give in full);
every number below is combined across `l=2,3,4`.

- **Two-group and symmetry-free together found 22,052 valid borrowed-identity circuits** — 21,920
  from two-group's full sweep (`l∈{2,3,4}`, `k≤7`, `n≤11`, skip parameters `s_total,s_O≤7` — the
  same sweep and count the paper itself reports) plus 132 from the symmetry-free search's own
  sweep (`l∈{2,3,4}`, `k≤6`, every output block-partition, checks`≤4`) — **combined in under three
  minutes** of search time on a laptop (two-group: well under a second; symmetry-free, which
  solves a targeted linear system rather than just sweeping, is the slower of the two at ~3
  minutes). That's the search itself; classifying every distinct shape found (`degree`,
  `essential_dim`, `t_count`, `decomposition` — a separate, slower step) is what turned those
  circuits into the 539 catalogued rows in `outputs/factory_catalogue_l{2,3,4}.csv`, including
  entangled- and mixed-output families no single earlier framework reached at once.
- **Among the many results are factories new to the literature** — not one isolated example, but
  a whole family at multiple sizes: the `T`-to-`CS` synthillation family `[[6m+6,2m,2]]`, found at
  both `[[18,4,2]]` (`m=2`) and `[[24,6,2]]` (`m=3`); and mixed-output factories with no two-group
  analogue at all, `[[18,5,2]]` (one `CS` block + one `CCZ` block from a single factory) and
  `[[26,6,2]]` (three different output types from one factory). Details in `outputs/README.md`.
- **Checked against a SAT/exhaustive classification's five "SAT-only" `l=3` targets** — out of
  22,052 factories found, only **3** are confirmed gaps neither search here reaches; a 4th
  (`[[16,6,2]]`) is dominated by a cheaper factory already found here, and a 5th (`[[12,3,2]]a`)
  turns out to already be in the two-group catalogue (a false negative from an earlier, looser
  classification). A small, specific miss against an exhaustive/SAT method, not a broad gap — see
  `outputs/README.md` for exactly which three and by how much.
- **The sequential search finds zero new factories** — every genuinely non-Clifford result it
  produces is already in the two-group catalogue. Its value is as the construction behind the
  paper's *malleable circuits* (one parent circuit, many distillation targets), not as a source of
  new records. Details in `searches/README.md`.

## An active, extensible tool

This repository is meant to keep being used, not just to reproduce the paper's figures. The
borrowed-identity condition is written generically in `l`, so finding a `Z(π/2^l)`-to-any-`D_l`
magic-state factory at a level or parameter range not yet swept can be done by widening the sweep
(see `searches/README.md`).
