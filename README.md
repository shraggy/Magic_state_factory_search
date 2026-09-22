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
| [`classification/`](classification/) | Turns a factory's gate list into `(degree, essential_dim, t_count, decomposition)`, and reconstructs a factory's explicit circuit from a catalogue row. **[README](classification/README.md)** |
| [`outputs/`](outputs/) | The factory catalogues (one CSV per level) and raw sweep logs. **[README](outputs/README.md)** |
| [`figures/`](figures/) | Generated plots and Quirk circuit diagrams. **[README](figures/README.md)** |



## Main findings

- **Two-group + symmetry-free together recover hundreds of distance-2 factories** across
  `l=2,3,4` (`outputs/factory_catalogue_l{2,3,4}.csv`), including entangled- and mixed-output
  families no single earlier framework reached at once.
- **The symmetry-free search found two factories new to the literature**: `[[18,4,2]]` (a `T`-to-`CS`
  synthillation family, `[[6m+6,2m,2]]`) and `[[18,5,2]]` (mixed `CS`+`CCZ` output from one factory).
- **The sequential search finds zero new factories** — every genuinely non-Clifford result it
  produces is already in the two-group catalogue. Its value is as the construction behind the
  paper's *malleable circuits* (one parent circuit, many distillation targets), not as a source of
  new records. Details in `searches/README.md`.
- **Checked against a SAT/exhaustive classification's five "SAT-only" `l=3` targets**: one
  (`[[12,3,2]]a`) turns out to already be in the two-group catalogue (a false negative from an
  earlier, looser classification) but isn't a good factory anyway; one (`[[16,6,2]]`) is dominated
  by a cheaper factory already found here; the other three are genuine gaps neither search here
  reaches. Details in `outputs/README.md`.

## An active, extensible tool

This repository is meant to keep being used, not just to reproduce the paper's figures. The
borrowed-identity condition is written generically in `l`, so finding a `Z(π/2^l)`-to-any-`D_l`
magic-state factory at a level or parameter range not yet swept can be done by widening the sweep
(see `searches/README.md`).
