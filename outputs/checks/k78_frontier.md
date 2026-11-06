# k = 7, 8 factories of sweep B against synthillation (searches/k78_frontier.py)

## Known-answer checks (all run by this code in this run)
- synthillation, 766 classes (k<=6): syn_hi with exact A equals synthillation_strongest_l3.json best_n on 766 of 766; syn_lo <= best_n on 766 of 766 and syn_lo = best_n on 766 of 766; every lower bound <= exact tau on 766 of 766; k essential outputs on 766 of 766; permutation-canonical form unchanged by a random relabelling on 766 of 766.
- least Delta per Table I head (exhaustive, delta_min): {'linear A0B0': 8, 'linear A0B1': 10, 'linear A1B0': 9, 'linear A1B1': 11, 'no linear A0B0': 0, 'no linear A0B1': 2, 'no linear A1B0': 1, 'no linear A1B1': 3}.
- T-count code: 20 of 20 k<=6 classes of staged_best.json; 142 of 142 layers of exact7_all_isd.json (14 needed the exact search).
- restriction lower bounds (run only at k >= 7): k=6 classes with T or CS terms in a random CNOT frame at k=7, 8: every bound <= tau on 120 of 120 (k=7 60 of 60, k=8 60 of 60); a restriction bound reaches tau on 120 of 120. A test, not a proof (the proof is Lemmas 10, 11 of Mazumder arXiv:2610.01024).

## Factories
6722 of 6722 re-verified (is_factory, distance 2, n), each also after relabelling its outputs into the permutation-canonical frame.

## Totals per k (targets = reduced polynomials up to qubit permutation)

| k | raw-frame keys | targets | dropped (< k essential) | kept | CCZ layers | fewer | same (every Table I) | same (built only) | fewer_than_built | more | open | on bound | lo = hi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | 362 | 337 | 0 | 337 | 124 | 14 | 300 | 0 | 0 | 23 | 0 | 124 | 337 |
| 8 | 4394 | 4019 | 2 | 4017 | 3713 | 0 | 505 | 0 | 20 | 60 | 3432 | 442 | 516 |

Verdicts split: CCZ layers k=7 {'same': 124}; k=8 {'same': 442, 'open': 3271}

Other targets k=7 {'same': 176, 'more': 23, 'fewer': 14}; k=8 {'open': 161, 'same': 63, 'more': 60, 'fewer_than_built': 20}

CCZ layers per GL(k,2) class (fingerprint classes; at k=7, 8 their number equals the Burnside count of outputs/checks/k78_class_count.json, 6 and 20, so they are the classes): k=7 6 classes {'same': 6}; k=8 20 classes {'same': 2, 'open': 18}

Other targets are not reduced under GL(k,2): 'more' targets sharing (mu*, linear terms, Thm 25, hi, lo) with a target of smaller n (possible CNOT-frame artefacts, not proved): k=7 20 of 23; k=8 46 of 60

CPU s: {'known answers: synthillation, bounds (766)': 3849.1, 'known answers: 20 k<=6 classes': 8.1, 'known answers: 142 exact7 layers': 977.0, 'known answers: restriction bounds (120)': 518.6, 'factory verification': 524.3, 'targets': 21970.3} (total 27847.6); wall 10778.7 s.

Verdicts. fewer: n < syn_lo, i.e. fewer inputs than every Campbell-Howard Thm 1 / Table I construction (n = col(A) + 2 col(B) + Delta, A a phase-polynomial synthesis of U, B a Lempel factor); this does not bound their general quasitransversal G-matrix protocols (Thm 4), a family that includes our own factories. same (every Table I): n = syn_hi = syn_lo with exact tau. same (built only): n = syn_hi with exact tau, but syn_lo < n. fewer_than_built: below the construction built here only (an upper-bound comparison). more: n > syn_hi. open: otherwise. T-count bounds hold for Hadamard-free {CNOT, T} circuits with clean ancillas.
