# The lambda(U) route (Hadamards, Pauli measurements, feedforward; no catalysts), l=3, k <= 6

Script: `searches/lambda_route.py`. Data: `outputs/checks/lambda_route.json`. Needs `outputs/checks/full_clifford_classes.json` (searches/full_clifford_classes.py).

## Method

```
The lambda(U) route at l=3, k <= 6: distil T (or CS, CCZ) states with the distance-2 factories of the catalogue, then
convert them into the target with Clifford gates, Pauli measurements and feedforward, no catalysts; against our factory.

lambda(U) = fewest |T> states from which stabilizer operations (Cliffords, Pauli measurements, feedforward) make |U>,
without catalysts. lambda(U) <= tau(U) (CNOT+T synthesis). What we use from Beverland, Campbell, Howard, Kliuchnikov,
arXiv:1904.01124 (full text read):
  (a) stabilizer nullity nu is monotone under stabilizer operations, measurements included, and additive (Sec. 2.1,
      Prop. 4.1), so lambda(U) >= nu(U). We compute nu exactly (Pauli spectrum, searches/full_clifford_classes.py).
  (b) stabilizer extent xi is monotone under stabilizer operations (Sec. 1.1, 2.2) and multiplicative on products of
      states of at most three qubits (Lemma 2.7(1)). Hence m T states make |psi> only if
      m >= log xi(psi) / log xi(T); Table 1 gives this ratio as 2.96818 for CS and 3.63356 for CCZ (W3: 3.63356,
      W4: 4.99907, W5: 5.93637). For a product of T, CS, CCZ on disjoint qubits, xi multiplies (each factor has <= 3
      qubits), so lambda >= ceil(a + 2.96818 b + 3.63356 c). The asterisk of Table 1 ("not known to hold for arbitrary
      catalysts") does not apply: we exclude catalysts, and the inputs T^m are products of 1-qubit states. For W4, W5 the
      extent itself is "computed to eight digits and reconstructed", not proved (App. A.4); those rows are flagged.
  (c) upper bounds without catalysts: tau(U) (exact T-count, staged_best 'tau'), minimized over the full-Clifford class
      (outputs/checks/full_clifford_classes.json): Clifford gates are free on this route. 4 T -> CCZ (Jones, [32] of
      Table 1) is tau of CS01,CS02, a member of the CCZ class, so it is covered by this minimum. For products of
      T, CS, CCZ on disjoint qubits also a + 3b + 4c (one T, three T for CS, four for CCZ, each without catalyst).
lambda is called known for a class when the best lower bound equals the best upper bound. Route cost at distance 2: the
fewest inputs for lambda T states from the cheapest distance-2 T^m factories of the catalogue
(outputs/large_unified_factories.csv, rows d=2 whose output is a product of T gates), batched; for products, also any
mix of directly distilled CS and CCZ states (pure CS^m, CCZ^m factories of the same catalogue) and T states for the
rest. To first order the error of the route is the sum over the distilled states, so it is a distance-2 route.
T supply, strongest distance-2 alternative we know besides the catalogue: Meier, Eastin, Knill, arXiv:1204.4221 (full
text read, Sec. II-III): the four-qubit-code routine takes ten |H> states to
two, accepts with a(p) = 1 - 10p + O(p^2), output error u(p) = 9p^2 + O(p^3) per output (every single error detected),
no catalyst. |H> = cos(pi/8)|0> + sin(pi/8)|1> is Clifford-equivalent to |T> = T|+> (checked numerically below), so it
is a T^2 factory with n = 10 for this route. Every T count below is the cheapest batching of the catalogue T^m factories
and MEK (10 inputs per 2 T states, surplus discarded). The CS and CCZ pieces use the cheaper of the catalogue pure
factories and T states (3 T per CS, 4 T per CCZ) from that supply.
Added 2026-11-06 (referee round 4), from Balada Gaggioli and Marecek, arXiv:2610.31115 (full text
read): (i) Prop. 7 and Table I, a nine-input distance-2 factory for one T
state (acceptance 1/128, error 27/2 p^2), and Prop. 5, no stabilizer protocol does it with eight, so T^1 costs 9 here
per run; per accepted output it costs more (that paper, Discussion: "per accepted output, the nine-input T factory
consumes more inputs than the ten-input protocol"), while the other factories used accept with 1 - n p + O(p^2). The
tallies are therefore also reported without it (tally_without_BGM, tally_full_clifford_without_BGM);
(ii) Prop. 10 (= Beverland et al. Prop. A.4): projecting one qubit of |CCZ> onto (|0> +- i|1>)/sqrt2 leaves CS|++> up
to a Clifford, so m CS states cost at most what m CCZ states cost (8 for one, [[8,3,2]]; Table I: n_min(CS,2) = 8).
Both are measurement routes without catalysts and are used for the T and CS pieces below (checked: cs_from_ccz_check);
the CS and CCZ pieces may also be drawn together from one batch of CCZ states (CCZ^m factories of the catalogue).
Full-Clifford groups are read from outputs/checks/full_clifford_classes.json; member labels are 'k:orbit:spec' and the
k=6 orbit labels contain colons ('6:0:1536:CCZ045,CCZ123'), so the orbit is everything between the first and last colon.
Invariants reported for the full-Clifford classes (all exact Clifford invariants, see invariants()): Pauli spectrum
(multiset of |<P>|), the number of Lagrangians inside the zero set, and the multiset of per-label Lagrangian counts.
Also recounted for every class: the nullity floor LB(nu) of outputs/nullity_l3.json (cheapest T/CS/CCZ supply of total
nullity >= nu; a lower bound for every such route, measurements included).

    python searches/lambda_route.py      -> outputs/checks/lambda_route.{json,md}
```

## Result

- Known answers reproduced: {'lambda(CCZ)=4': True, 'lambda(CS)=3': True, 'lambda(T)=1': True}.
- lambda is known (lower bound = upper bound) for **32 of 766** classes; by k: {1: 1, 2: 2, 3: 4, 4: 6, 5: 8, 6: 11}; lower bound that closes it: {'extent': 24, 'extent of W4 (Table 1; extent value numerical': 1, 'extent of W5 (Table 1; extent value numerical': 1, 'nullity': 6}.
- Gap upper - lower over all 766 classes (gap: classes): {0: 32, 1: 4, 2: 28, 3: 68, 4: 116, 5: 168, 6: 181, 7: 142, 8: 23, 9: 4}.
- Distance-2 factories with pure outputs used for the route: {'T': {1: 9, 2: 10, 3: 20, 4: 20, 5: 26, 6: 26}, 'CS': {1: 8, 2: 14, 3: 20, 4: 36}, 'CCZ': {1: 8, 2: 14, 3: 20}, 'T_catalogue': {1: 14, 2: 14, 3: 20, 4: 20, 5: 26, 6: 26}, 'T_MEK': {2: 10}, 'T_BGM': {1: 9}, 'CS_catalogue': {1: 12, 2: 18, 3: 24, 4: 36}}.

- Extent floor over any supply of distilled T, CS, CCZ states (any conversion by stabilizer operations, no catalysts; T supply includes MEK): our factory is below it for 0 of 32 of these classes, so this monotone proves no factory optimal against the measurement route; it only bounds the route from below.
- Rows where the route beats the factory (current table): T0 (route 9, factory 14; full-Clifford group ['1:1:T0']), T0,T1 (route 10, factory 14; full-Clifford group ['2:3:T0,T1']), CS01 (route 8, factory 12; full-Clifford group ['2:4:CS01']), T2,CS01 (route 17, factory 20; full-Clifford group ['3:12:T2,CS01']), CS01,CS02 (route 8, factory 12; full-Clifford group ['3:24:CS01,CS02', '3:64:CCZ012']), T0,T1,T2 (route 19, factory 20; full-Clifford group ['3:7:T0,T1,T2']), T3,CCZ012 (route 17, factory 20; full-Clifford group ['4:1032:T3,CCZ012', '4:56:T3,CS01,CS02']), CS03,CS12 (route 14, factory 18; full-Clifford group ['4:192:CS03,CS12']), T2,T3,CS01 (route 18, factory 22; full-Clifford group ['4:28:T2,T3,CS01']), T3,CS01,CS02 (route 17, factory 22; full-Clifford group ['4:1032:T3,CCZ012', '4:56:T3,CS01,CS02']), T3,T4,CS01,CS02 (route 18, factory 22; full-Clifford group ['5:120:T3,T4,CS01,CS02', '5:32792:T3,T4,CCZ012']), T3,T4,CCZ012 (route 18, factory 20; full-Clifford group ['5:120:T3,T4,CS01,CS02', '5:32792:T3,T4,CCZ012']), CS34,CCZ012 (route 14, factory 18; full-Clifford group ['5:49152:CS34,CCZ012', '5:896:CS03,CS04,CS12']), T4,CS03,CS12 (route 23, factory 26; full-Clifford group ['5:656:T4,CS03,CS12']), CS03,CS04,CS12 (route 14, factory 20; full-Clifford group ['5:49152:CS34,CCZ012', '5:896:CS03,CS04,CS12']), CS03,CS12,CCZ035,CCZ124 (route 14, factory 20; full-Clifford group ['6:0:1536:CCZ045,CCZ123', '6:2304:2304:CS03,CS12,CCZ035,CCZ124', '6:64:65544:CS01,CCZ015,CCZ234']), T4,CS03,CS12,CCZ035,CCZ124 (route 27, factory 28; full-Clifford group ['6:2320:2304:T4,CS03,CS12,CCZ035,CCZ124', '6:68:524289:T2,CS01,CCZ012,CCZ345']), T4,CS03,CS12,CCZ035 (route 23, factory 28; full-Clifford group ['6:2320:256:T4,CS03,CS12,CCZ035', '6:68:524288:T2,CS01,CCZ345']), T4,T5,CS03,CS12 (route 24, factory 28; full-Clifford group ['6:2352:0:T4,T5,CS03,CS12']), CS05,CS14,CS23 (route 20, factory 24; full-Clifford group ['6:41984:0:CS05,CS14,CS23']), CS01,CCZ015,CCZ234 (route 14, factory 18; full-Clifford group ['6:0:1536:CCZ045,CCZ123', '6:2304:2304:CS03,CS12,CCZ035,CCZ124', '6:64:65544:CS01,CCZ015,CCZ234']), T2,CS01,CCZ345 (route 23, factory 26; full-Clifford group ['6:2320:256:T4,CS03,CS12,CCZ035', '6:68:524288:T2,CS01,CCZ345']).
- Rows whose route cost the supplies outside the catalogue (MEK 10 -> 2 T, BGM 9 -> 1 T, CCZ + Y -> CS) lower: 21 of 32: ['T0', 'T0,T1', 'CS01', 'T2,CS01', 'T0,T1,T2', 'T3,CCZ012', 'CS03,CS12', 'T2,T3,CS01', 'T3,CS01,CS02', 'T3,T4,CS01,CS02', 'T3,T4,CCZ012', 'CS34,CCZ012', 'T4,CS03,CS12', 'CS03,CS04,CS12', 'T4,CS03,CS12,CCZ035,CCZ124', 'T4,CS03,CS12,CCZ035', 'T4,T5,CS03,CS12', 'T4,T5,CS03,CS12,CCZ034', 'CS05,CS14,CS23', 'T2,CS01,CCZ345', 'T2,CS01,CCZ012,CCZ345'].
- Lower bound that closes lambda: {'extent': 24, 'extent of W4 (Table 1; extent value numerical': 1, 'extent of W5 (Table 1; extent value numerical': 1, 'nullity': 6} (W4, W5 extents are numerical, App. A.4 of Beverland et al.).
- Per full-Clifford class (lambda known for 24 of 622; route and factory each minimised over the members): {'paper_v4': {'factory fewer': 3, 'route fewer': 14, 'same': 7}, 'current': {'factory fewer': 3, 'route fewer': 14, 'same': 7}}.

## Clifford invariants of the full-Clifford classes

- spectrum: 374 distinct values over 622 full-Clifford classes; groups with unequal values: 0.
- spectrum + n_lagrangians: 527 distinct values over 622 full-Clifford classes; groups with unequal values: 0.
- spectrum + n_lagrangians + multiset of per-label counts: 541 distinct values over 622 full-Clifford classes; groups with unequal values: 0.
- each count of distinct values is a lower bound on the number of full-Clifford classes; the exact count (union of label sets, full_clifford_classes.json) is the number of groups.

- Classes with lambda known, factory (paper_v4) vs cheapest lambda route: {'factory fewer': 3, 'route fewer': 22, 'same': 7} of 32.
- Same, without the nine-input T factory of arXiv:2610.31115 (acceptance 1/128): {'factory fewer': 3, 'route fewer': 20, 'same': 9} of 32; per full-Clifford class: {'factory fewer': 3, 'route fewer': 13, 'same': 8} of 24.
- Route vs the cheapest synthillation construction over each full-Clifford class (paper_v4): {'with_BGM': {'route fewer than synthillation': 14, 'same': 7, 'synthillation fewer': 3}, 'without_BGM': {'route fewer than synthillation': 13, 'same': 8, 'synthillation fewer': 3}} of 24.
- Classes with lambda known, factory (current) vs cheapest lambda route: {'factory fewer': 3, 'route fewer': 22, 'same': 7} of 32.
- Same, without the nine-input T factory of arXiv:2610.31115 (acceptance 1/128): {'factory fewer': 3, 'route fewer': 20, 'same': 9} of 32; per full-Clifford class: {'factory fewer': 3, 'route fewer': 13, 'same': 8} of 24.
- Route vs the cheapest synthillation construction over each full-Clifford class (current): {'with_BGM': {'route fewer than synthillation': 14, 'same': 7, 'synthillation fewer': 3}, 'without_BGM': {'route fewer than synthillation': 13, 'same': 8, 'synthillation fewer': 3}} of 24.
- All 766 classes, factory (paper_v4) vs nullity floor of every T/CS/CCZ route: {'above floor (open)': 763, 'at floor': 2, 'below floor': 1}.
- All 766 classes, factory (current) vs nullity floor of every T/CS/CCZ route: {'above floor (open)': 763, 'at floor': 2, 'below floor': 1}.

## Classes with lambda known

| k | spec | nu | tau | lambda | lower bound from | lambda route n (catalogue only) | extent floor, any T/CS/CCZ supply | factory n (paper) | factory n (current) | factory vs lambda route (paper) | (current) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | T0 | 1 | 1 | 1 | nullity | 9 (14) | 8 | 14 | 14 | route fewer | route fewer |
| 2 | T0,T1 | 2 | 2 | 2 | nullity | 10 (14) | 8 | 14 | 14 | route fewer | route fewer |
| 2 | CS01 | 2 | 3 | 3 | extent, product T^0 CS^1 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 8 (12) | 8 | 12 | 12 | route fewer | route fewer |
| 3 | T2,CS01 | 3 | 4 | 4 | extent, product T^1 CS^1 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 17 (20) | 14 | 20 | 20 | route fewer | route fewer |
| 3 | CS01,CS02 | 3 | 4 | 4 | extent, product T^0 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 8 (8) | 8 | 12 | 12 | route fewer | route fewer |
| 3 | CCZ012 | 3 | 7 | 4 | extent, product T^0 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 8 (8) | 8 | 8 | 8 | same | same |
| 3 | T0,T1,T2 | 3 | 3 | 3 | nullity | 19 (20) | 8 | 20 | 20 | route fewer | route fewer |
| 4 | T3,CCZ012 | 4 | 7 | 5 | extent, product T^1 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 17 (22) | 14 | 20 | 20 | route fewer | route fewer |
| 4 | T0,T1,T2,T3 | 4 | 4 | 4 | nullity | 20 (20) | 14 | 20 | 20 | same | same |
| 4 | CS03,CS12 | 4 | 6 | 6 | extent, product T^0 CS^2 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 14 (18) | 14 | 18 | 18 | route fewer | route fewer |
| 4 | CS01,CS03,CS12 | 4 | 5 | 5 | extent of W4 (Table 1; extent value numerical, App. A.4) | 26 (26) | 14 | 18 | 18 | factory fewer | factory fewer |
| 4 | T2,T3,CS01 | 4 | 5 | 5 | extent, product T^2 CS^1 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 18 (26) | 14 | 22 | 22 | route fewer | route fewer |
| 4 | T3,CS01,CS02 | 4 | 5 | 5 | extent, product T^1 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 17 (22) | 14 | 22 | 22 | route fewer | route fewer |
| 5 | T3,T4,CS01,CS02 | 5 | 6 | 6 | extent, product T^2 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 18 (22) | 14 | 22 | 22 | route fewer | route fewer |
| 5 | CS02,CS04,CS12,CS13 | 5 | 6 | 6 | extent of W5 (Table 1; extent value numerical, App. A.4) | 26 (26) | 14 | 18 | 18 | factory fewer | factory fewer |
| 5 | T0,T1,T2,T3,T4 | 5 | 5 | 5 | nullity | 26 (26) | 14 | 26 | 26 | same | same |
| 5 | T3,T4,CCZ012 | 5 | 8 | 6 | extent, product T^2 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 18 (22) | 14 | 20 | 20 | route fewer | route fewer |
| 5 | CS34,CCZ012 | 5 | 9 | 7 | extent, product T^0 CS^1 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 14 (20) | 14 | 18 | 18 | route fewer | route fewer |
| 5 | T2,T3,T4,CS01 | 5 | 6 | 6 | extent, product T^3 CS^1 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 26 (26) | 14 | 26 | 26 | same | same |
| 5 | T4,CS03,CS12 | 5 | 7 | 7 | extent, product T^1 CS^2 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 23 (32) | 14 | 26 | 26 | route fewer | route fewer |
| 5 | CS03,CS04,CS12 | 5 | 7 | 7 | extent, product T^0 CS^1 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 14 (20) | 14 | 20 | 20 | route fewer | route fewer |
| 6 | CCZ045,CCZ123 | 6 | 13 | 8 | extent, product T^0 CS^0 CCZ^2 (Lemma 2.7(1), Table 1 ratios) | 14 (14) | 14 | 14 | 14 | same | same |
| 6 | CS03,CS12,CCZ035,CCZ124 | 6 | 8 | 8 | extent, product T^0 CS^0 CCZ^2 (Lemma 2.7(1), Table 1 ratios) | 14 (14) | 14 | 20 | 20 | route fewer | route fewer |
| 6 | T4,CS03,CS12,CCZ035,CCZ124 | 6 | 7 | 7 | extent, product T^3 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 27 (28) | 14 | 28 | 28 | route fewer | route fewer |
| 6 | T4,CS03,CS12,CCZ035 | 6 | 8 | 8 | extent, product T^1 CS^1 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 23 (28) | 20 | 28 | 28 | route fewer | route fewer |
| 6 | T4,T5,CS03,CS12 | 6 | 8 | 8 | extent, product T^2 CS^2 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 24 (32) | 20 | 28 | 28 | route fewer | route fewer |
| 6 | T4,T5,CS03,CS12,CCZ034 | 6 | 7 | 7 | extent, product T^4 CS^1 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 28 (32) | 14 | 28 | 28 | same | same |
| 6 | T4,T5,CS03,CS12,CCZ035,CCZ124 | 6 | 6 | 6 | nullity | 26 (26) | 14 | 26 | 26 | same | same |
| 6 | CS05,CS14,CS23 | 6 | 9 | 9 | extent, product T^0 CS^3 CCZ^0 (Lemma 2.7(1), Table 1 ratios) | 20 (24) | 20 | 24 | 24 | route fewer | route fewer |
| 6 | CS01,CCZ015,CCZ234 | 6 | 10 | 8 | extent, product T^0 CS^0 CCZ^2 (Lemma 2.7(1), Table 1 ratios) | 14 (14) | 14 | 18 | 18 | route fewer | route fewer |
| 6 | T2,CS01,CCZ345 | 6 | 10 | 8 | extent, product T^1 CS^1 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 23 (28) | 20 | 26 | 26 | route fewer | route fewer |
| 6 | T2,CS01,CCZ012,CCZ345 | 6 | 9 | 7 | extent, product T^3 CS^0 CCZ^1 (Lemma 2.7(1), Table 1 ratios) | 27 (28) | 14 | 26 | 26 | factory fewer | factory fewer |
