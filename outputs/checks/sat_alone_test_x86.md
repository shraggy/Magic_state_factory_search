# Test A on x86-64: the three arm64 30-input factories

Same code (md5-identical chain_iii2/3/4, xor_sat2, staged, staged_class, sat_alone_test apart from the output folder), same class cache, seed and caps, same layer chain and 2e6 conflicts per query, one fresh process per query. LSF job 1230451 on ccc (x86_64, z3 5.3.2, CryptoMiniSat 5.17.0), 1 core per task.

- 2304:136 seed: x86 [24, 24, 22, 22, 22, 22, 22] (final 22); arm64 [22, 22, 22, 22, 22, 22, 22] (final 22); v6 [24, 24, 22, 22, 22, 22, 22]; pipeline n 22, synthillation 24
- 2320:9041 seed: x86 [30, 30, 30, 30, 30, 28, 28] (final 28); arm64 [30, 30, 30, 30, 30, 30, 30] (final 30); v6 [30, 30, 30, 30, 30, 28, 28]; pipeline n 28, synthillation 32
- 2352:1674 none: x86 [32, 32, 32, 32, 32, 32, 32] (final 32); arm64 [32, 32, 32, 32, 32, 30, 30] (final 30); v6 [32, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 32
- 2352:1674 seed: x86 [32, 32, 32, 32, 32, 32, 32] (final 32); arm64 [30, 30, 30, 30, 30, 30, 30] (final 30); v6 [32, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 32
- 2352:3240 none: x86 [34, 32, 32, 32, 32, 32, 32] (final 32); arm64 [32, 32, 32, 32, 32, 30, 30] (final 30); v6 [34, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 34
- 2352:3240 seed: x86 [34, 32, 32, 32, 32, 32, 32] (final 32); arm64 [32, 32, 32, 32, 32, 30, 30] (final 30); v6 [34, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 34
- 2352:3240 syn: x86 [34, 32, 32, 32, 32, 32, 32] (final 32); arm64 [32, 32, 32, 32, 32, 30, 30] (final 30); v6 [34, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 34
- 2352:8944 none: x86 [32, 32, 32, 32, 32, 32, 32] (final 32); arm64 [34, 34, 30, 30, 30, 30, 30] (final 30); v6 [34, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 34
- 2352:8944 seed: x86 [34, 32, 32, 32, 32, 32, 32] (final 32); arm64 [34, 34, 30, 30, 30, 30, 30] (final 30); v6 [34, 32, 32, 32, 32, 32, 32]; pipeline n 32, synthillation 34

Summary: {'known_answer_match': '2 of 2', 'three_classes_x86_final_32': '7 of 7', 'three_classes_arm64_final_30': '7 of 7', 'seeded_rows_match_v6': '5 of 5'}

Conclusion: on x86-64 the known answers reproduce v6 layer by layer and every run of the three classes ends at the pipeline's 32, so the 30-input factories found on arm64 are an effect of the arm64 build, not of running SAT without stages (i)-(ii). They stay valid, verified factories (outputs/checks/sat_alone_test/verify.json) kept outside the pipeline.
