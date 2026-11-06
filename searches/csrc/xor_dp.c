/* Min-cost XOR subset DP over a syndrome space of 2^bits (bits <= 34), uint8 saturating costs (255 = unreachable),
   one byte per state. Input (stdin, text): bits h m, then m lines "col cost" (col < 2^bits, cost 1..254).
   dp[s] = min total cost of a subset of columns XORing to s. Output (stdout): for every final state f = t << h
   (t < 2^(bits-h)) with dp[f] < 255, one line "t n c1 c2 ..." with the column indices of one optimum.
   Traceback uses the final table only: from s with value v pick any column c with dp[s^col_c] + cost_c == v; in an
   optimum no column repeats (two equal columns cancel at positive cost), so the chain ends at 0 with total v.
   Build: cc -O3 -march=native -o xor_dp xor_dp.c */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(void) {
    int bits, h, m;
    if (scanf("%d %d %d", &bits, &h, &m) != 3) return 1;
    uint64_t *col = malloc(sizeof(uint64_t) * m); int *cost = malloc(sizeof(int) * m);
    for (int i = 0; i < m; i++) if (scanf("%llu %d", (unsigned long long *)&col[i], &cost[i]) != 2) return 1;
    uint64_t S = (uint64_t)1 << bits;
    uint8_t *dp = malloc(S);
    if (!dp) { fprintf(stderr, "alloc failed\n"); return 2; }
    for (uint64_t s = 0; s < S; s++) dp[s] = 255;
    dp[0] = 0;
    for (int i = 0; i < m; i++) {
        uint64_t c = col[i]; if (!c) continue;
        int hb = 63 - __builtin_clzll(c); uint64_t B = (uint64_t)1 << hb; uint8_t w = (uint8_t)cost[i];
        for (uint64_t hi = 0; hi < S; hi += 2 * B)
            for (uint64_t lo = 0; lo < B; lo++) {
                uint64_t s = hi | lo, t = s ^ c;
                uint8_t a = dp[s], b = dp[t];
                unsigned aw = a + w, bw = b + w;
                if (aw > 255) aw = 255;
                if (bw > 255) bw = 255;
                dp[s] = (bw < a) ? (uint8_t)bw : a;
                dp[t] = (aw < b) ? (uint8_t)aw : b;
            }
    }
    int *path = malloc(sizeof(int) * (m + 1));
    uint64_t T = (uint64_t)1 << (bits - h);
    for (uint64_t t = 0; t < T; t++) {
        uint64_t s = t << h; int v = dp[s];
        if (v == 255) continue;
        int np = 0; uint64_t x = s; int left = v;
        while (x) {
            int found = -1;
            for (int i = 0; i < m; i++) if (col[i] && dp[x ^ col[i]] + cost[i] == left) { found = i; break; }
            if (found < 0) { fprintf(stderr, "traceback failed at %llu\n", (unsigned long long)s); return 3; }
            path[np++] = found; left -= cost[found]; x ^= col[found];
        }
        if (left) { fprintf(stderr, "traceback cost mismatch\n"); return 3; }
        printf("%llu %d", (unsigned long long)t, v);
        for (int j = 0; j < np; j++) printf(" %d", path[j]);
        printf("\n");
    }
    free(dp); return 0;
}
