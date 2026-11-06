/* Enumerate every Lagrangian subspace L of F_2^{2k} (k <= 6) with L \ {0} inside a given set Z.
   A vector u = (a << k) | b stands for the Pauli X^a Z^b. Each L is printed once, as its reduced basis
   (pivots = leading bits strictly decreasing, every basis vector zero on the other pivots).
   stdin: k, then 4^k characters '0'/'1' (Z membership of u = 0 .. 4^k - 1). stdout: one line per L, k integers.
   Used by searches/full_clifford_classes.py. Build: cc -O2 -o lagrangians lagrangians.c */
#include <stdio.h>
#include <stdlib.h>
static int K, M, NB; static unsigned char Z[4096];
static int sympl(int u, int v) { int m = M; return (__builtin_popcount(((u >> K) & v & m)) + __builtin_popcount((u & m) & (v >> K))) & 1; }
static int basis[6], piv[6], span[64];
static long nfound = 0;
static void rec(int d, int *cand, int nc) {
    if (d == K) { for (int i = 0; i < K; i++) printf(i ? " %d" : "%d", basis[i]); printf("\n"); nfound++; return; }
    int lo = d ? piv[d - 1] : 2 * K;
    int *next = malloc(sizeof(int) * (nc + 1));
    int ns = 1 << d;
    for (int ci = 0; ci < nc; ci++) {
        int v = cand[ci]; int lb = 31 - __builtin_clz(v);
        if (lb >= lo) continue;
        int bad = 0;
        for (int p = 0; p < d; p++) if ((v >> piv[p]) & 1) { bad = 1; break; }
        if (bad) continue;
        for (int p = 0; p < d; p++) if ((basis[p] >> lb) & 1) { bad = 1; break; }
        if (bad) continue;
        basis[d] = v; piv[d] = lb;
        for (int s = 0; s < ns; s++) span[ns + s] = span[s] ^ v;
        int nn = 0;
        for (int cj = 0; cj < nc; cj++) {
            int w = cand[cj]; if (sympl(v, w)) continue;
            int ok = 1;
            for (int s = ns; s < 2 * ns; s++) if (!Z[w ^ span[s]]) { ok = 0; break; }
            if (ok) next[nn++] = w;
        }
        rec(d + 1, next, nn);
    }
    free(next);
}
int main(void) {
    if (scanf("%d", &K) != 1 || K < 1 || K > 6) return 1;
    M = (1 << K) - 1; NB = 1 << (2 * K);
    char *buf = malloc(NB + 2); if (scanf("%s", buf) != 1) return 1;
    int *cand = malloc(sizeof(int) * NB); int nc = 0;
    for (int u = 0; u < NB; u++) { Z[u] = buf[u] == '1'; if (u && Z[u]) cand[nc++] = u; }
    span[0] = 0;
    rec(0, cand, nc);
    fprintf(stderr, "lagrangians %ld\n", nfound);
    return 0;
}
