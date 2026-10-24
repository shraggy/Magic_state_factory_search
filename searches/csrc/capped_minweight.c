// cap.c with 256-bit vectors (4 x uint64).
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>
typedef struct { uint64_t w[4]; } V;
static int d, wmax, cap, nthr; static V b, G[400];
static volatile int found = 0, best = 1 << 30; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static inline int pc(const V *v){ return __builtin_popcountll(v->w[0]) + __builtin_popcountll(v->w[1]) + __builtin_popcountll(v->w[2]) + __builtin_popcountll(v->w[3]); }
static inline V x(V a, const V *g){ for (int i = 0; i < 4; i++) a.w[i] ^= g->w[i]; return a; }
static void rec(V v, int start, int depth, int *lb){
    int c = pc(&v); if (c < *lb) *lb = c;
    if (c < cap) found = 1;
    if (depth == wmax || found) return;
    for (int i = start; i < d; i++) rec(x(v, &G[i]), i + 1, depth + 1, lb);
}
static void *work(void *arg){
    long t = (long)arg; int lb = 1 << 30;
    if (t == 0){ int c = pc(&b); if (c < lb) lb = c; if (c < cap) found = 1; }
    for (int i = t; i < d && !found; i += nthr) rec(x(b, &G[i]), i + 1, 1, &lb);
    pthread_mutex_lock(&mu); if (lb < best) best = lb; pthread_mutex_unlock(&mu); return NULL;
}
static V rd(void){ V v; for (int i = 0; i < 4; i++){ unsigned long long a; if (scanf("%llx", &a) != 1) exit(2); v.w[i] = a; } return v; }
int main(int argc, char **argv){
    nthr = atoi(argv[1]); cap = atoi(argv[2]); wmax = atoi(argv[3]);
    if (scanf("%d", &d) != 1) return 2; b = rd(); for (int i = 0; i < d; i++) G[i] = rd();
    pthread_t th[64]; for (long t = 0; t < nthr; t++) pthread_create(&th[t], NULL, work, (void*)t);
    for (int t = 0; t < nthr; t++) pthread_join(th[t], NULL);
    printf("%d %d\n", found, best); return 0;
}
