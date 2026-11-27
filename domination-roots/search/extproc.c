/* Exhaustive one-vertex extension test at x = -ROOT.
   Every connected graph of order n+1 is H + v with H connected of order n
   (delete a leaf of a spanning tree) and N(v) = T, a nonempty subset of V(H).
   D(H+v, x) = D(H,x) - U_H(V\T) + x * Q_{H,T}(x), where
     U_H(B)   = sum over dominating sets S of H with S subset of B of x^|S|,
     Q_{H,T}  = sum over S dominating V(H)\T of x^|S|.
   (Formula checked against brute force by verify.py, function search_small.)
   Used as a geng OUTPROC (-DOUTPROC=extproc) or via the standalone driver. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef ROOT
#define ROOT 4
#endif

long long ext_bases = 0, ext_tests = 0, ext_zeros = 0, ext_tpl = 0;
static void (*tplcb)(uint32_t T, int k, void *ctx) = 0;

/* returns number of T with value zero; prints them via cb */
int ext_check(int n, const uint32_t *closed, void (*cb)(uint32_t T, void *), void *ctx)
{
    static uint32_t cover[1 << 15];
    static int64_t unc[1 << 15], dom[1 << 15], w[16];
    uint32_t full = (1u << n) - 1, N = 1u << n;
    w[0] = 1;
    for (int k = 1; k <= n; ++k) w[k] = w[k - 1] * (-ROOT);
    for (uint32_t s = 0; s < N; ++s) unc[s] = dom[s] = 0;
    cover[0] = 0;
    unc[full] += 1;
    for (uint32_t s = 1; s < N; ++s) {
        uint32_t low = s & (~s + 1);
        cover[s] = cover[s ^ low] | closed[__builtin_ctz(low)];
        int64_t ws = w[__builtin_popcount(s)];
        unc[full ^ cover[s]] += ws;
        if (cover[s] == full) dom[s] = ws;
    }
    int64_t P = unc[0];
    for (uint32_t bit = 1; bit < N; bit <<= 1)
        for (uint32_t base = 0; base < N; base += 2 * bit)
            for (uint32_t j = base; j < base + bit; ++j) {
                unc[j + bit] += unc[j];
                dom[j + bit] += dom[j];
            }
    int z = 0;
    for (uint32_t T = 1; T < N; ++T) {
        int64_t v = P - dom[full ^ T] - (int64_t)ROOT * unc[T];
        if (v == 0) { ++z; if (cb) cb(T, ctx); }
#if ROOT == 4
        /* hub template: H plus hub h adjacent to T and to k private leaves,
           order n+1+k; D(-4)=0 iff 4^(k-1) P = 3^k Q_T (Q_T = unc[T]).
           Report every k with n+1+k <= 15 (k=5 at n=9 recovers G15). */
        {
            __int128 a = P, b = 3 * (__int128)unc[T];
            for (int k = 1; n + 1 + k <= 15; ++k) {
                if (a == b && tplcb) tplcb(T, k, ctx);
                a *= 4; b *= 3;
            }
        }
#endif
    }
    ext_bases++;
    ext_tests += N - 1;
    ext_zeros += z;
    return z;
}

#ifdef OUTPROC
#include "gtools.h"
static int registered = 0;
static int nseen = 0;
static void report(void)
{
    fprintf(stderr, "RESULT base_n=%d bases=%lld tests=%lld zeros=%lld tpl=%lld root=-%d\n",
            nseen, ext_bases, ext_tests, ext_zeros, ext_tpl, ROOT);
}
static graph *curg;
static void print_zero(uint32_t T, void *ctx)
{
    fprintf(stderr, "ZERO T=%u base=", T);
    writeg6(stderr, curg, 1, nseen);
}
static void print_tpl(uint32_t T, int k, void *ctx)
{
    ext_tpl++;
    fprintf(stderr, "TPL k=%d order=%d T=%u base=", k, nseen + 1 + k, T);
    writeg6(stderr, curg, 1, nseen);
}
void extproc(FILE *f, graph *g, int n)
{
    tplcb = print_tpl;
    uint32_t closed[32];
    if (!registered) { atexit(report); registered = 1; }
    nseen = n;
    curg = g;
    for (int i = 0; i < n; ++i) {
        uint32_t m = 1u << i;
        set *gi = GRAPHROW(g, i, 1);
        for (int j = 0; j < n; ++j)
            if (ISELEMENT(gi, j)) m |= 1u << j;
        closed[i] = m;
    }
    ext_check(n, closed, print_zero, NULL);
}
#else
/* standalone driver: reads "n m u1 v1 u2 v2 ..." from stdin */
static void pz(uint32_t T, void *ctx) { printf("zero at T=%u\n", T); }
int main(void)
{
    int n, m;
    uint32_t closed[32];
    if (scanf("%d %d", &n, &m) != 2) return 1;
    for (int i = 0; i < n; ++i) closed[i] = 1u << i;
    for (int e = 0; e < m; ++e) {
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 1;
        closed[u] |= 1u << v; closed[v] |= 1u << u;
    }
    int z = ext_check(n, closed, pz, NULL);
    printf("zeros=%d\n", z);
    return 0;
}
#endif
