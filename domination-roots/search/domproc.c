/* geng OUTPROC: count dominating sets by size for every generated graph and
   evaluate D(G,x) at x = -2..-10 (and the odd values as a sanity check).
   Independent of the repo: plain subset enumeration from closed neighbourhoods.
   Totals printed at exit: graphs seen, zeros per root, any zero graphs in graph6. */
#include "gtools.h"
#include <stdint.h>
#include <stdlib.h>

static long long ngraphs = 0;
static long long zeros[11];
static int registered = 0;
static int nseen = 0;

static void report(void)
{
    fprintf(stderr, "RESULT n=%d graphs=%lld", nseen, ngraphs);
    for (int r = 2; r <= 10; ++r) fprintf(stderr, " z%d=%lld", r, zeros[r]);
    fprintf(stderr, "\n");
}

void domproc(FILE *f, graph *g, int n)
{
    static uint32_t cover[1 << 14];
    uint32_t closed[32];
    long long cnt[33];
    if (!registered) { atexit(report); registered = 1; }
    nseen = n;
    ++ngraphs;
    for (int i = 0; i < n; ++i) {
        uint32_t m = 1u << i;
        set *gi = GRAPHROW(g, i, 1);
        for (int j = 0; j < n; ++j)
            if (ISELEMENT(gi, j)) m |= 1u << j;
        closed[i] = m;
    }
    for (int k = 0; k <= n; ++k) cnt[k] = 0;
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    cover[0] = 0;
    for (uint32_t s = 1; s <= full; ++s) {
        uint32_t low = s & (~s + 1);
        cover[s] = cover[s ^ low] | closed[__builtin_ctz(low)];
        if (cover[s] == full) cnt[__builtin_popcount(s)]++;
    }
    for (int r = 2; r <= 10; ++r) {
        __int128 v = 0, p = 1;
        for (int k = 0; k <= n; ++k) { v += (__int128)cnt[k] * p; p *= -r; }
        if (v == 0) {
            zeros[r]++;
            if (r >= 3) {
                fprintf(stderr, "ZERO r=-%d ", r);
                writeg6(stderr, g, 1, n);
            }
        }
    }
}
