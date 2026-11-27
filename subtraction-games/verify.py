#!/usr/bin/env python3
"""Verifier for "Necessity of Manabe's pure-periodicity criterion for three-move
subtraction games with b <= 2a or b <= 250" (paper.tex in this directory).

Standard library only.  Every check recomputes from the definitions: P-positions
from the subtraction-game recurrence, nim-values from the mex rule, and the
admissible angle sets from Manabe's Definition 19.  No stored result is read.

    python3 -B verify.py          # default mode, well under a minute on 16 cores
    python3 -B verify.py --full   # all certificates for b <= 250, see logs/verify-full.log

Any mismatch raises AssertionError (non-zero exit).  Each claim prints one PASS line.
"""
import os
import sys
import time
from math import gcd
from multiprocessing import Pool

# Numbers stated in paper.tex.  The full run recomputes and asserts each of them;
# both modes check that the paper prints exactly these values.
CLAIMS = {
    "box_admissible": 7557,          # Manabe's box, admissible (pure) triples
    "box_nonadmissible": 7828,       # Manabe's box, non-admissible triples
    "classes_covered": 3902781,      # primitive classes (a,b,r), b <= 250, admissible
    "classes_witnessed": 2568108,    # primitive classes, non-admissible, witness at c0 = 2m+r
    "short_covered": 5579832,        # primitive non-additive b < c < 2m, b <= 250, admissible
    "short_witnessed": 3020764,      # same range, non-admissible, witness found
}
BMAX_FULL = 250
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- games

def outcome_bits(a, b, c, length):
    """Integer whose bit x is 1 iff x is a P-position of {a,b,c}, for 0 <= x < length.

    Y(x) = 1 iff no legal move x -> x-s (s in {a,b,c}, s <= x) reaches a P-position.
    Positions [j, j+a) depend only on positions below j, so they are set together."""
    bits, chunk = 0, (1 << a) - 1
    for j in range(0, length, a):
        options = (bits << a) | (bits << b) | (bits << c)
        bits |= (~(options >> j) & chunk) << j
    return bits & ((1 << length) - 1)


def nim_values(moves, length):
    """Sprague-Grundy values by the mex rule, positions 0..length-1."""
    g = bytearray()
    for n in range(length):
        options = {g[n - s] for s in moves if s <= n}
        v = 0
        while v in options:
            v += 1
        g.append(v)
    return g


def nim_cycle(a, b, c):
    """(q, p): pre-period and least period of G_S, from the first repeated state.

    The state at t is G[t, t+c); it determines G(t+c), so the first repeat t' of an
    earlier state t gives q = t and least period p = t' - t exactly."""
    g, seen, n = bytearray(), {}, 0
    while True:
        options = {g[n - s] for s in (a, b, c) if s <= n}
        v = 0
        while v in options:
            v += 1
        g.append(v)
        if n >= c - 1:
            t = n - c + 1
            state = bytes(g[t:])
            if state in seen:
                return seen[state], t - seen[state]
            seen[state] = t
        n += 1


def find_witness(bits, length, a, b, c, lo, hi):
    """Least (y, s) with lo <= y < hi, s in {a,b,c}, h = a+b+c-s and
    Y(y) = Y(y+s) = 0, Y(y+h) = 1 (a witness, Lemma 2.4), or None."""
    assert hi + a + b + c <= length
    window = (1 << hi) - (1 << lo)
    best = None
    for s in (a, b, c):
        h = a + b + c - s
        hits = ~bits & ~(bits >> s) & (bits >> h) & window
        if hits:
            y = (hits & -hits).bit_length() - 1
            if best is None or y < best[0]:
                best = (y, s)
    return best


def holds(bits, a, b, c, y, s):
    h = a + b + c - s
    return not bits >> y & 1 and not bits >> (y + s) & 1 and bits >> (y + h) & 1


# ------------------------------------------------- Manabe, Definition 19

def interval(lo, hi):
    return (1 << hi) - (1 << lo) if hi > lo else 0


def rotate(x, t, m):
    """The circle translate x + t of a subset x of Z/m (bitmask)."""
    t %= m
    return ((x << t) | (x >> (m - t))) & ((1 << m) - 1)


def admissible(a, b):
    """Bitmasks on Z/(a+b) of Theta_{a+b}, Theta_{c+a}, Theta_{c+b} (Definition 19)."""
    m, delta = a + b, b - a
    z = 0                                   # Z = W0({a,b}) on [0, m), by the recurrence
    for n in range(m):
        if not any(n >= s and z >> (n - s) & 1 for s in (a, b)):
            z |= 1 << n
    eta, eps = divmod(delta, a)             # formula (4)
    if eta % 2 == 0:
        lam = interval((eta + 1) * a, (eta + 1) * a + eps)
    else:
        lam = interval(eta * a + eps, (eta + 1) * a)
    gaps = interval(a, b) & ~z              # D = [a, b) \ Z
    low = z & interval(0, delta)            # Z ∩ [0, delta)
    t0 = t1 = t2 = 0
    for r in range(m):
        if z & rotate(z, r, m) == 0:
            t0 |= 1 << r
        if low & rotate(z, delta - r, m) == 0:
            t1 |= 1 << r
        if rotate(lam, r, m) & z == 0 and rotate(gaps, r, m) & lam == 0:
            t2 |= 1 << r
    return t0, t1, t2


# ------------------------------------------------------------- checks

def check_generators():
    """Bit-parallel P-positions agree with per-position mex zero sets."""
    n = 0
    for a in range(1, 30):
        for b in range(a + 1, 40):
            for c in range(b + 1, 41):
                length = 3 * c + 7
                bits = outcome_bits(a, b, c, length)
                g = nim_values((a, b, c), length)
                assert all((bits >> x & 1) == (g[x] == 0) for x in range(length)), (a, b, c)
                n += 1
    return f"bit-parallel recurrence = mex zero set on {n:,} triples with c <= 40"


def check_box():
    """Manabe's box: exact nim cycles versus the angle criterion, and Lemma 2.4."""
    assert admissible(3, 8) == (sum(1 << i for i in (3, 8)), sum(1 << i for i in (2, 8)),
                                sum(1 << i for i in (3, 4, 5, 9, 10)))   # Example 1
    pure = nonpure = 0
    for a in range(2, 9):
        for b in range(a + 1, 41):
            m = a + b
            tests = admissible(a, b)
            for c in range(b + 1, 101):
                if c == m or gcd(gcd(a, b), c) > 1:
                    continue
                flags = [t >> (c % m) & 1 for t in tests]
                q, p = nim_cycle(a, b, c)
                assert (q == 0) == any(flags), (a, b, c, q, p)
                length = 6 * c
                bits = outcome_bits(a, b, c, length)
                found = find_witness(bits, length, a, b, c, 0, length - a - b - c)
                if any(flags):
                    assert p == min(P for P, f in zip((m, c + a, c + b), flags) if f)
                    assert found is None, (a, b, c)   # Lemma 2.4: pure games have none
                    pure += 1
                else:
                    assert found is not None, (a, b, c)
                    nonpure += 1
    assert (pure, nonpure) == (CLAIMS["box_admissible"], CLAIMS["box_nonadmissible"])
    return (f"Manabe's box: {pure:,} admissible pure with least period = least admissible, "
            f"{nonpure:,} non-admissible not pure, each with a Lemma 2.4 witness; "
            "no witness in any pure game")


def check_theorem1(amax):
    """Formula (3.1) and the witness of Theorem 3.2, plus Example 3.3."""
    n = 0
    for a in range(2, amax + 1):
        for b in range(a + 2, 2 * a + 1):
            m, delta = a + b, b - a
            for r in range(1, delta):
                hset = interval(a, b - r) | interval(b, 2 * a) | interval(m - r, m)
                for k in (1, 2, 3):
                    c = k * m + r
                    length = 2 * c + m
                    bits = outcome_bits(a, b, c, length)
                    reps = c // m + 1
                    periodic = sum(hset << (i * m) for i in range(reps)) & interval(0, c)
                    assert (bits >> c) & interval(0, c) == periodic, (a, b, c)
                    y = k * m - a - 1
                    assert holds(bits, a, b, c, y, c), (a, b, c)
                    n += 1
    bits = outcome_bits(5, 8, 28, 100)
    assert [bits >> x & 1 for x in (20, 48, 33)] == [0, 0, 1]
    return (f"formula (3.1) for all 0 <= u < c and witness y = km-a-1 on {n:,} cases "
            f"(a <= {amax}, k = 1,2,3); example (5,8,28)")


def transfer_task(args):
    a, b = args
    m, n = a + b, 0
    for r in range(m):
        seen = None
        for k in (2, 3, 4, 5):
            c = k * m + r
            bits = outcome_bits(a, b, c, 2 * c + 2 * m)
            windows = ((bits >> (c - m)) & interval(0, 3 * m),
                       (bits >> (2 * c - m)) & interval(0, 3 * m))
            assert seen is None or windows == seen, (a, b, r, k)
            seen = windows
        n += 1
    return n


def check_transfer(pool, bmax):
    tasks = [(a, b) for b in range(2, bmax + 1) for a in range(1, b)]
    n = sum(pool.imap_unordered(transfer_task, tasks, chunksize=8))
    return (f"Y(c+j), Y(2c+j), -m <= j < 2m, identical for k = 2..5 on all {n:,} "
            f"classes (a,b,r) with 1 <= a < b <= {bmax}")


def pair_task(args):
    """All certificate work for one base (a, b): Theta shapes, residue classes, c < 2m."""
    a, b, recheck = args
    m, delta = a + b, b - a
    t0, t1, t2 = admissible(a, b)
    covered = t0 | t1 | t2
    if b <= 2 * a:                                   # Lemma 3.1
        assert (t0, t1, t2) == (interval(a, b + 1), 1 | interval(b, m), interval(delta, a + 1))
        assert covered == interval(0, m) & ~interval(1, delta)
    elif b % a == 0 and (b // a) % 2 == 1:           # proof of Proposition 5.1
        assert t2 == interval(0, m)
    else:
        assert t0 == (1 << a) | (1 << b)
    out = dict(classes_covered=0, classes_witnessed=0, rechecked=0,
               short_covered=0, short_witnessed=0)
    g = gcd(a, b)
    for r in range(m):                               # Theorem 4.2, c >= 2m
        if gcd(g, r) > 1:
            continue
        if covered >> r & 1:
            out["classes_covered"] += 1
            continue
        c0 = 2 * m + r
        length = 2 * c0 + 2 * m
        bits = outcome_bits(a, b, c0, length)
        found = find_witness(bits, length, a, b, c0, c0 - m, c0 + m)
        assert found is not None, ("no class witness", a, b, r)
        out["classes_witnessed"] += 1
        j, kind = found[0] - c0, (a, b, c0).index(found[1])
        for k in recheck:                            # direct check, without Lemma 4.1
            c = k * m + r
            bits_k = outcome_bits(a, b, c, 2 * c + 2 * m)
            assert holds(bits_k, a, b, c, c + j, (a, b, c)[kind]), ("recheck", a, b, r, k)
            out["rechecked"] += 1
    for c in range(b + 1, 2 * m):                    # Theorem 4.2, b < c < 2m
        if c == m or gcd(g, c) > 1:
            continue
        if covered >> (c % m) & 1:
            out["short_covered"] += 1
            continue
        length = 4 * c + 2 * m
        while True:
            bits = outcome_bits(a, b, c, length)
            found = find_witness(bits, length, a, b, c, 0, length - a - b - c)
            if found is not None:
                break
            assert length < 64 * c, ("no witness found", a, b, c)
            length *= 2
        out["short_witnessed"] += 1
    return out


def check_certificates(pool, bmax, recheck):
    tasks = [(a, b, recheck) for b in range(bmax, 2, -1) for a in range(2, b)]
    total = dict.fromkeys(("classes_covered", "classes_witnessed", "rechecked",
                           "short_covered", "short_witnessed"), 0)
    for out in pool.imap_unordered(pair_task, tasks, chunksize=4):
        for key, value in out.items():
            total[key] += value
    if bmax == BMAX_FULL:
        for key in ("classes_covered", "classes_witnessed", "short_covered", "short_witnessed"):
            assert total[key] == CLAIMS[key], (key, total[key], CLAIMS[key])
    return total


def check_paper_numbers():
    with open(os.path.join(HERE, "paper.tex"), encoding="utf-8") as f:
        text = f.read()
    for key, value in CLAIMS.items():
        assert f"{value:,}".replace(",", "{,}") in text, ("paper.tex lacks", key, value)
    return f"paper.tex states the {len(CLAIMS)} computed counts"


def main():
    if not __debug__:
        raise SystemExit("run without -O: the checks are assertions")
    full = sys.argv[1:] == ["--full"]
    assert full or not sys.argv[1:], "usage: verify.py [--full]"
    bmax = BMAX_FULL if full else 100
    start = time.time()

    def report(claim, message, t):
        print(f"PASS {claim}: {message} [{time.time() - t:.1f}s]", flush=True)

    t = time.time()
    report("definitions", check_generators(), t)
    t = time.time()
    report("consistency with Manabe Thm 34", check_box(), t)
    t = time.time()
    report("Theorem 3.2 (b <= 2a)", check_theorem1(40 if full else 25), t)
    with Pool(os.cpu_count()) as pool:
        t = time.time()
        report("Lemma 4.1 (transfer)", check_transfer(pool, 60 if full else 30), t)
        t = time.time()
        total = check_certificates(pool, bmax, (3,) if full else (3, 4))
    scope = "all" if full else "subset:"
    report("Lemma 3.1 and Prop 5.1 angle sets",
           f"Theta formulas for b <= 2a, Theta_(a+b) = {{a,b}} or Theta_(c+b) = Z/m "
           f"for b > 2a, every 2 <= a < b <= {bmax}", t)
    report("Theorem 4.2, c >= 2m",
           f"{scope} b <= {bmax}: {total['classes_covered']:,} admissible primitive classes, "
           f"{total['classes_witnessed']:,} non-admissible classes each with a witness at "
           f"c0 = 2m+r, -m <= j < m"
           + f"; {total['rechecked']:,} witnesses rechecked directly at k = "
           + ("3" if full else "3, 4"), t)
    report("Theorem 4.2, b < c < 2m",
           f"{scope} b <= {bmax}: {total['short_covered']:,} admissible triples, "
           f"{total['short_witnessed']:,} non-admissible triples each with a witness", t)
    t = time.time()
    report("paper numbers", check_paper_numbers(), t)
    print(f"ALL PASS ({'full' if full else 'default'} mode, {time.time() - start:.1f}s, "
          f"{os.cpu_count()} processes, Python {sys.version.split()[0]})")


if __name__ == "__main__":
    main()
