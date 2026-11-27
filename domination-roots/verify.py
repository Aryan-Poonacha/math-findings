#!/usr/bin/env python3
"""Recompute every computer-assisted claim of paper.pdf from the definitions.

    python3 -B verify.py                      default mode, a few seconds
    python3 -B verify.py --full NAUTY_DIR     also rerun the exhaustive search
                                              through order 12 (about an hour
                                              on 15 threads)

NAUTY_DIR is a built nauty 2.8.9 source tree:
    curl -O https://pallini.di.uniroma1.it/nauty2_8_9.tar.gz   # sha256 c97ab42b...52c4ee
    tar xzf nauty2_8_9.tar.gz && cd nauty2_8_9 && ./configure && make
--full needs gcc. Each finished slice is kept in rerun/ with the sha256 of the
binary that produced it and its arguments, so an interrupted run resumes where
it stopped; a slice whose identity differs is recomputed, and every record,
fresh or kept, passes the same validation.

Standard library only. Input: data/graphs.json (the graphs) and
data/G10-edges.txt (the 373-vertex graph as a plain edge list). Nothing is read
as a result: every polynomial, value and count is recomputed here, and the
generated file tables.tex (included by paper.tex) must equal what this script
produces. Any mismatch raises an exception.

D(G,x) is the sum of x^|S| over all dominating sets S of G (S dominates G when
every vertex is in S or has a neighbour in S).
"""
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))

# OEIS A001349: number of connected graphs on n unlabelled vertices, n = 0..11.
CONNECTED = [1, 1, 1, 2, 6, 21, 112, 853, 11117, 261080, 11716571, 1006700565]


def check(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------- graphs

def graph(n, edges):
    """Open neighbourhoods as bit masks; rejects loops, repeated edges, bad labels."""
    adj = [0] * n
    for u, v in edges:
        check(0 <= u < n and 0 <= v < n and u != v, f'bad edge {u}-{v}')
        check(not adj[u] >> v & 1, f'repeated edge {u}-{v}')
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def edge_list(adj):
    return [(u, v) for u in range(len(adj)) for v in range(u + 1, len(adj)) if adj[u] >> v & 1]


def connected(adj):
    seen, stack = 1, [0]
    while stack:
        new = adj[stack.pop()] & ~seen
        seen |= new
        stack.extend(v for v in range(len(adj)) if new >> v & 1)
    return seen == (1 << len(adj)) - 1


def is_support(adj, v):
    """v is adjacent to a vertex of degree one."""
    return any(adj[v] >> u & 1 and adj[u].bit_count() == 1 for u in range(len(adj)))


def hub_graph(spec):
    """Vertex 0 is the hub h, vertices 1..leaves its pendant leaves, then the
    blocks in data order, each copy relabelled consecutively; h is joined to the
    attach vertices of every copy."""
    edges = [(0, i) for i in range(1, spec['leaves'] + 1)]
    n = spec['leaves'] + 1
    for block in spec['blocks']:
        for _ in range(block['copies']):
            edges += [(n + u, n + v) for u, v in block['edges']]
            edges += [(0, n + t) for t in block['attach']]
            n += block['order']
    return graph(n, edges)


def join_supports(parts):
    """Disjoint union of (adj, support) pairs plus a path through the supports."""
    adj, supports = [], []
    for part, s in parts:
        offset = len(adj)
        adj += [mask << offset for mask in part]
        supports.append(offset + s)
    for a, b in zip(supports, supports[1:]):
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    return adj


def star(t):
    """K_{1,t}, centre 0."""
    return graph(t + 1, [(0, i) for i in range(1, t + 1)])


# ------------------------------------------------------------ polynomials

def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    return [c + (b[i] if i < len(b) else 0) for i, c in enumerate(a)]


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, c in enumerate(a):
        if c:
            for j, d in enumerate(b):
                out[i + j] += c * d
    return out


def ppow(a, e):
    out = [1]
    for _ in range(e):
        out = pmul(out, a)
    return out


def value(p, x):
    out = 0
    for c in reversed(p):
        out = out * x + c
    return out


def derivative(p):
    return [i * c for i, c in enumerate(p)][1:]


def multiplicity(p, r):
    """Multiplicity of the root r of the nonzero integer polynomial p."""
    m = 0
    while value(p, r) == 0:
        quotient, carry = [], 0          # synthetic division by (x - r)
        for c in reversed(p[1:]):
            carry = carry * r + c
            quotient.append(carry)
        p = quotient[::-1]
        m += 1
    return m


# -------------------------------------------- three ways to compute D(G,x)

def brute_force(adj, external=0):
    """Coefficients of sum x^|S| over all S whose closed neighbourhood covers
    every vertex outside `external`. external=0 gives D(G,x); external=T gives
    Q_{G,T}(x). Direct enumeration of all 2^n subsets."""
    n = len(adj)
    full = (1 << n) - 1
    closed = [adj[v] | 1 << v for v in range(n)]
    cover = [0] * (1 << n)
    counts = [0] * (n + 1)
    counts[0] = int(external == full)
    for s in range(1, 1 << n):
        low = s & -s
        cover[s] = cover[s ^ low] | closed[low.bit_length() - 1]
        if cover[s] | external == full:
            counts[s.bit_count()] += 1
    return counts


IN, DOM, UND = 0, 1, 2   # in S / not in S but dominated / not yet dominated


def frontier_dp(adj):
    """D(G,x) from the adjacency alone, by a dynamic program over a vertex order.

    Vertices are processed one at a time. A processed vertex with an unprocessed
    neighbour is 'open' and carries one of the statuses IN, DOM, UND; when its
    last neighbour has been processed it closes, and states in which it is still
    UND are discarded. A polynomial sum c_k x^k (0 <= c_k < 2^n) is stored as
    the integer sum c_k 2^(B k), B = n + 1, so multiplying by x is a shift.
    Returns (coefficients, largest number of states)."""
    n = len(adj)
    nbrs = [[u for u in range(n) if adj[v] >> u & 1] for v in range(n)]
    B = n + 1
    remaining = [len(nb) for nb in nbrs]
    done = [False] * n
    is_open = [False] * n
    open_list, states, widest = [], {(): 1}, 1
    candidates = {max(range(n), key=lambda v: len(nbrs[v]))}
    for _ in range(n):
        if not candidates:
            candidates = {next(v for v in range(n) if not done[v])}

        def growth(v):   # change in the number of open vertices, then tie-breaks
            closes = sum(1 for u in nbrs[v] if is_open[u] and remaining[u] == 1)
            return (int(remaining[v] > 0) - closes, -sum(done[u] for u in nbrs[v]), v)
        v = min(candidates, key=growth)
        candidates.discard(v)
        done[v] = True
        position = [i for i, u in enumerate(open_list) if adj[v] >> u & 1]
        new = defaultdict(int)
        for st, w in states.items():
            new[st + ((DOM if any(st[i] == IN for i in position) else UND),)] += w
            chosen = list(st)
            for i in position:
                if chosen[i] == UND:
                    chosen[i] = DOM
            new[tuple(chosen) + (IN,)] += w << B
        open_list.append(v)
        for u in nbrs[v]:
            remaining[u] -= 1
            if not done[u]:
                candidates.add(u)
        keep = [i for i, u in enumerate(open_list) if remaining[u] > 0]
        shut = [i for i, u in enumerate(open_list) if remaining[u] == 0]
        states = defaultdict(int)
        for st, w in new.items():
            if all(st[i] != UND for i in shut):
                states[tuple(st[i] for i in keep)] += w
        open_list = [open_list[i] for i in keep]
        for u in range(n):
            is_open[u] = False
        for u in open_list:
            is_open[u] = True
        widest = max(widest, len(states))
    total = states[()]
    mask = (1 << B) - 1
    return [(total >> (B * k)) & mask for k in range(n + 1)], widest


def hub_formula(spec):
    """Lemma 1: D(G,x) = x^l prod D(F_i,x)^c_i + x (1+x)^l prod Q_i(x)^c_i,
    with D(F_i,x) and Q_i = Q_{F_i,T_i} counted by brute force on each block."""
    leaves = spec['leaves']
    absent = [0] * leaves + [1]
    present = pmul([0, 1], ppow([1, 1], leaves))
    for block in spec['blocks']:
        f = graph(block['order'], block['edges'])
        attach = sum(1 << t for t in block['attach'])
        absent = pmul(absent, ppow(brute_force(f), block['copies']))
        present = pmul(present, ppow(brute_force(f, attach), block['copies']))
    return padd(absent, present)


def extension_values(adj, x):
    """Lemma 3 at the point x for every T: the list whose entry T (T > 0) is
    D(H + v_T, x) = D(H,x) - U_{H,T}(x) + x Q_{H,T}(x), where v_T is a new
    vertex with neighbourhood T. Two subset-sum transforms, as in extproc.c."""
    n = len(adj)
    full = (1 << n) - 1
    closed = [adj[v] | 1 << v for v in range(n)]
    cover = [0] * (1 << n)
    q = [0] * (1 << n)           # indexed by the set of undominated vertices
    u = [0] * (1 << n)           # indexed by the dominating set itself
    q[full] = 1
    for s in range(1, 1 << n):
        low = s & -s
        cover[s] = cover[s ^ low] | closed[low.bit_length() - 1]
        w = x ** s.bit_count()
        q[full ^ cover[s]] += w
        if cover[s] == full:
            u[s] = w
    d = q[0]
    for i in range(n):
        bit = 1 << i
        for s in range(1 << n):
            if s & bit:
                q[s] += q[s ^ bit]
                u[s] += u[s ^ bit]
    return [None] + [d - u[full ^ t] + x * q[t] for t in range(1, 1 << n)]


# --------------------------------------------------------- small orders

def canonical(adj):
    """Canonical form: colour refinement, then the least adjacency code over all
    orderings that list the colour classes in increasing order."""
    n = len(adj)
    colours = [adj[v].bit_count() for v in range(n)]
    while True:
        signature = [(colours[v], tuple(sorted(colours[u] for u in range(n) if adj[v] >> u & 1)))
                     for v in range(n)]
        rank = {s: i for i, s in enumerate(sorted(set(signature)))}
        refined = [rank[s] for s in signature]
        if len(set(refined)) == len(set(colours)):
            break
        colours = refined
    cells = [[v for v in range(n) if refined[v] == c] for c in range(len(rank))]
    best = None
    for parts in itertools.product(*(itertools.permutations(c) for c in cells)):
        order = [v for part in parts for v in part]
        code = 0
        for i in range(n):
            row = adj[order[i]]
            for j in range(i + 1, n):
                code = 2 * code + (row >> order[j] & 1)
        if best is None or code < best:
            best = code
    return (n, tuple(sorted(refined)), best)


def connected_graphs(max_order):
    """All connected graphs of order 1..max_order up to isomorphism: every
    connected graph of order n+1 is H + v with H connected of order n (delete a
    leaf of a spanning tree)."""
    levels = {1: [[0]]}
    for n in range(1, max_order):
        found = {}
        for h in levels[n]:
            for t in range(1, 1 << n):
                g = [mask | (t >> i & 1) << n for i, mask in enumerate(h)] + [t]
                found.setdefault(canonical(g), g)
        levels[n + 1] = list(found.values())
    return levels


# ------------------------------------------------------------- the claims

def load():
    with open(os.path.join(HERE, 'data', 'graphs.json')) as f:
        data = json.load(f)
    graphs = {name: graph(data[name]['order'], data[name]['edges']) for name in ('G4', 'G6', 'G8')}
    centres = {name: data[name]['center'] for name in ('G4', 'G6', 'G8')}
    edges = []
    with open(os.path.join(HERE, 'data', 'G10-edges.txt')) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                u, v = map(int, line.split())
                edges.append((u, v))
    return data, graphs, centres, edges


def factorised(m):
    sign, m, parts, p = ('-' if m < 0 else ''), abs(m), [], 2
    while p * p <= m:
        e = 0
        while m % p == 0:
            m //= p
            e += 1
        if e:
            parts.append(f'{p}^{{{e}}}' if e > 1 else str(p))
        p += 1
    if m > 1:
        parts.append(str(m))
    return sign + r'\cdot '.join(parts)


def theorem1(data, edges, tex):
    spec = data['G10']
    adj = hub_graph(spec)
    check(edge_list(adj) == sorted(tuple(sorted(e)) for e in edges),
          'data/G10-edges.txt is not the graph described in data/graphs.json')
    file_adj = graph(373, edges)
    n, m = len(file_adj), len(edges)
    check((n, m) == (373, 673) and connected(file_adj), 'G10: order, size or connectivity')
    poly, widest = frontier_dp(file_adj)                  # method A: adjacency only
    check(poly == hub_formula(spec), 'G10: Lemma 1 disagrees with the frontier DP')  # method B
    check(poly[-1] == 1 and all(c >= 0 for c in poly), 'G10: coefficients')
    check(value(poly, -10) == 0, 'G10: D(G,-10) != 0')
    roots = {r: multiplicity(poly, r) for r in range(-11, 0)}
    check({r: e for r, e in roots.items() if e} == {-10: 1, -2: 16}, f'G10: integer roots {roots}')
    gamma = next(k for k, c in enumerate(poly) if c)
    # Table 1 and the identity 10^(l-1) prod P_i^c_i = 9^l prod Q_i^c_i (Lemma 1 at x=-10).
    k, leaves = 10, spec['leaves']
    lhs, rhs, rows, blocks = k ** (leaves - 1), (k - 1) ** leaves, [], 0
    for i, block in enumerate(spec['blocks'], 1):
        f = graph(block['order'], block['edges'])
        check(connected(f), f'block {i} is not connected')
        p = value(brute_force(f), -k)
        q = value(brute_force(f, sum(1 << t for t in block['attach'])), -k)
        lhs *= p ** block['copies']
        rhs *= q ** block['copies']
        blocks += block['copies']
        rows.append(f"{i} & {block['copies']} & {block['order']} & ${factorised(p)}$ & ${factorised(q)}$ \\\\")
    check(lhs == rhs, 'G10: product identity of Table 1 fails')
    for prime, name in ((2, 'TenTwo'), (3, 'TenThree'), (5, 'TenFive')):
        e = 0
        while lhs % prime ** (e + 1) == 0:
            e += 1
        tex[name] = e
    tex['TenOrder'], tex['TenEdges'], tex['TenLeaves'] = n, m, leaves
    tex['TenBlocks'], tex['TenTypes'], tex['TenGamma'] = blocks, len(spec['blocks']), gamma
    tex['BlockRows'] = '\n'.join(rows)
    tex['BlockEdgeRows'] = '\n'.join(
        f"{i} & $\\{{{','.join(map(str, b['attach']))}\\}}$ & {' '.join(f'{u}{v}' for u, v in b['edges'])} \\\\"
        for i, b in enumerate(spec['blocks'], 1))
    print(f'PASS Theorem 1: G10 connected, n=373, m=673; D(G10,x) by frontier DP (max {widest} states) '
          f'equals Lemma 1; -10 simple, -2 of multiplicity 16, no other root in [-11,-1]; Table 1 identity holds')
    return poly


def corollary2(seeds):
    k2 = brute_force(graph(2, [(0, 1)]))
    check(k2 == [0, 2, 1], 'K2')
    for adj, poly in [(graph(2, [(0, 1)]), k2)] + list(seeds.values()):
        check(connected(adj) and value(poly, 1) % 2 == 1, 'parity sanity check')
    for k, (adj, poly) in seeds.items():
        check(value(poly, -k) == 0, f'no root -{k}')
    print('PASS Corollary 2: connected witnesses for -2 (K2), -4, -6, -8 (MathArena), -10 (G10); '
          'D(G,1) odd for each, as Brouwer-Csorba-Schrijver predict')


def theorem3_small():
    """Orders <= 8 in pure Python (the paper's order-12 search, scaled down)."""
    levels = connected_graphs(7)
    check([len(levels[n]) for n in range(1, 8)] == CONNECTED[1:8], 'graph counts differ from A001349')
    zeros_at_minus2, transform_zeros = [], 0
    for n in range(1, 8):
        z = 0
        for g in levels[n]:
            poly = brute_force(g)
            check(all(value(poly, r) != 0 for r in (-4, -6, -8, -10)), f'order {n} graph with a root')
            z += value(poly, -2) == 0
            if n <= 5:                                   # Lemma 3 and its subset transform against brute force
                full = (1 << n) - 1
                fast = {x: extension_values(g, x) for x in (-2, -4)}
                for t in range(1, 1 << n):
                    avoiding = [0] * (n + 1)             # U_{H,T}: dominating sets of H disjoint from T
                    for s in range(1 << n):
                        covered = 0
                        for v in range(n):
                            if s >> v & 1:
                                covered |= g[v] | 1 << v
                        if not s & t and covered == full:
                            avoiding[s.bit_count()] += 1
                    formula = padd(padd(poly, [-c for c in avoiding]), [0] + brute_force(g, t))
                    h = [mask | (t >> i & 1) << n for i, mask in enumerate(g)] + [t]
                    direct = brute_force(h)
                    check(formula == direct, 'Lemma 3 fails')
                    for x, values in fast.items():
                        check(values[t] == value(direct, x), 'subset transform disagrees with brute force')
                        transform_zeros += values[t] == 0
        zeros_at_minus2.append(z)
    check(transform_zeros > 0, 'the subset transform found no zero')
    for g in levels[7]:
        check(all(v != 0 for v in extension_values(g, -4)[1:]), 'order-8 graph with root -4')
    print('PASS Theorem 3 (scaled down): connected graphs of orders 1-7 = A001349; none has a root at '
          '-4, -6, -8, -10; none of order 8 has root -4 (Lemma 3 by subset transform; Lemma 3 and the '
          f'transform at -2 and -4 equal brute force on every extension of order <= 6, {transform_zeros} of '
          'them zeros)')
    return zeros_at_minus2


def theorem4(data, g4, tex):
    spec = data['G16']
    adj = hub_graph(spec)
    poly = brute_force(adj)
    check(len(adj) == 16 and len(edge_list(adj)) == 25 and connected(adj), 'G16 structure')
    check(poly == hub_formula(spec), 'G16: Lemma 1')
    check(value(poly, -4) == 0 and value(derivative(poly), -4) != 0, 'G16: root -4 not simple')
    block = spec['blocks'][0]
    f = graph(block['order'], block['edges'])
    check(connected(f) and block['attach'], 'G16: block not connected or attachment set empty')
    p = value(brute_force(f), -4)
    q = value(brute_force(f, sum(1 << t for t in block['attach'])), -4)
    check(4 ** (spec['leaves'] - 1) * p == 3 ** spec['leaves'] * q, 'G16: Lemma 1 condition')
    tex['SixteenProd'] = 3 ** spec['leaves'] * q
    g4poly = brute_force(g4)
    check(is_support(g4, data['G4']['center']) and is_support(star(1), 0), 'supports')
    check(len(g4) == 15 and connected(g4) and value(g4poly, -4) == 0
          and value(derivative(g4poly), -4) != 0, 'G4')
    for s in range(2, 6):                       # orders 17..20: G4 joined to a star on s vertices
        joined = join_supports([(g4, data['G4']['center']), (star(s - 1), 0)])
        jp, _ = frontier_dp(joined)
        check(len(joined) == 15 + s and connected(joined), 'star join structure')
        check(jp == pmul(g4poly, brute_force(star(s - 1))), 'Lemma 2 fails on a star join')
        check(multiplicity(jp, -4) == 1, 'star join: root -4 not simple')
    tex['SixteenP'], tex['SixteenQ'] = p, q
    tex['SixteenEdges'] = r',\ '.join(f'{u}{v}' for u, v in block['edges'])
    tex['SixteenT'] = r'\{' + ','.join(map(str, block['attach'])) + r'\}'
    tex['FourEdges'] = ', '.join(f'$\\{{{u},{v}\\}}$' for u, v in edge_list(g4))
    terms = [(f'{c}' if c != 1 else '') + (f'x^{{{i}}}' if i > 1 else 'x' if i == 1 else '')
             for i, c in reversed(list(enumerate(poly))) if c]
    tex['SixteenPoly'] = ' + '.join(terms[:7]) + r' \\ + ' + ' + '.join(terms[7:])
    tex['SixteenDeriv'] = value(derivative(poly), -4)
    tex['FourDeriv'] = value(derivative(g4poly), -4)
    print("PASS Theorem 4: G16 connected, n=16, m=25, D(G16,-4)=0 with D'(G16,-4)="
          f"{tex['SixteenDeriv']} (brute force over 2^16 sets, equals Lemma 1); G4 joined to stars "
          'gives orders 17-20 with a simple root -4 and D multiplying (Lemma 2)')


def theorem5(data, seeds, centres, tex):
    p = 1009
    rows = []
    leaves = {k: sum(adj[centres[k]] >> u & 1 and adj[u].bit_count() == 1 for u in range(len(adj)))
              for k, (adj, _) in seeds.items()}
    check(leaves == {4: 5, 6: 10, 8: 16, 10: 22}, f'leaves at the designated supports: {leaves}')
    for k, (adj, poly) in seeds.items():
        check(is_support(adj, centres[k]), f'designated vertex of G{k} is not a support')
        check(multiplicity(poly, -k) == 1, f'G{k}: root -{k} is not simple')
        residues = []
        for j in seeds:
            v = value(poly, -j)
            check((v == 0) == (j == k), f'D_{k}(-{j})')
            check(j == k or v % p != 0, 'choose another modulus')
            residues.append(v % p)
        d = value(derivative(poly), -k) % p
        check(d != 0, 'choose another modulus')
        rows.append(f"$G_{{{k}}}$ & {len(adj)} & " + ' & '.join(map(str, residues)) + f' & {d} \\\\')
    tex['ResidueRows'] = '\n'.join(rows)
    # the construction on a small instance: (m4, m6, m8, m10) = (2, 1, 1, 1)
    target = {4: 2, 6: 1, 8: 1, 10: 1}
    parts = [(seeds[k][0], centres[k]) for k in seeds for _ in range(target[k])]
    joined = join_supports(parts)
    poly, _ = frontier_dp(joined)
    product = [1]
    for k in seeds:
        product = pmul(product, ppow(seeds[k][1], target[k]))
    check(connected(joined) and poly == product, 'Lemma 2 fails on the multiplicity instance')
    check({k: multiplicity(poly, -k) for k in seeds} == target, 'multiplicities')
    print(f'PASS Theorem 5: G4, G6, G8, G10 have a support, a simple root at their own point and '
          f'none at the other three (residues mod {p}); the {len(joined)}-vertex join for (2,1,1,1) '
          'has D = product and multiplicities exactly (2,1,1,1)')


def self_test():
    """The frontier DP against brute force on random graphs."""
    rng = random.Random(1)
    for _ in range(300):
        n = rng.randint(1, 11)
        density = rng.random()
        adj = graph(n, [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < density])
        check(frontier_dp(adj)[0] == brute_force(adj), 'frontier DP disagrees with brute force')
    print('PASS self-test: frontier DP equals brute force on 300 random graphs of order <= 11')


def tables(tex):
    out = ['% Generated by verify.py --write-tables from data/; do not edit.']
    for key, val in tex.items():
        body = f'%\n{val}' if key.endswith('Rows') else val
        out.append(f'\\newcommand{{\\{key}}}{{{body}}}')
    return '\n'.join(out) + '\n'


# ------------------------------------------------------------ --full mode

# The fields, in order, of the RESULT line each search program prints at exit.
RESULT_FIELDS = {'geng_dom': ['n', 'graphs'] + [f'z{r}' for r in range(2, 11)],
                 'geng_ext4': ['base_n', 'bases', 'tests', 'zeros', 'tpl', 'root']}


def parse_result(lines, program, n, res, mod):
    """Counts of one search slice from its record lines, fresh or checkpointed:
    exactly one line, a RESULT line with the expected fields, order and root,
    nonnegative counts, and tests = bases (2^n - 1) for the extensions."""
    where = f'{program} n={n} slice {res}/{mod}'
    check(len(lines) == 1 and lines[0].startswith('RESULT '), f'{where}: expected one RESULT line, got {lines}')
    pairs = [field.split('=') for field in lines[0].split()[1:]]
    check([p[0] for p in pairs] == RESULT_FIELDS[program]
          and all(len(p) == 2 and re.fullmatch(r'-?[0-9]+', p[1]) for p in pairs),
          f'{where}: malformed record {lines[0]!r}')
    counts = {key: int(val) for key, val in pairs}
    check(counts.pop('n' if program == 'geng_dom' else 'base_n') == n, f'{where}: wrong order')
    if program == 'geng_ext4':
        check(counts.pop('root') == -4, f'{where}: wrong root')
        check(counts['tests'] == counts['bases'] * (2 ** n - 1), f'{where}: tests != bases (2^n - 1)')
    check(all(c >= 0 for c in counts.values()), f'{where}: negative count')
    return counts


def full(nauty, small_zeros, g4):
    """Rebuild domproc/extproc against nauty 2.8.9 and rerun the exhaustive
    search: brute force for orders 1..11, Lemma 3 extensions for order 12.
    Each slice's record is checkpointed in rerun/ under the sha256 of the
    binary and the slice arguments; rerun the command to resume. A checkpoint
    with another identity is regenerated; one with this identity whose record
    is malformed stops the run."""
    build, out = os.path.join(HERE, 'build'), os.path.join(HERE, 'rerun')
    os.makedirs(build, exist_ok=True)
    flags = ['gcc', '-O3', '-mpopcnt', '-march=native']
    # Positive control for the C extension check: G4 minus its leaf 14 gains a root -4 by
    # restoring that leaf at the support 2, and every zero agrees with Lemma 3 in Python.
    subprocess.run(flags + ['-o', os.path.join(build, 'ext_standalone'), os.path.join(HERE, 'search', 'extproc.c')],
                   check=True)
    base = [mask & ~(1 << 14) for mask in g4[:14]]
    edges = edge_list(base)
    stdin = ' '.join(map(str, [14, len(edges)] + [v for e in edges for v in e]))
    stdout = subprocess.run([os.path.join(build, 'ext_standalone')], input=stdin, capture_output=True,
                            text=True, check=True).stdout
    c_zeros = {int(line.split('=')[1]) for line in stdout.splitlines() if line.startswith('zero at T=')}
    python_zeros = {t for t, v in enumerate(extension_values(base, -4)) if v == 0}
    check(g4[14] == 1 << 2 and 1 << 2 in c_zeros and c_zeros == python_zeros
          and stdout.splitlines()[-1] == f'zeros={len(c_zeros)}',
          'extproc.c misses the G4 extension or disagrees with Lemma 3')
    print(f'  extproc.c positive control: {len(c_zeros)} zero(s) on G4 minus a leaf, including G4 itself; '
          'equal to Lemma 3 in Python', flush=True)
    objects = [os.path.join(nauty, o) for o in
               ('gtoolsW.o', 'nautyW1.o', 'nautilW1.o', 'naugraphW1.o', 'schreier.o', 'naurng.o')]
    fingerprint = {}
    for name, proc in (('geng_dom', 'domproc'), ('geng_ext4', 'extproc')):
        subprocess.run(flags + ['-DMAXN=WORDSIZE', '-DWORDSIZE=32', f'-DOUTPROC={proc}', '-I', nauty,
                                '-o', os.path.join(build, name), os.path.join(nauty, 'geng.c'),
                                os.path.join(HERE, 'search', proc + '.c'), *objects],
                       check=True)
        with open(os.path.join(build, name), 'rb') as f:
            fingerprint[name] = hashlib.sha256(f.read()).hexdigest()
    jobs = max(1, (os.cpu_count() or 2) - 1)

    def run(program, n, res, mod):
        identity = f'# {program} sha256={fingerprint[program]} n={n} res={res} mod={mod}'
        log = os.path.join(out, f'{program}-n{n}', f'res-{res}.log')
        if os.path.exists(log):
            with open(log) as f:
                lines = f.read().splitlines()
            if lines[:1] == [identity]:
                return parse_result(lines[1:], program, n, res, mod)
        os.makedirs(os.path.dirname(log), exist_ok=True)
        args = [os.path.join(build, program), '-cq', str(n)] + ([f'{res}/{mod}'] if mod > 1 else [])
        result = subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, check=True)
        lines = [l for l in result.stderr.splitlines() if l.startswith(('RESULT', 'ZERO'))]
        counts = parse_result(lines, program, n, res, mod)
        with open(log + '.tmp', 'w') as f:
            f.write('\n'.join([identity] + lines) + '\n')
        os.replace(log + '.tmp', log)
        return counts

    def totals(program, n, mod, label):
        start = time.time()
        with ThreadPoolExecutor(jobs) as pool:
            slices = list(pool.map(lambda r: run(program, n, r, mod), range(mod)))
        tot = {key: sum(s[key] for s in slices) for key in slices[0]}
        print(f'  {label}: {tot}  ({time.time() - start:.0f} s)', flush=True)
        return tot

    for n in range(1, 12):
        tot = totals('geng_dom', n, 1 if n <= 9 else 64, f'order {n}')
        check(tot['graphs'] == CONNECTED[n], f'order {n}: geng count differs from A001349')
        check(all(tot[f'z{r}'] == 0 for r in (3, 4, 5, 6, 7, 8, 9, 10)), f'order {n}: a root found')
        if n <= 7:
            check(tot['z2'] == small_zeros[n - 1], f'order {n}: C and Python disagree at -2')
    tot = totals('geng_ext4', 11, 512, 'order 12 (one-vertex extensions of order 11)')
    check(tot['bases'] == CONNECTED[11] and tot['tests'] == CONNECTED[11] * 2047, 'order 12 coverage')
    check(tot['zeros'] == 0, 'order 12: a graph with root -4')
    print('PASS Theorem 3 (full): no connected graph of order <= 11 has a root at -4, -6, -8, -10; '
          'no connected graph of order 12 has root -4')


def main():
    if not __debug__:
        raise SystemExit('run without -O')
    if not (len(sys.argv) == 1 or sys.argv[1:] == ['--write-tables']
            or (len(sys.argv) == 3 and sys.argv[1] == '--full')):
        raise SystemExit(__doc__)
    start = time.time()
    data, graphs, centres, edges = load()
    tex = {}
    self_test()
    g10 = theorem1(data, edges, tex)
    seeds = {4: graphs['G4'], 6: graphs['G6'], 8: graphs['G8']}
    seeds = {k: (adj, frontier_dp(adj)[0]) for k, adj in seeds.items()}
    seeds[10] = (graph(373, edges), g10)
    corollary2(seeds)
    small_zeros = theorem3_small()
    theorem4(data, graphs['G4'], tex)
    theorem5(data, seeds, {4: centres['G4'], 6: centres['G6'], 8: centres['G8'], 10: 0}, tex)
    generated = tables(tex)
    path = os.path.join(HERE, 'tables.tex')
    if sys.argv[1:] == ['--write-tables']:
        with open(path, 'w') as f:
            f.write(generated)
    with open(path) as f:
        check(f.read() == generated, 'tables.tex is stale: rerun with --write-tables')
    print('PASS tables.tex equals the values recomputed here')
    print(f'default mode done in {time.time() - start:.1f} s')
    if sys.argv[1:2] == ['--full']:
        full(os.path.abspath(sys.argv[2]), small_zeros, graphs['G4'])


if __name__ == '__main__':
    main()
