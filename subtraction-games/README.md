# Manabe's pure-periodicity criterion for three-move subtraction games

Paper: [paper.pdf](paper.pdf) (source [paper.tex](paper.tex)). Archived at [doi:10.5281/zenodo.24606431](https://doi.org/10.5281/zenodo.24606431). Not peer reviewed; AI use is acknowledged in the paper.

Setting: the subtraction game S = {a, b, c} with 2 ≤ a < b < c, gcd(a, b, c) = 1, c ≠ a + b, and Manabe's angle tests on ρ = c mod (a + b) ([arXiv:2609.05358v2](https://arxiv.org/abs/2609.05358v2), Conjecture 52). Sufficiency of the tests is Manabe's Theorem 34; the new part is necessity.

New here:
- **Theorem 3.2.** Conjecture 52 holds for every b ≤ 2a and every c (necessity is new; with Manabe's Theorem 34). If ρ fails the tests, three explicit positions rule out every pure period.
- **Lemma 4.1 (transfer) and Theorem 4.2.** For fixed a, b and c ≥ 2(a + b), the P-positions at c + j and 2c + j, −(a + b) ≤ j < 2(a + b), depend only on c mod (a + b); with a finite witness search, Conjecture 52 holds for every b ≤ 250 and every c (computer-assisted).
- **Proposition 5.1.** For 2 ≤ a < b < c with c ≥ 2b (no gcd or additivity condition), if ρ fails the tests, the P-positions and the nim-sequence have no period p ≤ c + b from 0.

Known, cited not claimed: Conjecture 57 (Levinstein, Zenodo, 2026-11-25; H. Zhang, GitHub, 2026-11-26) and Levinstein's partial period exclusions; move deletion (Miklós–Post Prop. 2.9, Manabe Prop. 25); the line b = 2a (S. Zhang Prop. 3.1); Moriwaki's exclusion of the three candidate periods for b ≤ 2a; exact words for families inside Theorems 3.2 and 4.2, such as {4,6,10k+1} (X. Zeng), {6,17,c} (Pan) and {3,8,11k+1} (C. Lyons).

Verify (Python 3, standard library only; prints one PASS line per check, fails on any mismatch):

    python3 -B verify.py          # default (fast), about 10 s on 16 cores
    python3 -B verify.py --full   # full, about 1 min on 16 cores (log in logs/)

Default mode: the generator against the mex rule; exact nim cycles against the criterion on Manabe's box and on 2,000 random admissible games with 40 < b ≤ 250; finite-range tests of Lemma 3.1, Theorem 3.2's witness, Lemma 4.1 and the angle facts of Proposition 5.1; the Theorem 4.2 certificate for b ≤ 100. Full mode: the same with larger ranges and the complete Theorem 4.2 certificate for b ≤ 250, asserting its four counts. Lemma 2.4, Theorem 3.2, Lemma 4.1 and Proposition 5.1 are proved by hand; the checks on them are finite-range tests only.
