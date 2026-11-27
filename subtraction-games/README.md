# Manabe's pure-periodicity criterion for three-move subtraction games

Paper: [paper.pdf](paper.pdf) (source [paper.tex](paper.tex)). Not peer reviewed; AI use is disclosed in Section 8.

Setting: the subtraction game S = {a, b, c} with 2 ≤ a < b < c, gcd(a, b, c) = 1, c ≠ a + b, and Manabe's angle tests on ρ = c mod (a + b) ([arXiv:2609.05358v2](https://arxiv.org/abs/2609.05358v2), Conjecture 52).

New here:
- **Theorem 3.2.** Conjecture 52 holds for every b ≤ 2a and every c. If ρ fails the tests, three explicit positions rule out every pure period.
- **Lemma 4.1 (transfer) and Theorem 4.2.** For fixed a, b and c ≥ 2(a + b), the P-positions at c + j and 2c + j, −(a + b) ≤ j < 2(a + b), depend only on c mod (a + b); with a finite witness search, Conjecture 52 holds for every b ≤ 250 and every c (computer-assisted).
- **Proposition 5.1.** For any 1 ≤ a < b < c with c ≥ 2b (no gcd or additivity condition), if ρ fails the tests, the P-positions and the nim-sequence have no period p ≤ c + b from 0.

Known, cited not claimed: Conjecture 57 (Levinstein, Zenodo, 2026-11-25; H. Zhang, GitHub, 2026-11-26); move deletion (Miklós–Post Prop. 2.9, Manabe Prop. 25); the line b = 2a (S. Zhang Prop. 3.1); Moriwaki's exclusion of the three candidate periods for b ≤ 2a.

Verify (Python 3, standard library only; prints one PASS line per claim, fails on any mismatch):

    python3 -B verify.py          # about 10 s on 16 cores: all proofs on finite ranges, certificates for b ≤ 100
    python3 -B verify.py --full   # about 70 s on 16 cores: the complete certificate for b ≤ 250 (log in logs/)
