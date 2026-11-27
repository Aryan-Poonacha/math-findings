# math-findings

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24606431.svg)](https://doi.org/10.5281/zenodo.24606431)

Short research notes on open problems in combinatorics, each with a full proof and a verifier that recomputes every computer-assisted claim. Archived at [doi:10.5281/zenodo.24606431](https://doi.org/10.5281/zenodo.24606431) (concept DOI: [10.5281/zenodo.24606430](https://doi.org/10.5281/zenodo.24606430)).

| Note | New results | Paper | Default check (fast) | Full check |
|---|---|---|---|---|
| [Integer domination roots](domination-roots/) | −10 is a domination root, so the integer roots in [−11, −1] are exactly −2, −4, −6, −8, −10; roots −4, −6, −8, −10 together with any multiplicities | [PDF](domination-roots/paper.pdf) | `python3 -B domination-roots/verify.py` (about 1 s) | `python3 -B domination-roots/verify.py --full NAUTY_DIR` (about 1 h; reruns the order ≤ 12 search) |
| [Three-move subtraction games](subtraction-games/) | Manabe's Conjecture 52 for all b ≤ 2a, and for all b ≤ 250 with any c (computer-assisted) | [PDF](subtraction-games/paper.pdf) | `python3 -B subtraction-games/verify.py` (about 10 s; certificate for b ≤ 100) | `python3 -B subtraction-games/verify.py --full` (about 1 min; certificate for b ≤ 250) |

## Acknowledgements

AI systems were used to help find the constructions and proofs and to run the computations. The author verified all arguments and computations and takes full responsibility for all content.

## Checks

Each `verify.py` recomputes its paper's computer-assisted claims from the definitions with the Python standard library and stops on any mismatch. In the domination note, D(G10, x) is computed in two independent ways; its search through order 12 is one C implementation, run twice. Each paper states which parts are computer-assisted and credits the earlier work it builds on.

Corrections are welcome: open an issue.

## Licence

Text: CC BY 4.0. Code and data: MIT. See `LICENSE`.
