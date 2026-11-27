# math-findings

Short research notes on open problems in combinatorics, each with a full proof and a verifier that recomputes every computer-assisted claim.

| Note | New results | Paper | Check |
|---|---|---|---|
| [Integer domination roots](domination-roots/) | −10 is a domination root, so the integer roots in [−11, −1] are exactly −2, −4, −6, −8, −10; no graph on at most 12 vertices has root −4, so the least order is 13, 14 or 15 | [PDF](domination-roots/paper.pdf) | `python3 -B domination-roots/verify.py` |
| [Three-move subtraction games](subtraction-games/) | Manabe's Conjecture 52 for all b ≤ 2a, and for all b ≤ 250 with any c (computer-assisted) | [PDF](subtraction-games/paper.pdf) | `python3 -B subtraction-games/verify.py` |

## How these were produced

The results were found with AI systems (language-model reasoning and computer search) directed by the author. A result is included only after an independent re-implementation of every computation, an adversarial review of every proof by a different model, and a search for prior work. Each paper states which parts are computer-assisted and credits the earlier work it builds on.

Corrections are welcome: open an issue.

## Licence

Text: CC BY 4.0. Code and data: MIT. See `LICENSE`.
