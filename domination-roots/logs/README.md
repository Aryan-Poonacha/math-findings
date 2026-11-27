Logs of the exhaustive search (Remark 2 of paper.pdf). One C implementation (`search/domproc.c`, `search/extproc.c`), run twice.

- `audit-orders-1-10.log`, `audit-order-11.log`, `audit-order-12.log`: the first run
  (2026-11-27, nauty 2.8.9, AMD Ryzen 9 4900H, 15 parallel jobs;
  order 11 took 535 s, order 12 took 3518 s).
- `full-run.log`: the second run, `python3 -B verify.py --full NAUTY_DIR` on the same machine,
  made by an earlier revision of `verify.py` (its PASS labels use the old theorem numbers).
- `full-run-slices.txt`: the per-slice RESULT lines of that run (no ZERO lines occurred).
  These runs predate the checkpoint identity line that `verify.py` now writes. Every
  line passes the current `parse_result` for its program and slice, and with these
  records as checkpoints under the identity of the binaries below, `verify.py --full`
  reproduces the totals of `full-run.log` without rerunning.

Both runs used the C sources with the hashes below. The current `search/extproc.c` differs
from them only in the comment on line 7 (sha256 2219af66…a854e); compiled with the
commands below it gives the same `geng_ext4` binary, byte for byte (checked 2026-11-27).

Build (gcc 13.3.0, Ubuntu 24.04, x86-64; N = nauty 2.8.9 source tree after
`./configure && make`), the commands `verify.py --full` runs:

    gcc -O3 -mpopcnt -march=native -DMAXN=WORDSIZE -DWORDSIZE=32 -DOUTPROC=domproc -I N \
        -o geng_dom N/geng.c search/domproc.c \
        N/gtoolsW.o N/nautyW1.o N/nautilW1.o N/naugraphW1.o N/schreier.o N/naurng.o
    gcc -O3 -mpopcnt -march=native -DMAXN=WORDSIZE -DWORDSIZE=32 -DOUTPROC=extproc -I N \
        -o geng_ext4 N/geng.c search/extproc.c \
        N/gtoolsW.o N/nautyW1.o N/nautilW1.o N/naugraphW1.o N/schreier.o N/naurng.o

SHA-256:

    c97ab42bf48796a86a598bce3e9269047ca2b32c14fc23e07208a244fe52c4ee  nauty2_8_9.tar.gz
    d1fcaa297797e5409b4b9060b376d686d6cbf97cbefd5c398d8b999764a2f152  search/domproc.c
    20f1814ec0ebf8c3643abdedbf6d767fb5daf1f9e42ac4bd32786b2df7b5ebac  search/extproc.c (as run)
    93b77215979bb9e42e2641d3e302b87126acc2326766cd3c5fd228862041fd69  geng_dom
    7d1f7fcfb9ed266d433eedfe57b05e48788fa58d4825fd4b3c5d78388a0e43d1  geng_ext4

The binary hashes depend on the compiler and on `-march=native` (here AMD Ryzen 9 4900H).
On that machine the commands reproduce them byte for byte, equal to the binaries of the first run.
A fresh, independent rerun must start with no `rerun/` directory.
