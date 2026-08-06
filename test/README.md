# Test data for insize

A small hand-written alignment set (22 records, 2 contigs) that covers every
branch of the filtering logic in `src/insize.cpp`.

Run it with:

    cd src/
    make test

or directly:

    ./test/run_tests.sh

`run_tests.sh` converts `test.sam` to BAM with samtools (it prefers
`lib/samtools/samtools`, which the build already produces, and falls back to
whatever is on `$PATH`), runs `insize` over a matrix of options and diffs the
output against `expected/`. If no samtools is found it runs on the SAM directly
— htslib autodetects the format — and skips the BAM round-trip.

Override the binary or the converter if needed:

    INSIZE=/path/to/insize SAMTOOLS=/path/to/samtools ./test/run_tests.sh

## What is in test.sam

Paired reads report `abs(TLEN)` from **read1 only**; single-end reads report the
length of `SEQ`.

| read | what it is | reported |
|---|---|---|
| `pair_a` | proper pair, read1 forward, TLEN +100 | 100 |
| `pair_b` | proper pair, read1 reverse, TLEN -150 | 150 (sign is dropped) |
| `pair_c` | mapped pair *not* flagged proper, TLEN 200 | 200, dropped by `-p` |
| `pair_d` | mates on chr1/chr2, TLEN 0 | nothing (TLEN 0 is always skipped) |
| `pair_e` | proper pair, TLEN 30 | 30, below `-l 50` |
| `pair_f` | proper pair, TLEN 500 | 500, above `-L 300` |
| `pair_g` | read1 mapped, mate unmapped, TLEN 0 | nothing |
| `single_a` | mapped single-end, 75 bp | 75 |
| `single_b` | unmapped single-end, 60 bp | 60, dropped by `-m` |
| `single_c` | mapped single-end, reverse, 45 bp | 45 |
| `dup_se` | PCR duplicate (0x400) | nothing |
| `qcfail_se` | QC fail (0x200) | nothing |
| `secondary_se` | secondary (0x100) | nothing |
| `supp_se` | supplementary (0x800) | nothing |
| `pair_h` | proper pair also flagged duplicate | nothing |

So the default run yields `100 150 200 30 500 75 60 45`, in file order.

`test.sam` is committed as text on purpose: it is readable in a diff and easy to
extend. To add a case, append a record and regenerate the affected expected
file, e.g.

    ./src/insize -m test/test.sam > test/expected/mapped.txt

## Note on the `-v` output

`expected/verbose_90_110.txt` has 22 fields for a 21-value window
(`-l 90 -L 110`). That is what the code currently does: it prints `count[0]`
once, then loops from index 0 again, so the first bin is emitted twice. The
expected file records today's behaviour rather than endorsing it — if that
first column is ever removed, update this file along with the fix.

Also note `-v` allocates one counter per length in `[-l, -L]`, so it should
always be given an explicit `-L`; with the default maximum it would try to
allocate an array of `INT_MAX` entries.
