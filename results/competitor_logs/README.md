# Competitor GNU-time logs

Deposited 2026-09-03 for issue #30 and `quarantine.md` §6.20/§6.13.

Until now every competitor cost cell in the paper traced to these files, but the
files themselves lived only on the benchmarking cluster at
`/data/gpfs/assoc/pgl/filip/bwtandem_results/benchmarking_results/<tool>/logs/`.
A reader outside this filesystem had to take the transcription on trust. They are
small, so there was no reason for that.

Each file is `/usr/bin/time -v` output, copied verbatim, and records the command
line, wall clock, maximum resident set size and exit status. Names are
`<tool>__<original filename>`; the original names carry the assembly.

## What they establish

The human GCA_000001405.15 row of `comparator_baselines.md` recomputes from these
exactly: ULTRA 29:46:49 / 1,758,540 KiB, TRF 33:43:46 / 1,518,668 KiB, tantan
51:46.94 / 281,316 KiB, TRASH 107:37:36 / 15,298,244 KiB, mreps 54:41.08 /
6,691,220 KiB. Divide the kibibyte counts by 1024² for the GiB figures printed in
the tables.

## Two things they also show

The three Col-CEN TRASH runs are separate executions with separate costs —
CEN159 6:18:29 / 1,389,316 KiB, CEN178 25:22:40 / 2,761,680 KiB, de novo
5:47:30 / 1,357,252 KiB. Table 2's two TRASH rows currently carry the CEN159
figure for both, while the template row is scored from the union of all three
(`quarantine.md` §6.1).

The human TRASH log ends with exit status 1 after a Circos plotting error. The
failure appears to fall after data export, but the repository presented the run
as an ordinary completion (`quarantine.md`, Codex round-1 finding 12).
