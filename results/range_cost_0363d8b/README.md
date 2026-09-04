# Range-cost re-measurement at the release build (C-1)

Three paired runs of the human primary FASTA at commit `0363d8b`, `--max-period` 100
against 2000, four threads, the manuscript's F configuration. Both arms of a pair run
sequentially inside one SLURM job, so a pair always shares a node.

## Why this was re-run

Every published range-cost pair — the original and all three replicates — ran at commit
`07ad6fa`. That is the only commit where Tier 3's search window was fixed at 100 bp-100 kb
regardless of `--max-period`, so the narrow arm performed the full long-period search and
discarded the result. At `0363d8b` the ceiling bounds Tier 3 as well
(`bwtandem/finder.py:116`, whose own comment reads "This is NOT output-neutral").

## Results

| rep | job | node | p100 | p2000 | ratio | p100 MaxRSS | p2000 MaxRSS |
|---|---|---|--:|--:|--:|--:|--:|
| 1 | 6147698 | cpu-15 | 4:01:05 | 7:18:21 | 1.818 | 17.21 GiB | 17.37 GiB |
| 2 | 6147699 | cpu-28 | 4:06:38 | 7:16:41 | 1.771 | 17.22 GiB | 17.37 GiB |
| 3 | 6147700 | cpu-28 | 4:01:22 | 7:07:45 | 1.772 | 17.22 GiB | 17.38 GiB |

Mean ratio **1.79**, range 1.77-1.82.

Call counts are identical across all three replicates: **3,993,151** at p100 and
**4,014,108** at p2000. The p100 count matches the native `--max-period 100` run behind
Table 1d's F row, which is an independent check that the same configuration produced the
same output. All six arms exited 0.

## What changed against the published pairs

| | published (`07ad6fa`) | this run (`0363d8b`) |
|---|--:|--:|
| p100 | 6.61 / 6.58 / 5.50 h | 4.02 / 4.11 / 4.02 h |
| p2000 | 8.61 / 8.58 / 7.75 h | 7.31 / 7.28 / 7.13 h |
| ratio | 1.30 / 1.30 / 1.41 (mean 1.34) | 1.82 / 1.77 / 1.77 (mean **1.79**) |

The narrow arm is what moved: 6.6 h to 4.0 h. At `07ad6fa` it was paying for a
long-period search whose output it then threw away; at `0363d8b` it does not run that
search. The wide arm fell less, 8.6 h to 7.3 h. So the published 1.3-1.4x compared a
wide run against a narrow run that was not actually narrow.

## Caveats

- Replicates 2 and 3 shared cpu-28; replicate 1 had cpu-15 to itself, and it is the
  highest of the three at 1.82. The published replicates were co-scheduled more heavily
  still (reps 1 and 2 ran all four arms simultaneously on cpu-51), so this design is
  cleaner than the one it replaces but is not fully isolated.
- Both arms of a pair run sequentially here; in the published pairs they ran
  concurrently. The spread is correspondingly narrower, 1.77-1.82 against 1.30-1.41.
- Memory is effectively unchanged by the ceiling: 17.21 GiB at p100, 17.37 GiB at p2000.
