# Quarantine — runtime cells only (2026-08-05, revised twice)

**Scope has been narrowed twice.** The first version quarantined every published
BWTandem figure. The second narrowed it to the runtime *and* memory cells. This
version narrows it again, to the **runtime cells alone** — the provenance of every
published cost number has now been traced to a specific job, and the memory cells
turned out not to be a provenance failure at all.

Read this before reusing any cost figure.

## The evidence that settled it

`filip/bwt/run.sbatch` echoes `HOST $(hostname)` and `date`, and the per-array SLURM
logs survive in `filip/bwt/logs/`. The row counts in those logs act as fingerprints,
so each published number can be attributed to an execution:

| array | date | nodes | hg38 rows | colcen rows | maize rows |
|---|---|---|--:|--:|--:|
| 5912536 | Jul 16 | cpu-29, cpu-15 | 3,364,455 | 88,034 | 1,640,603 |
| 5924706 | Jul 19 | cpu-37, cpu-7 | 4,441,404 | 170,949 | 2,840,913 |
| **5935102** | **Jul 21** | **cpu-64, cpu-78** | **3,994,477** | **161,187** | **2,396,215** |

The published BEDs are the **Jul 21** array. The published runtimes are the **Jul 16**
array's wall-clock: 7:21:43, 18:27 and 6:36:18 — i.e. the 7.4 h / 0.31 h / 6.6 h of
the tables, to the digit.

## 1. Runtime — a real defect, and the correct value was in the same directory

The runtime cells were paired with the wrong execution. The Jul 21 run recorded its
own cost via `/usr/bin/time -v`, in the same log as the BED it produced:

| genome | log | `Total repeats found` | BED lines | wall clock |
|---|---|--:|--:|--:|
| hg38 | `bwt_hg38.log` | 3,994,477 | 3,994,477 | 7:49:08 |
| colcen | `bwt_colcen.log` | 161,187 | 161,187 | 16:44.62 |
| maize | `bwt_maize.log` | 2,396,215 | 2,396,215 | 9:26:34 |

> ⚠️ An earlier version of this file claimed "each published run log records a call
> count that does not match the BED". **That was wrong** — it read the Jul 16 array
> log (`slurm_5912536_*.log`, `rows=3364455`) and the Jul 21 pipeline log
> (`bwt_hg38.log`, `Total repeats found: 3994477`) as if they were one file. The
> pipeline logs are internally consistent with their BEDs on all three genomes.

## 2. Memory — not a defect. Two metrics, one of which undercounts.

`src/main.py:162` uses `ProcessPoolExecutor`, so `--threads` spawns **processes**,
not threads. `/usr/bin/time -v` reports `getrusage(RUSAGE_CHILDREN)`, which is the
**maximum over children**, not their sum — it therefore equals the true footprint
only when one worker is running:

| workers | job | `/usr/bin/time -v` | `sacct MaxRSS` |
|---|---|--:|--:|
| 1 | 5982392 | 13.98 GB | 13.32 GB ← agree |
| 2 | 5983793 | 12.99 GB | **21.86 GB** |
| 4 | 5983772 | 12.99 GB | **35.29 GB** |

The published 12.99 GB is a `/usr/bin/time` reading, i.e. **one worker's share of a
two-worker run**, and it reproduces to 0.05 % when re-measured the same way
(13,622,540 kB published vs 13,620,700 kB re-measured). Nothing failed to reproduce.

**But it is the wrong number to publish.** `sacct MaxRSS` is the concurrent cgroup
total, which is what a user must provision, and it is what makes the manuscript's own
claim ("grows with thread count") true — 12.99 GB is flat in the number of workers.
Competing tools are single-process and multi-threaded, so for them the two metrics
coincide; comparing their figure against BWTandem's `sacct` total is the fair
direction.

## 3. Why the re-measured runtime is larger

| | partition | nodes | generation | workers |
|---|---|---|---|--:|
| BWTandem published | `cpu-s2-core-0` | cpu-64 / cpu-78 | **intelv5** | 2 |
| BWTandem re-measured | `cpu-s1-pgl-0` | cpu-48–51 | **intelv4** | 2 |
| all competitors | `cpu-s1-pgl-0` | cpu-48–51 | **intelv4** | 1–2 |

Same input, same worker count, same `/usr/bin/time`, 196 % vs 197 % CPU: 7:49:08 on
intelv5 against 12:05:37 on intelv4. **The 1.55× is node generation.** The
re-measured runs are the only BWTandem numbers taken on the same hardware as the
competitors.

## What is quarantined

**Do not reuse** the published runtime cells — 7.4 h (human), 6.6 h (maize),
0.31 h (Col-CEN). They belong to the Jul 16 array, whose BEDs no longer exist.

**Do not reuse** the published memory cells either — 12.99 / 14.37 / 1.31 GB. They
reproduce exactly, but they are per-worker readings and understate the footprint.

Both are superseded by jobs 5983792/93/94, where BED, runtime and memory come from
one recorded execution with the environment printed in the log:

| genome | published (quarantined) | in use now | why |
|---|---|---|---|
| Col-CEN | 0.31 h / 1.31 GB | **0.51 h / 1.95 GB** | wrong array / per-worker |
| human | 7.4 h / 12.99 GB | **12.1 h / 21.86 GB** | wrong array / per-worker |
| maize | 6.6 h / 14.37 GB | **15.4 h / 22.41 GB** | wrong array / per-worker |

`CLAUDE.md_stale_figures.md` in this directory holds the block removed from
CLAUDE.md; its accuracy figures are fine, its cost figures are not.

**Not quarantined:** every accuracy figure. Verified by re-measurement on all three
genomes.

## 4. A scoring bug of mine, which made the defect look far worse than it was

`period_of()` in `score_colcen.py`, `score_maize_3a.py` and
`analyze_unique_regions.py` read the period as `int(float(parts[4]))`. filip's
converted BEDs are 6 columns with an integer period in column 5; a BED written by
BWTandem itself is 8 columns whose column 5 is a **copy count** like `46.3`, with no
period column at all (there the period is `len(motif)`). `int(float("46.3"))`
silently returns 46, so every period-banded metric counted copy numbers.

That is what produced the alarming "CEN180 count collapses 1,380 -> 21" and
"maize SSR bp 41 M -> 131 M". With `int()` instead, Col-CEN re-measurement returns
**1,380 exactly** and maize Table 3A returns 41,314,898 bp against a published
41,007,584. Fixed 2026-08-05; the other eight scorers were already correct.
**If you write a new scorer, use `int()`.**

## Also affected by the scoring bug, and re-run

- WP3.1 unique-region **period distribution** (length, entropy and catalog-overlap
  figures never used the period and are unaffected; the "not low-complexity"
  conclusion stands). Re-run: job in `logs/wp31_refix_*.out`.
- The gap-fill ablation's banded columns. Corrected values: base 1,294 CEN180 calls
  and 91.72 % banded recall, off 991 / 85.38 %, away 993 / 85.57 %. The ablation's
  conclusion is unchanged — disabling the pass costs ~5.8 pp coverage, moving the
  window outside CEN180 costs 0.54 pp.
- Maize **Table 3A-b** carried the pre-re-measurement BWTandem row until 2026-08-05:
  41,007,584 bp / 1,498,365 regions / 35,017 TAG arrays, now **41,314,898 /
  1,508,821 / 34,247**. Table 3A had been updated and 3A-b had not, so the same
  quantity appeared twice with different values.
