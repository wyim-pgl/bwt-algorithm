# Comparator baselines — fixed historical outputs and their search ranges (#16)

**Policy (adopted 2026-08-27):** the regeneration campaign reruns BWTandem
only. Competitor outputs are retained as **fixed historical baselines**; they
are not silently comparable at equal search ranges, so every table that uses
them must disclose the ranges below. Where a range-corrected rerun exists it
is listed beside the published run and the tables should prefer it. The one
conspicuous gap that cannot be closed is human TRF at 2,000 bp: its rerun was
cancelled incomplete after 6.6 days (manifest row 66), 4.7x the 500 bp
runtime.

Sources: `results/manifest.tsv` (rows named), `manuscript.md` §2.2/§4.

**Where the inherited competitor cost figures come from.** The earlier
benchmarking round's GNU-time logs survive on the cluster, one per tool and
assembly, under
`/data/gpfs/assoc/pgl/filip/bwtandem_results/benchmarking_results/<tool>/logs/`.
Each records the full `singularity exec ./bwtbench.sif …` command line, the
elapsed wall clock and the maximum resident set size, so every competitor
runtime and memory cell can be checked against a preserved file even though the
SLURM accounting entries for those jobs are gone. Verified for human
GCA_000001405.15 on 2026-09-02:

| tool | elapsed | max RSS (kB) | as reported |
|---|---|--:|---|
| ULTRA | 29:46:49 | 1,758,540 | 29.8 h, 1.68 GB |
| TRF | 33:43:46 | 1,518,668 | 33.7 h, 1.45 GB |
| tantan | 51:46.94 | 281,316 | 0.9 h, 0.27 GB |
| TRASH | 107:37:36 | 15,298,244 | 107.6 h, 14.59 GB |
| mreps | 54:41.08 | 6,691,220 | not reported — the human run is the chr4-only defect |

Every cost cell Table 1a prints for a competitor reproduces from these logs to the
precision the table uses; none was found to disagree.

A second log exists per tool for the RefSeq-flavoured `GCF_000001405.26`
assembly; the manuscript tables use the GCA runs.

## Human GRCh38 (Tables 1a/1b/1c) — BWTandem searches 1–2,000 bp

| Tool | Period range | Provenance | Note |
|---|---|---|---|
| TRF | max 500 (no min-period parameter) | published | 2,000 bp rerun infeasible: cancelled at 6.6 d incomplete |
| ULTRA | 1–100 (default) | published | default `-p 100`; 20x below BWTandem's max; 2,000 bp rerun infeasible: cancelled at 1 d 22 h with ~4% of the assembly done (manifest `ULTRA-p2000-attempt`) |
| tantan | window 100 (default `-w`) | published | structurally zero in the 101–2,000 bp stratum |
| tantan | 1–2,000 (re-run) | job 6085144 | range-corrected arm; prefer for Table 1c |
| mreps | unbounded | published | |
| TRASH | n/a (motif discovery) | published | |

## Arabidopsis Col-CEN (Table 2) — CEN180 band 150–200 bp

| Tool | Period range | Provenance | Note |
|---|---|---|---|
| TRF | max 200 (no min-period parameter) | published | |
| ULTRA | 1–100 published; **1–500 re-run** | published / job 5981977 | published run cannot see CEN180 at all |
| tantan | window 100 published; **1–500 re-run** | published / job 6085141 | 6,075 calls land in 150–200 band vs 0 before |
| mreps | **150–400 re-run** | job 5981987 | |
| NCRF, TRASH-dn/tpl | motif-guided / n/a | published | |

## Zea mays Mo17 (Tables 3A/3B/3C)

| Tool | Period range | Provenance | Note |
|---|---|---|---|
| TRF | 6 (Table 3A STR band) | published | |
| ULTRA | 6 (3A) | published (recovered) | |
| tantan | window 100 default | published | enters 3B/3C only through boundary calls at period 100; the 3A/3B/3C published files are byte-identical (manifest rows 33, 39, 56) |
| NCRF | motif-guided | published | |
| TRASH-dn/tpl | n/a | published | experiment/template-specific BEDs; see manifest rows 35, 40, 43, 54, 55 |

## Rules for table construction after regeneration

1. Regenerated BWTandem rows are labelled with their provenance JSON
   (job id + commit `0363d8b`); competitor rows keep their historical labels.
2. A claim comparing counts across tools MUST name both ranges or use a
   matched band (e.g. 150–200, 1–100).
3. Range-corrected reruns replace published rows where they exist
   (tantan human/colcen, ULTRA colcen, mreps colcen); the published row stays
   in the supplement for continuity.
