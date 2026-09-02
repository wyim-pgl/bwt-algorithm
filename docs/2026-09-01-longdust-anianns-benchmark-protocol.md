# longdust + AniAnn's benchmark protocol (preregistered 2026-09-01)

Author decision: add the two 2026 tools to the benchmark (prerequisite for
the *Bioinformatics* venue path; see
`2026-09-01-venue-and-positioning-analysis.md`). Decisions below are fixed
BEFORE any run is scored, per this project's discipline.

## Tools and versions

- **longdust** v1.4-1-g9491215 (github.com/lh3/longdust, cloned
  2026-09-01), built with the bwtandem-env gcc against the micromamba
  zlib: `gcc -O3 -std=c99 -o longdust main.c longdust.c kalloc.c -lz -lm`.
  Binary + build line recorded in the run logs.
- **AniAnn's** 0.7.1 (github.com/marbl/anianns @ 0d79851, cloned
  2026-09-01), installed with `pip install .` into a venv on the
  bwtandem-env Python 3.11.

## Inputs — identical to the regenerated BWTandem runs

| Genome | FASTA |
|---|---|
| human | `/data/gpfs/assoc/pgl/devel/exp1_human/data/hg38_primary.fa` (3.1 GB) |
| maize | `/data/gpfs/assoc/pgl/filip/bwt/zmays/zmays.fna` (2.2 GB) |
| Col-CEN | `/data/gpfs/assoc/pgl/filip/bwt/colcen/colcen.fa` (134 MB) |

Unlike the 2024-era competitor baselines (which consumed the ~5% larger
accession-flavoured human FASTA), these runs consume the same inputs as
the regenerated BWTandem rows, so their cost numbers are directly
comparable to ours and their outputs carry chr names natively.

## Arms

| Arm | Command | Rationale |
|---|---|---|
| longdust default | `longdust FASTA` | recommended defaults |
| longdust long-motif | `longdust -k8 -w20000 FASTA` | the README's own recommendation "for longer motifs" — analogous to the tantan window re-runs |
| AniAnn's default | `anianns annotate -f FASTA -j 2` | recommended defaults, 2 threads (the human competitor-run thread convention) |

Each longdust arm runs on all three genomes. AniAnn's runs on all three
too: Col-CEN and maize are its target domain (satellites); the human run
is included so the cost row exists — if it does not complete within the
job limits that is reported as such, exactly as TRF's 2,000 bp attempt
was. `--classify` is NOT used (no satellite k-mer database): every metric
we score is coordinate-based, and class labels play no role in any table.

## What each tool can fairly be scored on

- **longdust** emits 3-column intervals with no repeat units (the paper's
  own stated limitation). Scoreable: human region/bp metrics (Table 1a
  pipeline), Col-CEN coverage/CEN180 monomer recall (unfiltered rule
  only — no period column exists), maize unfiltered coverage. Every
  period-conditional column is structurally empty for it, which is the
  structured-output point, not a scoring failure.
- **AniAnn's** emits array intervals with a periodicity score (monomer
  length) and no per-copy structure. Scoreable: the same interval-level
  metrics, plus period-banded rules using its periodicity column where a
  table's rule reads a period. Its own paper states the window design
  cannot detect arrays substantially smaller than the window (≥1,000 bp),
  so its absence from short-STR strata is expected and will be stated,
  not counted silently against it.
- Cost: `/usr/bin/time -v` per arm plus `sacct` for the single-arm jobs,
  same as the existing competitor rows.

## What this does NOT change

The deposited 2024-baseline rows and every existing table stay untouched;
longdust/AniAnn's enter as NEW rows clearly dated 2026 and scored by the
same frozen scorers (`score_table1_regen.py`, `score_colcen.py`,
`score_maize_*`). Related-work text will cite longdust's verbatim
"unable to report the repeat units" and AniAnn's window/database
limitations from their own texts.
