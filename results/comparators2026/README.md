# 2026 tool additions — longdust + AniAnn's (Section 2.2.4)

Preregistered protocol: `docs/2026-09-01-longdust-anianns-benchmark-protocol.md`
(fixed before any run was scored). These deposits feed the 2026 rows of
Tables 1a, 2, 3B and 3C.

## Versions and runs

| Tool | Version | Run job | Arms |
|---|---|---|---|
| longdust | 1.4 (git g9491215), built on an el7 node: conda gcc, `-O3 -std=c99 -I/usr/include -L/usr/lib64 -lz -lm`, binary sha256 `3469655a…` | 6146419 | default; `-k8 -w20000` × human/maize/Col-CEN |
| AniAnn's | 0.7.1 (git 0d79851), shared conda env `anianns` (bioconda pysam 0.24.0), no `--classify` DB | 6146422 (human), 6146483 (Col-CEN), 6146484 (maize) | defaults, `-j 2` |

Inputs are the chromosome-only FASTAs of the regenerated BWTandem runs, so
cost cells are comparable to the BWTandem rows and NOT to the 2024
competitor rows (which consumed the ~5% larger accession-flavoured FASTA).
Reproducer warning: AniAnn's exits 0 having annotated nothing when its
pysam `.fai` build fails on an unwritable FASTA directory (first attempt,
jobs 6146420/6146421) — symlink the FASTA beside a pre-built index.

## Costs (GNU time, per arm)

| Arm | Wall clock | Peak RSS |
|---|---|---|
| longdust human default / `-k8 -w20000` | 27:56 / 2:00:36 | 0.47 GB |
| longdust maize default / `-k8 -w20000` | 18:00 / 3:53:40 | 0.58 GB |
| longdust Col-CEN default / `-k8 -w20000` | 1:56 / 10:23 | 0.07 GB |
| AniAnn's human / maize / Col-CEN | 2:19:13 / 1:34:26 / 7:18 | 0.53 / 0.55 / 0.50 GB |

## Conversion and scoring

`scripts/scoring/convert_2026_tools.py` maps raw outputs onto the scorers'
period-column contract — longdust's converted BED is deliberately 3-column
(its paper: "longdust is unable to report the repeat units in case of
tandem repeats"), so every period-conditional metric is structurally empty
rather than fabricated; AniAnn's carries its array-level periodicity as
the period column and no motif. Scoring job **6146742** (repo `7b2113d`,
0 dirty; converted-BED hashes in the job log) loaded the FROZEN scorers
(`score_table1.py`, `score_colcen.py` via its LABEL:PATH interface,
`score_maize_postmerge.py`) and overrode only their source lists
(`scripts/scoring/score_2026_tools.py`), so the new rows come from exactly
the code that produced the published rows.

## Headline results (full outputs in the three .txt files here)

- **Human (adotto, full range, one-base rule)**: longdust 41.16% region
  recall / 57.16% precision (half of BWTandem's 80.53% recall at similar
  precision); `-k8 -w20000` 35.42/60.08; AniAnn's 0.10/83.79 — expected,
  its ≥1 kb window cannot see short arrays (its own stated limitation).
- **Col-CEN**: AniAnn's posts the table's highest centromere coverage
  (86.97%) from 60 array-level intervals; longdust matches the leaders'
  monomer recall (99.86%) with no repeat units (CEN180 count 0,
  structural).
- **Maize satellites**: longdust emits **zero calls inside any curated
  knob180/TR-1/CentC array** in either arm — these satellites are not
  low-complexity under its model. AniAnn's leads array-level coverage
  (knob180 87.61%, TR-1 67.98%, CentC 81.42% versus BWTandem's
  79.79/50.34/58.55) with boundary offsets of 14.9–83.7 kb versus
  BWTandem's 253 bp (knob180, gap 0).

The complementary profile these numbers draw — interval-only speed
(longdust), array-level coverage with coarse boundaries (AniAnn's),
per-copy structured records at base-pair boundaries (BWTandem) — is the
positioning stated in the Related Work.
