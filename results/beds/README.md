# Deposited BWTandem interval files

The three whole-genome BWTandem outputs behind Tables 1a, 1c, 2 and 3, gzipped.
Every accuracy figure attributed to BWTandem in those tables is computed from one
of these files, so a reader can re-run the scoring scripts in `../../scripts/scoring/`
without access to our cluster.

| file | genome | calls | assembly | configuration |
|---|---|--:|---|---|
| `bwtandem_human.bed.gz` | human GRCh38 (UCSC hg38, chr1–22, X, Y) | 4,014,108 | 3.1 Gb | relaxed short-period gate + catch-all id 0.72, `MIN_COPIES=3`, 1–2,000 bp, 2 threads |
| `bwtandem_maize.bed.gz` | *Zea mays* Mo17 T2T (GCA_022117705.1) | 2,406,800 | 2.18 Gb | same as human |
| `bwtandem_colcen.bed.gz` | *A. thaliana* Col-CEN v1.2 | 161,330 | 132 Mb | relaxed short-period gate, catch-all **off**, 1–2,000 bp, 2 threads |

The exact environment for each configuration is Supplementary Methods S2 of the
manuscript; `../manifest.tsv` carries the job identifier, elapsed time and cgroup
peak for each run, and now points at these files by repository-relative path.

## Format

BED, tab-separated, one call per line, 0-based half-open coordinates:

```
chrom  start  end  motif  copies  tier  mismatch_rate  strand
```

Column 5 is a copy count with one decimal, not a period — the period is the
length of the motif in column 4. Column 6 records the pass that produced the
call numerically (`1`–`5`).

## Verifying

```bash
sha256sum -c SHA256SUMS
gzip -dc bwtandem_colcen.bed.gz | head
```

## Provenance

All three files were regenerated from a clean checkout of commit
`0363d8bdc83cb8fb8e62c19d665479ce878e9f70`. The run-specific JSON provenance
records and scoring outputs are named in `../manifest.tsv`.
