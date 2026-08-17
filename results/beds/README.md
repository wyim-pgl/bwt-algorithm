# Deposited BWTandem interval files

The three whole-genome BWTandem outputs behind Tables 1a, 1c, 2 and 3, gzipped.
Every accuracy figure attributed to BWTandem in those tables is computed from one
of these files, so a reader can re-run the scoring scripts in `../../scripts/scoring/`
without access to our cluster.

| file | genome | calls | assembly | configuration |
|---|---|--:|---|---|
| `bwtandem_human.bed.gz` | human GRCh38 (UCSC hg38, chr1–22, X, Y) | 4,009,261 | 3.1 Gb | v2.2 gate + catch-all id 0.72, `MIN_COPIES=3`, 1–2,000 bp, 2 threads |
| `bwtandem_maize.bed.gz` | *Zea mays* Mo17 T2T (GCA_022117705.1) | 2,405,401 | 2.18 Gb | same as human |
| `bwtandem_colcen.bed.gz` | *A. thaliana* Col-CEN v1.2 | 161,187 | 132 Mb | v2.2 gate, catch-all **off**, 1–2,000 bp, 2 threads |

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
call (`1`, `2`, `3` for the three tiers; the deposited files carry the text
label `satellite` for gap-fill calls, which the released code renders
numerically — see Data and Code Availability in the manuscript).

## Verifying

```bash
sha256sum -c SHA256SUMS
gzip -dc bwtandem_colcen.bed.gz | head
```

## Provenance caveat

These files were produced by the working tree described in the manuscript's Data
and Code Availability section (repository HEAD `0e17d1a` plus uncommitted `src/`
modifications), not by the released commit `294f8ac`. Re-running the released
code on the human genome yields 198 fewer calls, all of them motifs containing
characters other than A, C, G and T. A re-run on the released build is in
progress; when it lands it will be deposited alongside these files rather than
replacing them, so the numbers in the published tables stay traceable to the
files that produced them.
