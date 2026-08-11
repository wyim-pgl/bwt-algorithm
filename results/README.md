# Provenance for the reported figures

`manifest.tsv` maps each table row in the manuscript to the interval file it was
scored from, the scoring script and its hash, the scoring rule, and the
repository commit. `figures/` holds Figure 1 and the values behind it.
`ground_truth/` holds the curated coordinate sets every accuracy figure is scored
against.

## ground_truth/

| File | What it is |
|---|---|
| `colcen_cen180.bed` | 66,683 CEN180 monomers on Col-CEN v1.2, columns 5–6 are percent identity and strand |
| `CEN178_consensus.fa` | the 177 bp query used to produce them |
| `colcen_centromeres.bed` | the five Arabidopsis centromere intervals |
| `mo17_knob180_arrays.bed`, `mo17_tr1_arrays.bed`, `mo17_centc_arrays.bed` | the 25, 17 and 17 curated maize arrays from Chen et al. (2023) |

The CEN180 coordinates are blastn hits of `CEN178_consensus.fa` against the
assembly, filtered from 68,840 raw hits to the 66,683 deposited here; the raw hit set itself was not retained. **The blastn
version and command line were not retained.** The deposited coordinates are
therefore the reproducible artefact, not the procedure that made them — every tool
in Table 2 is scored against this one file, so no comparison in the paper depends
on the lost invocation, but regenerating the set from scratch requires choosing
parameters afresh.

The human ground truth is not deposited here: it is adotto v1.2.1, obtainable
unmodified from https://zenodo.org/records/13987414 and used restricted to
chr1–22, X and Y.

## What the manifest does and does not cover

**Accuracy figures — fully covered.** Every recall, precision, coverage,
detection-count, offset and fragmentation number in the manuscript traces to a
row here: source file, byte count, line count, a SHA-256 of the first megabyte,
the scorer, the scorer's hash, and the rule applied. Re-running the named script
on the named file reproduces the cell.

**BWTandem cost figures — fully covered.** Runtime and peak memory for every
BWTandem row come from a single recorded SLURM execution. The `producer_job`,
`elapsed`, `peak_rss_gb_sacct` and `node` columns identify it. Peak memory is the
cgroup maximum from `sacct MaxRSS`, not GNU time, because BWTandem distributes
chromosomes across worker processes and GNU time reports the maximum over
children rather than their sum.

**Competitor cost figures — not covered, and cannot be.** Rows whose
`producer_job` is the literal string `published` come from an earlier
benchmarking round. The output files and the command lines survive; the job
accounting records do not. Those rows therefore have empty `threads`, `elapsed`,
`peak_rss_gb_sacct` and `node` fields, and the runtime and memory printed in
Tables 1a, 1d, 2, 3A, 3B and 3C for those tools are the GNU time values recorded at
the time. They are reported because they are the best record available, and they
should be read as approximate. The accuracy figures computed from the same files
do reproduce.

For scale on how approximate: two BWTandem executions of identical settings on
the same genome at the same thread count differ by 29% in wall clock
(9 h 21 m against 12 h 06 m, both 4,009,261 calls). Cost differences below that
magnitude between any two rows here are not interpretable.

## The `repo_commit` column

Every row reads `0e17d1a`, which was the repository HEAD when the runs were
executed. The working tree at that moment carried uncommitted modifications under
`src/`; those modifications are now **commit `294f8ac`**, and that is the commit to
check out to obtain the code behind these results.

Two behavioural differences separate `294f8ac` from the tree that actually ran.
First, it carries the fix that stops the satellite gap-filling pass from emitting
motifs containing characters other than A, C, G and T; re-running the released code
on the human genome therefore yields 198 fewer calls (4.57 Mb) than
`remeas_human.bed` contains. Second, the released code renders the column recording
each call's originating pass numerically, where the deposited outputs carry a text
label on satellite calls — a byte-level formatting change that alters no coordinate,
motif or statistic. Every other output should reproduce.

## Known issue: a memory-layout-dependent flaky test

`tests/test_ground_truth.py::TestAdjacentGroundTruth::test_sensitivity` fails in
roughly 60% of full-suite runs on a small synthetic fixture. It is memory-layout
dependent rather than a logic defect, and it does **not** reach chromosome-scale
output: an 18.8 Mb chromosome produces byte-identical results under three different
process memory layouts. No figure in this manifest depends on it. See CLAUDE.md for
the full evidence table.

## Row classes

| `table` value | What it is |
|---|---|
| `1a`, `1b`, `1c`, `1d`, `1e` | human GRCh38 against adotto v1.2.1 |
| `sweep` | the single-parameter catch-all identity sweep, full range |
| `range-cost` | the 4-thread 1–2,000 bp arm; paired with row `1b` BWTandem for the range-cost claim. Cost-only, not scored |
| `2` | Arabidopsis Col-CEN v1.2 |
| `ablation` | Col-CEN satellite gap-fill ablation, three arms from one job |
| `3A`, `3B`, `3C`, `3B-b/3C-b` | Zea mays Mo17 |

Rows whose `source_bed` reads `(same source BEDs as …)` are alternative scorings
of files already listed, not new outputs.
