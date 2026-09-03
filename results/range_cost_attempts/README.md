# Matched-ceiling competitor attempts (range-cost, human GRCh38)

BWTandem searches 1–2,000 bp on every genome. Two post-hoc attempts were made
to run the competitors at that same ceiling on human, so that the range-cost
comparison of Results 3.1 would have a competitor arm. Neither completed. Both
are recorded in `results/manifest.tsv` under `table = range-cost`, in Methods
2.2 and Results 3.1 of the manuscript, and in Supplementary Table S1.

| Attempt | Job | Threads | Elapsed at termination | Progress at termination | Manifest row |
|---|---|---|---|---|---|
| TRF 4.10.0rc2, `MAXP 2000` | 6076847 | 1 | 6 d 13 h 57 m (6.6 d; 4.7× its 33.7 h 500 bp run) | partial `-ngs` output, 379,077 lines | `TRF-p2000-attempt` |
| ULTRA 1.2.1, `-t 2 -p 2000` | 6145581 | 2 | 1 d 22 h 15 m (1.55× its 29.8 h 100 bp run, which covered the whole file) | 138,425 calls, all on chr1 (NC_000001.11) up to 124,786,615 bp, ~4% of the 3.25 Gb file; output did not grow during the final 5 h | `ULTRA-p2000-attempt` |

Both ratios in that table divide a **terminated partial** run by a **completed** one,
and both denominators (33.7 h, 29.8 h) are inherited GNU-time figures whose job
records no longer exist (manuscript Section 2.2.1). They bound the cost from below;
they are not like-for-like speed ratios.

## `ultra_p2000/`

- `run_ultra_human_p2000.sbatch`: the submitted script. It matches the published
  human ULTRA invocation (`ultra -t 2 -o OUT.tsv FASTA`, ULTRA 1.2.1, the same
  GCA_000001405.15 FASTA) except for `-p 2000`. **One difference beyond the
  period, verified rather than assumed:** this run used a local installation of
  the binary, whereas the published run ran inside the Singularity sandbox. The
  two are version- and input-matched, not environment-matched.

  The published run's GNU-time log survives at
  `/data/gpfs/assoc/pgl/filip/bwtandem_results/benchmarking_results/ultra/logs/GCA_000001405.15_GRCh38_genomic_run.log`
  and records `singularity exec ./bwtbench.sif ultra -t 2 -o …
  GCA_000001405.15_GRCh38_genomic.fna` with an elapsed time of **29:46:49**,
  which is the 29.8 h the manuscript quotes. (A second human ULTRA log, for the
  RefSeq-flavoured `GCF_000001405.26`, records 27:24:19; the manifest and the
  tables use the GCA run.) `bwtbench.sif` is a writable sandbox **directory**,
  not an image file, so its binary can be compared directly:

  | | path | size | SHA-256 of first MB |
  |---|---|--:|---|
  | published | `…/filip/bwt/bwtbench.sif/opt/ULTRA/ultra` | 638,912 | `f11b614f2858f9f0…` |
  | this attempt | `~/micromamba/bin/ultra` | 573,672 | `a475843c1879aeb5…` |

  Different files, both self-reporting version 1.2.1. An earlier draft of this
  README and of Methods 2.2 called them "the same binary"; that was wrong and was
  corrected in `24cd12a`, and the table above is the evidence.
- `ultra_h_p2000_6145581.log`: provenance header plus an hourly
  `PROGRESS ... out_bytes=` marker (output-file size; ULTRA writes through a
  4 KB stdio buffer, so the size is a coarse progress proxy). Growth was steady
  at roughly 100–470 KB per hour for the first 41 h, then zero from 10:05 until
  the cancellation at 15:20 on 2026-09-02. The last emitted calls
  (124.74–124.79 Mb, periods 169–340 bp) place the run inside the chromosome 1
  centromeric alpha-satellite region.
- `ultra_human_p2000.tsv.settings`: ULTRA's own parameter dump for the run
  (`max_period: 2000`, `threads: 2`, everything else default).

The partial output itself (`ultra_human_p2000.tsv`, 9,781,248 bytes, 138,426
lines including the header, first-megabyte SHA-256 `c30692f4351e9a0d…`) stays
on the cluster at the path recorded in the manifest and is hashed in
`results/external_evidence.sha256`; it is not scored. Note that the manifest's
`lines` field for this row is the raw line count, 138,426, one more than the
138,425 calls, because this row points at ULTRA's own TSV rather than at a
converted BED as every other competitor row does.

SLURM accounting (`sacct -j 6145581`): CANCELLED by the user at
2026-09-02T15:20:13, Elapsed 1-22:15:07, batch-step MaxRSS 17,972,740 K
(17.14 GB by the manifest's K/1024² convention; the published 100 bp run
recorded 1.68 GB under GNU time). The cancellation was a decision, not a
failure: at the pre-stall rate of about 3 Mb per hour the 3.2 Gb assembly
extrapolates to well over a month, beyond the 14-day partition limit, and the
run had produced no new output for five hours.

No artefacts of the TRF attempt survive beyond its partial `-ngs` file, which
is likewise on the cluster path recorded in the manifest.
