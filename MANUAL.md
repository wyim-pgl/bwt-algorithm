# BWTandem Manual

BWTandem is a de novo tandem repeat finder for assembled genomes (FASTA). It
detects repeats with periods from 1 bp up to 100 kb in a single pass through a
three-tier BWT/FM-index pipeline, and reports them as BED, VCF, TRF `.dat` or
STRfinder `.csv`.

This manual is the authoritative reference: command line, tier behaviour,
output schemas, environment-variable configuration, testing, and how to
reproduce the manuscript numbers. Installation lives in the
[README](README.md#installation).

---

## 1. Execution contexts

Three ways to run the same code — know which one you are in:

| Context | Command | Notes |
|---|---|---|
| Installed (pip/conda) | `bwtandem FASTA [options]` | the console entry point; also `python3 -m bwtandem.main` |
| Repository checkout | `python3 -m bwtandem.main FASTA [options]` | run from the repo root; needs numpy + pydivsufsort in the environment and a built accelerator for real workloads |
| Container | `docker run <image> FASTA [options]` | the `Dockerfile` image; `ENTRYPOINT` is `bwtandem` |

Checkout-only resources: the bundled test data (`arabadopsis_chrs/`,
`tests/fixtures/`), `examples/quickstart.sh`, the test suite, and everything
under `results/` and `scripts/`. The installed wheel ships only the
`bwtandem` package itself.

### Quick validation

```bash
# installed command, any FASTA
bwtandem your.fa --format bed -o check -v

# from a checkout: 30-second smoke test on bundled data
bash examples/quickstart.sh
```

## 2. Command-line reference

```
bwtandem FASTA [options]
```

| Option | Default | Meaning |
|---|---|---|
| `fasta_file` | — | input FASTA (multi-sequence allowed; each sequence indexed and scanned independently) |
| `--min-period N` | 1 | minimum repeat period (bp) |
| `--max-period N` | 2000 | maximum repeat period (bp). Bounds both search and reporting: Tier 3 searches up to `min(100000, max_period)`, so raise this to reach very long satellite units |
| `--min-array-bp N` | none | drop arrays shorter than N bp (period × copies) |
| `--max-array-bp N` | none | drop arrays longer than N bp |
| `--tiers LIST` | all | comma list of `tier1,tier2,tier3`, or `all`. Invalid names are rejected (fail-closed) |
| `-o, --output PREFIX` | input name | output file prefix; the extension follows the format |
| `--format {bed,vcf,trf,strfinder}` | bed | output format |
| `-t, --threads N` | 1 | parallel processes, one sequence per worker. A worker failure aborts the whole run before partial output is written (fail-closed) |
| `--mask {none,soft,hard,both}` | none | skip lowercase (soft), `N` (hard), or both — see §3 |
| `--tier3-mode {fast,balanced,sensitive}` | balanced | Tier 3 speed/accuracy preset — see §4 |
| `-v, --verbose` | off | progress per sequence and tier |
| `--profile` | off | cProfile hotspots; also writes `PREFIX.profile.prof` (inspect with `python -m pstats`) |

Memory grows with thread count (one FM-index per in-flight sequence): on
human, roughly a 5 GB base plus ~8 GB per thread; a 250 Mb chromosome needs
~6 bytes per base for the index alone. `--threads` parallelises across
sequences, so it has no effect on a single-sequence FASTA.

## 3. Input handling and masking

Sequences are uppercased internally; a `$` sentinel is appended for BWT
construction and excluded from detection. Coordinates are 0-based internally;
VCF and STRfinder output converts to 1-based.

| `--mask` | Lowercase (acgt) | `N` runs | Use case |
|---|---|---|---|
| `none` | uppercased and analysed | kept | default: analyse everything |
| `soft` | replaced with `N` (skipped) | kept | exclude RepeatMasker-style soft-masked interspersed repeats |
| `hard` | uppercased and analysed | skipped | exclude assembly gaps only |
| `both` | replaced with `N` | skipped | exclude both |

## 4. The three-tier pipeline

The coordinator (`TandemRepeatFinder`) builds one FM-index per sequence, then
runs the enabled tiers in order; each tier skips regions already claimed by
earlier tiers. Post-processing merges adjacent same-motif calls (gap ≤
max(10, motif length)), keeps the higher-scoring call among overlaps, and
applies the array-length bounds.

- **Tier 1 — short perfect repeats, motifs 1–9 bp** (≥3 copies). Sequences
  < 10 Mbp: FM-index backward search enumerates all canonical motifs;
  larger sequences use a sliding-window scanner with adaptive step.
- **Tier 2 — medium/imperfect repeats, motifs ≥ 10 bp.** A long-unit strict
  phase (LCP array via Kasai's algorithm, periods ≥ 20 bp) plus a general
  scanning phase (BWT k-mer seeding, periods 10–50 bp), both extended with
  mismatch tolerance (up to ~20% mismatches, ~10% indels).
- **Tier 3 — long repeats, periods 100 bp–100 kb** (bounded by
  `--max-period`). Sparse BWT k-mer seeding (k≈20, stride≈100), arithmetic-
  progression detection in occurrence positions. Ultra-long arrays (> 100
  copies or > 10 kb) use anchor-based boundary verification instead of full
  DP refinement.

Reported motifs are canonical: the lexicographically smallest rotation over
both strands, reduced to the primitive period (ATAT → AT; exact and
≤2%-error approximate tests).

### Tier 3 adaptive parameters and `--tier3-mode`

Tier 3 recomputes its search parameters per sequence from length, GC content
and prior-tier coverage: `kmer_size` 12–28, `stride` 20–300,
`allowed_mismatch_rate` 0.15–0.20, `tolerance_ratio` 0.02–0.04 (scales with
`--max-period`), `max_occurrences` 200–1500, anchor scan windows. The preset
shifts the speed/sensitivity balance: `fast` (larger k-mer/stride, for
screening big genomes), `balanced` (default), `sensitive` (smaller
k-mer/stride, for small sequences or fine mapping). Sequences over 100 Mbp
are forced to coarse settings (stride ≥ 150, k ≥ 20), sequences under 100 kb
to fine ones.

## 5. Output formats

All formats report the same calls; only the schema differs. Real output for
a toy 35 bp input (five copies of TCATCGG) is shown for each.

### BED (default) — `PREFIX.bed`, 0-based

```
repeatTCATCGG_5	0	35	TCATCGG	5.0	1	0.000	+
```

`chrom  start  end  motif  copies  tier  mismatch_rate  strand`

- `copies` — array length divided by period, one decimal.
- `tier` — 1, 2 or 3 for the detecting tier; 4 marks satellite gap-fill
  calls and 5 marks catch-all periodicity calls when those passes run.
- `mismatch_rate` — fraction of mismatching bases against the consensus.

### VCF — `PREFIX.vcf`, 1-based

```
##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the repeat">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
repeatTCATCGG_5	1	.	T	<STR>	.	.	END=35;MOTIF=TCATCGG;CONS_MOTIF=TCATCGG;COPIES=5.0;TIER=1;CONF=1.00;MM_RATE=0.000;MAX_MM_PER_COPY=0;N_COPIES_EVAL=5;STRAND=+
```

INFO fields emitted per record: `END`, `MOTIF` (canonical), `CONS_MOTIF`
(consensus), `COPIES`, `TIER`, `CONF` (confidence 0–1), `MM_RATE`,
`MAX_MM_PER_COPY`, `N_COPIES_EVAL`, `STRAND`.

### TRF `.dat` — `PREFIX.dat`, 1-based

```
1 35 7 5.0 7 100 0 70 14 29 29 29 1.95 TCATCGG TCATCGGTCATCGGTCATCGGTCATCGGTCATCGG
```

`Start End Period CopyNumber ConsensusSize PercentMatches PercentIndels
Score A C G T Entropy ConsensusPattern Sequence` — TRF-compatible, so
downstream TRF parsers work unchanged.

### STRfinder `.csv` — `PREFIX.csv`, 1-based, 12 columns

```csv
STR_marker,STR_position,STR_motif,STR_genotype_structure,STR_genotype,STR_core_seq,Allele_coverage,Alleles_ratio,Reads_Distribution,STR_depth,Full_seq,Variations
STR_repeatTCATCGG_5_0,repeatTCATCGG_5:1-35,[TCATCGG]n,"7[TCATCGG]5,0",5,TCATCGGTCATCGGTCATCGGTCATCGGTCATCGG,100%,-,5:5,5,TCATCGGTCATCGGTCATCGGTCATCGGTCATCGG,-
```

`STR_genotype_structure` is `period[MOTIF]complete_copies,truncated_bases`
(CSV-quoted because of the comma). Read-oriented columns
(`Allele_coverage`, `Reads_Distribution`, `STR_depth`) are filled with
assembly-appropriate placeholders since the input is a genome, not reads.

## 6. Environment-variable configuration

Two capabilities are env-var only; unset always means the built-in default,
and **all defaults reproduce the baseline behaviour exactly**.

### FM-index cache — `BWT_INDEX_CACHE=DIR`

The index depends only on the sequence, so a parameter sweep otherwise
rebuilds the identical suffix array once per configuration:

```bash
BWT_INDEX_CACHE=/scratch/idxcache bwtandem genome.fa ...
```

Keyed by a sequence hash (an edited or renamed FASTA cannot hit a stale
entry), re-verified on load, atomic writes; a cache that cannot be written
warns and does not fail the run. Size ≈ 6 bytes per base. Worth it for
sweeps, not single runs (construction is 0.4–7% of a whole-genome run).

### Sensitivity levers

Every tier threshold has an override (`TIER1_MIN_ARRAY_LEN`,
`TIER1_MIN_SCORE`, `TIER1_FMSCAN`, `CATCHALL_SCAN`, `CATCHALL_MIN_IDENTITY`,
`SAT_FILL_*`, …):

```bash
TIER1_MIN_ARRAY_LEN=20 TIER1_MIN_SCORE=20 bwtandem genome.fa
```

The dominant short-STR recall levers are `TIER1_MIN_ARRAY_LEN` and
`TIER1_MIN_SCORE`. `TIER1_FMSCAN=1` adds an FM-enumeration pass that
recovers diverged short STRs, and `CATCHALL_SCAN=1` adds a periodicity pass
that trades base-pair precision for region recall (the manuscript's
operating points differ in exactly these levers). The complete list with
defaults and measured effects is in
[CLAUDE.md](CLAUDE.md#tuning-detection-sensitivity-env-vars);
`tests/test_env_var_docs.py` fails if that list and the code diverge.

`BWT_DISABLE_NATIVE=1` forces the pure-Python fallback paths (byte-identical
output, much slower) — for debugging a suspected accelerator issue.

## 7. Worked examples

```bash
# STRs only, VCF output
bwtandem input.fa --tiers tier1 --format vcf -o strs

# periods 1-50 bp only
bwtandem input.fa --tiers tier1,tier2 --min-period 1 --max-period 50 -o short

# long repeats in a big genome, fast preset, arrays >= 500 bp
bwtandem genome.fa --tiers tier3 --tier3-mode fast --min-period 100 --min-array-bp 500 -o long

# whole genome, 4 workers, skip soft-masked interspersed repeats
bwtandem genome.fa -t 4 --mask soft -v

# where is the time going?
bwtandem input.fa --tiers tier1 --profile   # prints hotspots, writes input.profile.prof
```

## 8. Testing (checkout only)

```bash
python3 -m pytest tests/ -q     # 196 tests, ~1 min with the compiled accelerator
bash examples/quickstart.sh     # 30-second smoke test on bundled data
```

Without the compiled accelerator the cross-implementation parity tests skip
(there is nothing to compare against); everything else runs on the faithful
pure-Python fallbacks. Developer details are in
[CONTRIBUTING.md](CONTRIBUTING.md).

## 9. Reproducing the manuscript numbers (checkout only)

Requires a full repository checkout — the provenance evidence and scoring
code are not shipped in the wheel. Every benchmarked run is an
environment-override configuration on top of the defaults; the complete
per-run blocks are in **Supplementary Methods S2** of
[manuscript.md](manuscript.md). The provenance chain:

- [`results/manifest.tsv`](results/manifest.tsv) — every reported table cell
  → source BED, SLURM job, elapsed, sacct peak, scorer + hash, commit.
- [`results/manifest.sha256`](results/manifest.sha256) — hashes of the
  deposited evidence set (`sha256sum -c results/manifest.sha256`).
- [`results/regen/`](results/regen/) — per-run provenance JSON and score
  reports; [`results/beds/`](results/beds/) — the three whole-genome BEDs;
  [`results/one_to_one/`](results/one_to_one/) — strict one-to-one scoring.
- Scoring and audit code: [`scripts/scoring/`](scripts/scoring/);
  benchmark harness: [`scripts/benchmark/`](scripts/benchmark/)
  (`run_with_provenance.sh` refuses a dirty tree and records rusage +
  SHA-256 of every output).

Example — the whole-genome human configuration (operating point F):

```bash
export TIER1_FMSCAN=1 TIER1_FMSCAN_MIN_DENSITY=0.45 TIER1_FMSCAN_MIN_LLR=6.0
export TIER1_MIN_ARRAY_LEN=20 TIER1_MIN_SCORE=20 TIER1_MIN_COPIES=2
export TIER1_COPYBASE=6 TIER1_COPYADD=2 TIER1_EXT_COPIES=2
export TIER2_MISMATCH=0.30 TIER2_SHORT_REQ_COPIES=2
export TIER1_SHORT_PERIOD_MAX=9 TIER1_SHORT_MIN_ARRAY_LEN=17 TIER1_SHORT_MIN_SCORE=17
export CATCHALL_SCAN=1 CATCHALL_MIN_IDENTITY=0.72 CATCHALL_MIN_COPIES=3
bwtandem hg38_primary.fa --min-period 1 --max-period 2000 \
    --threads 2 --format bed -o bwt_human -v
```

## 10. Known limitations

Stated plainly, as in the paper: shared-range accuracy on the human catalog
is non-leading (ULTRA ranks first in region recall under every overlap rule
tested); calls unique to BWTandem are predominantly unsupported under manual
audit, so recall-favouring settings need downstream filtering where
precision matters; peak memory is the largest of the compared de novo tools
on large genomes and grows with threads; satellite output is fragmented
(hundreds of calls per array before any coordinate merging), and reported
periods on satellites are less stable than competitors' under period-band
filtering.

## 11. Troubleshooting

- **`bwtandem: command not found`** — the entry point installs with
  `pip install .`; from a bare checkout use `python3 -m bwtandem.main`.
- **Slow runs** — check that the accelerator compiled: `python3 -c "import
  bwtandem._accelerators"` must succeed, and the four
  `bwtandem/c_extensions/*.so` build on first import, which needs `gcc` on
  `PATH` at runtime. Without either, faithful pure-Python fallbacks run —
  same results, far slower.
- **`pydivsufsort` install fails** (very old glibc) — use the conda route in
  the README; the pinned environment compiles it locally.
- **Output location** — files are written as `PREFIX.<ext>` next to the
  input unless `-o` gives another prefix.
- **N-heavy assemblies** — ambiguous bases are never counted as periodic
  evidence; runs of `N` are skipped with `--mask hard`/`both`.
