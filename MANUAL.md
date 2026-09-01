# BWTandem Manual

BWTandem is a de novo tandem repeat finder for assembled genomes (FASTA). It
detects repeats with periods from 1 bp to 100 kb in a single pass through a
three-tier BWT/FM-index pipeline, and reports them as BED, VCF, TRF `.dat`
or STRfinder `.csv`.

This manual covers the command line, the output formats, configuration via
environment variables, and how to reproduce the manuscript numbers.
Installation lives in the [README](README.md#installation).

---

## 1. Command line

```
python3 -m bwtandem.main FASTA [options]
```

| Option | Default | Meaning |
|---|---|---|
| `fasta_file` | — | input FASTA (multi-sequence allowed; each sequence indexed and scanned independently) |
| `--min-period N` | 1 | minimum repeat period (bp) |
| `--max-period N` | 2000 | maximum reported period (bp). Tier 3's internal search window is fixed at 100 bp–100 kb regardless; this bounds what is reported |
| `--min-array-bp N` | none | drop arrays shorter than N bp |
| `--max-array-bp N` | none | drop arrays longer than N bp |
| `--tiers LIST` | all | comma list of `tier1,tier2,tier3`, or `all`. Invalid names are rejected (fail-closed) |
| `-o, --output PREFIX` | input name | output file prefix |
| `--format {bed,vcf,trf,strfinder}` | bed | output format |
| `-t, --threads N` | 1 | parallel processes, one sequence per worker. A worker failure aborts the whole run before partial output is written (fail-closed) |
| `--mask {none,soft,hard,both}` | none | skip lowercase (soft), `N` (hard), or both |
| `--tier3-mode {fast,balanced,sensitive}` | balanced | Tier 3 speed/accuracy preset |
| `-v, --verbose` | off | progress per sequence and tier |
| `--profile` | off | cProfile hotspots; also writes `PREFIX.profile.prof` |

Memory grows with thread count (one FM-index per in-flight sequence): on
human, roughly a 5 GB base plus ~8 GB per thread; a 250 Mb chromosome needs
~6 bytes per base for the index alone.

## 2. Output formats

Coordinates are 0-based internally; VCF and STRfinder output is converted to
1-based. The motif written is the canonical rotation: lexicographically
smallest among all cyclic rotations of both strands, reduced to the primitive
period (ATAT → AT, exact and ≤2%-error approximate tests).

### BED (default), one line per array

```
chrom  start  end  motif  copies  tier  mismatch_rate  strand
```

- `copies` — array length divided by period, one decimal.
- `tier` — 1, 2 or 3 for the detecting tier; 4 marks satellite gap-fill
  calls and 5 marks catch-all periodicity calls when those passes run.
- `mismatch_rate` — fraction of mismatching bases against the consensus.

### VCF

One record per array with `INFO` fields `END`, `MOTIF`, `COPIES`, `TIER`,
`MISMATCH`; `POS` is 1-based.

### TRF `.dat`

TRF-compatible table (1-based start) so downstream TRF parsers work
unchanged.

### STRfinder `.csv`

1-based CSV with the same fields as the BED plus the repeat sequence.

## 3. Configuration by environment variable

Two capabilities are env-var only:

- **`BWT_INDEX_CACHE=DIR`** — cache the FM-index across runs of the same
  sequence (keyed by a sequence hash, re-verified on load, atomic writes).
  Worth it for parameter sweeps, not single runs.
- **Sensitivity levers** — every tier threshold has an override
  (`TIER1_MIN_ARRAY_LEN`, `TIER1_MIN_SCORE`, `TIER1_FMSCAN`,
  `CATCHALL_SCAN`, `CATCHALL_MIN_IDENTITY`, `SAT_FILL_*`, …). Unset means
  the built-in default; **all defaults reproduce the baseline exactly**.
  The complete list with defaults and measured effects is in
  [CLAUDE.md](CLAUDE.md#tuning-detection-sensitivity-env-vars), and
  `tests/test_env_var_docs.py` fails if that list and the code diverge.

The optional passes matter for sensitivity work: `TIER1_FMSCAN=1` adds an
FM-enumeration pass that recovers diverged short STRs, and `CATCHALL_SCAN=1`
adds a periodicity pass that trades base-pair precision for region recall
(the manuscript's operating points P/B/F/H differ in exactly these levers).

## 4. Reproducing the manuscript numbers

Every benchmarked run is an environment-override configuration on top of the
defaults; the complete per-run blocks are in **Supplementary Methods S2** of
[manuscript.md](manuscript.md). The provenance chain:

- [`results/manifest.tsv`](results/manifest.tsv) — every reported table cell
  → source BED, SLURM job, elapsed, sacct peak, scorer + hash, commit.
- [`results/manifest.sha256`](results/manifest.sha256) — hashes of the
  deposited evidence set (`sha256sum -c results/manifest.sha256`).
- [`results/regen/`](results/regen/) — per-run provenance JSON and score
  reports; [`results/beds/`](results/beds/) — the three whole-genome BEDs.
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
python3 -m bwtandem.main hg38_primary.fa --min-period 1 --max-period 2000 \
    --threads 2 --format bed -o bwt_human -v
```

## 5. Tests

```bash
python3 -m pytest tests/ -q     # 179 tests
bash examples/quickstart.sh     # 30-second smoke test on bundled data
```

## 6. Known limitations

Stated plainly, as in the paper: shared-range accuracy on the human catalog
is non-leading (ULTRA ranks first in region recall under every overlap rule
tested); calls unique to BWTandem are predominantly unsupported under manual
audit, so recall-favouring settings need downstream filtering where
precision matters; peak memory is the largest of the compared de novo tools
on large genomes and grows with threads; satellite output is fragmented
(hundreds of calls per array before any coordinate merging), and reported
periods on satellites are less stable than competitors' under period-band
filtering.
