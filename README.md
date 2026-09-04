# BWTandem — BWT-based Tandem Repeat Finder

De novo tandem repeat detection in assembled genomes (FASTA). One pass finds
everything from 1 bp microsatellites to 100 kb satellite units through a
three-tier Burrows-Wheeler Transform / FM-index pipeline, reported as BED,
VCF, TRF `.dat` or STRfinder `.csv`.

![How BWT finds tandem repeats](figures/Video%20Project.gif)

- **One wide-range pass** — periods 1–2,000 bp by default (search reaches
  100 kb via `--max-period`), so short STRs and long satellites come from a
  single run instead of per-range re-runs.
- **Three-tier pipeline** — exact FM-index enumeration for short perfect
  repeats, LCP/k-mer seeding with mismatch tolerance for medium imperfect
  repeats, sparse seeding with anchor verification for long arrays.
- **Four output formats** — BED (default), VCF, TRF-compatible `.dat`,
  STRfinder CSV.
- **Honest benchmarks** — every published number is provenance-tracked to a
  hashed input, script, and commit ([below](#benchmarks)).

## Installation

BWTandem is a standard Python package: `pip install .` compiles the Cython
accelerator into the wheel and installs a `bwtandem` command. Requirements:
Python ≥ 3.11 and a C compiler (`gcc`).

### Route A — micromamba / conda (recommended; brings its own compiler)

```bash
git clone https://github.com/wyim-pgl/bwt-algorithm && cd bwt-algorithm

# environment.core.lock.yml pins the exact core the published numbers were
# measured on (Python 3.11, numpy, Cython, gcc, pydivsufsort)
micromamba create -n bwtandem -f environment.core.lock.yml
micromamba run -n bwtandem pip install .

# verify (~30 s smoke test on bundled data)
micromamba run -n bwtandem bwtandem --help
micromamba run -n bwtandem bash examples/quickstart.sh
```

With conda instead of micromamba, the same file works:
`conda env create -f environment.core.lock.yml && conda activate bwtandem && pip install .`

### Route B — plain pip (any Linux/macOS with Python ≥ 3.11 and gcc)

```bash
git clone https://github.com/wyim-pgl/bwt-algorithm && cd bwt-algorithm
pip install .
bwtandem --help
```

`pip install .` pulls numpy and pydivsufsort automatically. On very old
systems (glibc < 2.24) the pydivsufsort binary wheel is not installable —
use Route A, which compiles it locally.

### What needs a compiler, exactly

- **At install time**: the Cython accelerator is compiled into the wheel, so
  `pip install .` requires `gcc`. This is the part that makes real runs fast.
- **At first run**: four small C libraries (`bwtandem/c_extensions/*.c`,
  shipped in the wheel) compile themselves on first import and need `gcc` on
  `PATH` at runtime — the Route A environment provides it. If they cannot
  build, faithful pure-Python fallbacks take over: **the same repeats are
  detected, just far more slowly** (`tests/test_accel_parity.py` pins the
  two paths to byte-identical output).

### Editable developer install

```bash
micromamba run -n bwtandem pip install -e .[dev]
micromamba run -n bwtandem python3 -m pytest tests/ -q   # 196 tests
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow.

### Docker

```bash
docker build -t bwtandem .
docker run --rm -v "$PWD:/data" bwtandem /data/input.fa --format bed -o /data/out
```

The image installs the pinned environment and pre-compiles every native
extension, so it also converts cleanly to a read-only Singularity/Apptainer
image (`singularity build bwtandem.sif docker-daemon://bwtandem:latest`).

## First run

From a repository checkout, on the bundled test data:

```bash
# 35 bp toy sequence: five copies of TCATCGG (runs in <1 s)
bwtandem arabadopsis_chrs/test_seq1.fa --format bed -o /tmp/quick
cat /tmp/quick.bed
# repeatTCATCGG_5    0    35    TCATCGG    5.0    1    0.000    +

# a real 367 kb sequence (Arabidopsis mitochondrial genome, ~2 s, 45 calls)
bwtandem arabadopsis_chrs/ChrM.fa --format bed -o /tmp/chrm -v

# or the packaged smoke test
bash examples/quickstart.sh
```

On your own data (any context):

```bash
# defaults: all tiers, periods 1-2,000 bp, BED written next to the input
bwtandem input.fa
```

`python3 -m bwtandem.main` is the module-form equivalent of the `bwtandem`
command (useful from an uninstalled checkout).

## What the output looks like

```
chr1    100    145    AT       22.5    1    0.022    +
chr1    500    620    AATGG    24.0    2    0.083    +
```

BED columns: `chrom start end motif copies tier mismatch_rate strand`
(0-based; the motif is the canonical rotation reduced to the primitive
period). Complete schemas for all four formats, with real examples, are in
[MANUAL.md §5](MANUAL.md#5-output-formats).

## Common options

| Option | Default | Meaning |
|---|---|---|
| `--min-period` / `--max-period` | 1 / 2000 | period range to search and report (Tier 3 reaches 100 kb if you raise the max) |
| `--tiers tier1,tier2,tier3` | all | run a subset, e.g. `--tiers tier1` for STRs only |
| `--format bed\|vcf\|trf\|strfinder` | bed | output format |
| `-t, --threads N` | 1 | one worker per sequence (helps multi-sequence FASTA only) |
| `--mask soft\|hard\|both` | none | skip soft-masked (lowercase) and/or `N` regions |
| `-o PREFIX` | input name | output prefix |

The full CLI reference, masking semantics, Tier 3 presets, and the
environment-variable sensitivity levers are in [MANUAL.md](MANUAL.md).

## How it works

One FM-index is built per sequence, then three tiers run in order, each
skipping regions the previous tiers already claimed:

1. **Tier 1** (motifs 1–9 bp): exact FM-index backward-search enumeration of
   all canonical short motifs (sliding-window scan on sequences > 10 Mbp).
2. **Tier 2** (motifs ≥ 10 bp): LCP-array candidates plus BWT k-mer seeding,
   extended with mismatch/indel tolerance for imperfect repeats.
3. **Tier 3** (periods 100 bp–100 kb): sparse k-mer seeding with adaptive
   parameters; ultra-long arrays use anchor-based boundary verification
   instead of full alignment.

Details, including the post-processing and the optional sensitivity passes,
are in [MANUAL.md §4](MANUAL.md#4-the-three-tier-pipeline).

## Benchmarks

All manuscript benchmarks were regenerated from clean commit `0363d8b` with
full provenance (SLURM job, elapsed time, sacct peak memory, SHA-256 of every
output). **The authoritative numbers live in [`manuscript.md`](manuscript.md)
and [`results/manifest.tsv`](results/manifest.tsv)** — every reported cell is
linked there to its source BED, scoring script and hash. Deposited evidence:
[`results/regen/`](results/regen/) (score reports, per-run provenance),
[`results/beds/`](results/beds/) (the three whole-genome BWTandem BEDs),
[`results/audit11/`](results/audit11/) (blinded specificity audit),
[`results/one_to_one/`](results/one_to_one/) (strict one-to-one scoring),
[`results/figures/`](results/figures/) (operating-point figure and data).

Headline results (human GRCh38 vs the GIAB adotto catalog, Arabidopsis
Col-CEN, maize Mo17; see the manuscript for every caveat):

| Claim | Measurement |
|---|---|
| One wide-range pass | periods 1–2,000 bp in 12.6 h / 25.3 core-hours on GRCh38 (ULTRA: 29.8 h / 59.6 core-hours capped at 100 bp; not range-matched) |
| Near-flat range cost | widening the reportable maximum 100→2,000 bp costs 1.30–1.41× in paired runs; we did not obtain completed 2,000 bp competitor runs on human (we stopped the TRF and ULTRA attempts terminated incomplete at 6.6 d and 1 d 22 h) |
| Shared-range accuracy | non-leading: ULTRA ranks first in region recall (81.62%); BWTandem 78.87% at the whole-genome configuration, restricted the same way its competitors are, 81.60% at a permissive setting with lower precision |
| Long-period stratum | 3.43% of catalog regions above period 100 (ULTRA: none by construction) |
| Plant satellites | Col-CEN CEN180 monomer recall 99.72% in 40 min; maize unfiltered coverage in the leading group |
| Specificity audit | 4 of 400 blinded BWTandem-only calls supported (single reader) — unmatched calls are predominantly over-calls |
| Cost of the design | 28.08 GB peak memory on human (ULTRA: 1.68 GB) and fragmented satellite calls |

Reproducing any number requires the environment overrides in Supplementary
Methods S2 of the manuscript — the benchmarked operating points differ from
the built-in defaults ([MANUAL.md §9](MANUAL.md#9-reproducing-the-manuscript-numbers-checkout-only)).

## Documentation

- [MANUAL.md](MANUAL.md) — the reference: CLI, masking, tiers, output
  schemas, environment variables, testing, reproduction, troubleshooting
- [CONTRIBUTING.md](CONTRIBUTING.md) — development setup, CI expectations,
  repository artifact policy
- [CLAUDE.md](CLAUDE.md) — architecture notes and the complete sensitivity
  env-var list (kept in sync with the code by `tests/test_env_var_docs.py`)
- [manuscript.md](manuscript.md) — the paper, with Supplementary Methods S2
  giving the exact configuration of every benchmarked run
- [CHANGELOG.md](CHANGELOG.md) — release history

## Citation

If you use BWTandem, please cite the repository (see
[CITATION.cff](CITATION.cff)) until the paper is published; the citation
will be updated on acceptance.

## License

MIT — see [LICENSE](LICENSE).
