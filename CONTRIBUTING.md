# Contributing

## Development setup

```bash
pip install -e .[dev]        # or: conda env create -f environment.yml
# build the Cython accelerator (required for a real run; see CLAUDE.md)
python -c "
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
ext = [Extension('src._accelerators', ['src/_accelerators.pyx'], include_dirs=[np.get_include()])]
setup(script_args=['build_ext', '--inplace'], ext_modules=cythonize(ext, compiler_directives={'language_level': '3'}))
"
python -m pytest tests/ -q   # run from the repo root
```

CI runs the native build, the full suite, a pure-Python fallback parity leg,
a CLI smoke test, and a clean-environment wheel-install check. Pull requests
must pass all of them, and code changes need accompanying tests — the parity
tests (`test_accel_parity.py`, `test_align_parity.py`) pin the C and Python
paths to identical output, so a change to one side usually means a change to
both.

Commit messages follow Conventional Commits (`feat:`, `fix:`, `docs:`,
`test:`, `chore:`, `perf:`, `ci:`).

## Repository hygiene and artifact policy

A fresh clone is expected to have a **clean `git status`** and to pass the
test suite after the accelerator build above. What lives where:

| Class | Destination | Examples |
|---|---|---|
| Source, tests, docs | Git | `src/`, `tests/`, `docs/`, `scripts/` |
| Benchmark provenance | Git, under `results/` | `results/manifest.tsv` (+ `manifest.sha256`) maps every reported table row to the interval file, scoring script, hash, and commit that produced it |
| Deposited evidence | Git, under `results/` | the three whole-genome BWTandem BEDs (`results/beds/`), curated ground-truth coordinate sets (`results/ground_truth/`), audit sheets (`results/audit11/`) |
| Superseded / quarantined claims | Git, under `archive/` | `archive/2026-08-05-unreproducible/` documents which historical cost figures were re-measured and why — read it before reusing any cost number |
| Small test data | Git | `arabadopsis_chrs/`, `tests/fixtures/`, `examples/` |
| Genomes, comparator outputs, run logs | **Not in Git** — external benchmark workspace | referenced by absolute path + SHA-256 in `results/manifest.tsv`; the manifest hash is the ground truth when a path is unreachable |
| Generated output | Not in Git (ignored) | builds (`build/`, `dist/`, `*.so`), caches (`__pycache__/`), run output (`*.dat`, `slurm-*.out`), session working notes (`resume.md`) |

Rules of thumb:

- Never track generated files; if a new generated pattern shows up, extend
  `.gitignore` in the same change that introduces it.
- Anything a reported number depends on either lives in `results/` with a
  hash, or is listed in `results/manifest.tsv` with a hash. No exceptions —
  a figure whose input cannot be hashed is not reportable.
- Retiring a claim is done by moving its write-up under `archive/` (kept
  tracked, clearly marked superseded), not by deleting it.
- Large artifacts (> a few MB) that must be preserved belong in the external
  workspace or an archival deposit (e.g. Zenodo), never in Git history.
