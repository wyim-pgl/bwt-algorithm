# 2026-08-10 — Flaky adjacent-test instrumentation campaign: every tool clean, reproduction blocked by a system-level ASLR change

Target: `tests/test_ground_truth.py::TestAdjacentGroundTruth::test_sensitivity`,
which failed ~60% of full-suite runs on 2026-08-09/10 (always the same two
arrays, `ACG@20075` and `GGT@45371`, 81.8% vs the 95% floor). The 2026-08-10
diagnosis established memory-layout dependence (one extra env var of any kind
suppresses it; `setarch -R` suppresses it) but did not locate the read. This
campaign attempted to finish the job with the tools that were missing then.

## What was run (all on build 9d0d669, all via sbatch unless noted)

| probe | scope | result |
|---|---|---|
| valgrind memcheck 3.27.1, `--track-origins=yes`, `PYTHONMALLOC=malloc`, production `.so` | `test_ground_truth.py` alone, then the full suite (115 tests) | **0 errors, 0 contexts** in both; 9.2M/3.19GB heap traffic confirms the pipeline really ran |
| ASan-instrumented Cython `_accelerators` (`-fsanitize=address -fsanitize-address-use-after-scope`, `detect_stack_use_after_return=1`), the one native component the 2026-08-10 ASan pass skipped | full suite ×3 | **0 reports** |
| bare reproduction: `env -i` on login-0 (live repo) and cpu-48 (worktree) | 3 + 15 runs | all pass |
| env-block-size sweep: `env -i PAD=<L>` for 34 lengths 0–4096 | 34 runs | all pass |
| dense env-size sweep: 256 sizes (0–2040, 8-byte step) × 5 reps, 16-task array, per-task isolated copies | **1,280 runs** | **all pass** |
| mmap-layout sweep: `ulimit -s` ∈ {1M…unlimited, 9 values} × 16 env sizes × 2 (ASLR=0 makes mmap base a function of the stack rlimit) | 288 runs | all pass |
| current-session full environment, plain `pytest tests/ -q` | 5 runs | all pass |

Total: **~1,640 full-suite executions, zero adjacent-test failures**, plus the
two instrumentation engines reporting nothing.

## The explanation: the cluster's ASLR state changed

Every host reachable today — login-0, login-1, cpu-48/49/50/51 — has
`/proc/sys/kernel/randomize_va_space = 0`. The original diagnosis contained a
`setarch -R` arm (0 failures in 4 runs) which is only a meaningful A/B if the
baseline had ASLR **on**; had it already been off, `setarch -R` is a no-op and
that arm would have failed at the same ~60% rate as bare. So ASLR was on when
the failures were observed and is off now: the entire cluster currently sits in
the diagnosis table's "ASLR off → 0 failures" row, which is exactly what
~1,640 green runs across three independent layout levers look like.

Consequences:

- The stochasticity of the original 8/13 bare failures came from per-run ASLR
  placement (stack/mmap/brk), not from anything the test controls. With ASLR
  off, layout is a deterministic function of (env block, argv, rlimits), and no
  combination we can reach deterministically lands in a failing configuration.
- Re-enabling randomization requires root (`randomize_va_space` is a sysctl;
  `setarch` can only disable, never enable). **Reproduction is impossible from
  user space on the current cluster.** The one path left is asking the admins
  whether/when `randomize_va_space` was changed and having it re-enabled on a
  single node for a reproduction window.

## What today's clean instrumentation does and does not prove

- memcheck clean (production binary, failing-context suite) rules out
  uninitialised heap/stack reads and heap overruns *on the code paths executed
  in today's layouts*. ASan on `_accelerators` additionally rules out global
  and stack overruns in the Cython extension in those layouts.
- It does **not** prove the bug gone: the offending read plausibly executes
  (or mis-resolves) only in layouts that ASLR-on produced. Those layouts are
  currently unreachable, so the instrumented runs never traverse the bad
  configuration.
- Publication impact remains nil on separate grounds: chromosome-scale output
  was byte-identical across three environment layouts when the bug *was*
  reproducible (Chr4, 25,766 calls, same sha256), and today's ~1,640 green
  suites are the strongest stability evidence yet for the current environment.

## Ready-made divergence tracer (for the day reproduction returns)

A stage-level tracer was built and smoke-tested for the bisection that was cut
short. Design constraints learned the hard way: it is enabled by a **marker
file** (`TRACE_ON` at the repo root), never an env var (an extra env var
changes the environment-block size and masks the bug), and it hooks nine
pipeline points in `finder.find_all` (tier1_raw, tier2_long_raw, tier2_med_raw,
tier3_raw, pre_merge, post_merge, post_filter, post_satellite, final), writing
one summary sha + per-repeat detail lines per stage to `traces/stagetrace_<pid>.tsv`.
Diffing a failing against a passing run pinpoints the first stage where
`ACG@20075` / `GGT@45371` vanish.

The module (`src/_stagetrace.py`) and the nine-line `finder.py` patch:

```python
"""Stage-level pipeline tracer for the flaky-adjacent-test divergence hunt."""
import hashlib
import os

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
ENABLED = os.path.exists(os.path.join(_ROOT, "TRACE_ON"))
_fh = None


def log(stage, chromosome, repeats):
    global _fh
    if not ENABLED:
        return
    if _fh is None:
        tdir = os.path.join(_ROOT, "traces")
        os.makedirs(tdir, exist_ok=True)
        _fh = open(os.path.join(tdir, f"stagetrace_{os.getpid()}.tsv"), "a")
    rows = sorted(
        (r.start, r.end, getattr(r, "motif", ""), round(getattr(r, "mismatch_rate", 0.0), 6))
        for r in repeats
    )
    sig = hashlib.sha1(repr(rows).encode()).hexdigest()[:16]
    _fh.write(f"S\t{chromosome}\t{stage}\t{len(rows)}\t{sig}\n")
    for row in rows:
        _fh.write(f"D\t{chromosome}\t{stage}\t{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")
    _fh.flush()
```

Hook placement (one `_stagetrace.log(...)` call after each of): the
`find_strs` return, the `find_long_unit_repeats_strict` return, the
`find_long_repeats` (tier 2) return, the tier-3 `find_long_repeats` return,
the post-sort `all_repeats`, `_merge_adjacent_repeats` return,
`_filter_overlaps` return, the satellite-pass loop exit, and the final
bounds-filtered list — plus `from . import _stagetrace` in the imports.

## Sundry facts recorded for the next attempt

- valgrind is now installed cluster-wide as a shared env:
  `/data/gpfs/assoc/pgl/bin/conda/conda_envs/valgrind/bin/valgrind` (3.27.1,
  conda-forge; catalogued in the lab micromamba README).
- The live repo's `_accelerators.so` was rebuilt 2026-08-10 13:17 — after the
  failures were observed — so today's binary is not bit-identical to the one
  that failed (same source, same compiler; relevant only if someone chases
  code-generation-level explanations).
- `pytest -p no:cacheprovider` changes plugin composition; the original
  failures were observed under plain `pytest tests/ -q`. Instrumented runs here
  used both forms; neither reproduced.
- Full-suite wall time is ~42 s native, ~70 s under ASan, ~7.5 min under
  memcheck; a 1,280-run sweep is an afternoon as a 16-task array.
