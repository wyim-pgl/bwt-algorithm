# Adversarial review of untracked submission deliverables

Commit inspected: `0363d8bdc83cb8fb8e62c19d665479ce878e9f70`.

## `scripts/scoring/score_one_to_one.py` -- fixes required

1. **High, lines 64--90:** descending-overlap greedy is neither maximum-cardinality
   nor maximum-total-overlap bipartite matching. One long candidate can consume
   the only match of one truth while a slightly smaller edge would permit two
   matches. Consequently sensitivity/precision depend on a heuristic, and ties
   depend on reverse `(truth index, prediction index)` ordering. Fix: build each
   chromosome's eligible bipartite graph, use maximum-cardinality matching with
   overlap as the secondary objective (e.g. min-cost max-flow), and document the
   objective. At minimum rename every output/claim `greedy_1to1` and add a
   counterexample regression test.
2. **Medium, lines 69--82:** candidate generation is Cartesian within each
   chromosome, O(T*P) time and potentially O(T*P) memory. Fix with a sweep-line
   active set; only materialize overlapping edges.
3. **Medium, lines 112--119:** `--min-overlap` accepts negative, zero, NaN, or
   values >1. Negative/zero silently reduces the rule to any positive overlap;
   NaN also defeats both comparisons. Concrete guard after parsing:
   `if not math.isfinite(x) or not 0 < x <= 1: ap.error(...)`.
4. **Medium, lines 19--20 and 46--52:** the same column-5 parser is used for truth
   and prediction and is described as copies. That happens to fit current
   8-column BWTandem output but will misread comparator BEDs whose column 5 is
   period. Fix with explicit `--truth-period-col/--truth-copy-col` and
   `--pred-period-col/--pred-copy-col` (1-based), retaining motif length only as
   an explicit fallback.
5. **Low, lines 149--150:** truth copy count zero is skipped by truthiness, while
   negative values can be divided. Use `is not None` and require truth copies >0.

## `scripts/scoring/sample_specificity_audit.py` -- fixes required

1. **High, lines 96--100:** the claimed “4x100” sample silently becomes smaller
   when a stratum has <100 records. That changes the preregistered design. Fail
   closed unless an explicit `--allow-underfilled` is supplied.
2. **High, lines 85 and 89--92:** an existing output directory and sealed key are
   overwritten. Use `os.mkdir(outdir)` (must not exist), write temporary files
   with mode `x`, fsync, then rename; set the key to mode 0600.
3. **Medium, lines 32--57 and 87:** `.fai` offsets/line widths are byte offsets,
   but FASTA is opened in text mode. Newline translation and multibyte decoding
   make `seek()` non-portable; CRLF leaves `\r` because only `\n` is removed.
   Open `rb`, strip both `b"\n"` and `b"\r"`, decode ASCII, and assert the return
   length equals the clipped requested length.
4. **Medium, lines 45--47 and 103--109:** missing chromosome raises an unlabelled
   `KeyError`; malformed/inverted BED coordinates are not validated; clipped
   flanks are not recorded. Fail with source line/chromosome, require `0<=s<e`,
   and add `left_flank_bp/right_flank_bp` columns.
5. **Medium, lines 52--57:** the loop slices after removing linefeeds but does not
   reject premature EOF. A truncated FASTA can yield a short context that is
   silently accepted. The exact-length assertion above fixes this.
6. **Low, lines 108--109:** lowercase and non-ACGT/gap content are passed through
   without normalization or disclosure. Uppercase the sequence and add masked
   base count/fraction so reviewers do not mistake gaps for periodic evidence.
7. **Low/design, lines 91 and 106:** `period_bp` reveals the stratum, so the sheet
   is not blinded to period stratum (only tier/score/source are blinded). Either
   describe this accurately as source-blinded, or put proposed period in a
   separate randomized prompt used only during scoring.

## `scripts/benchmark/acceptance_gate.sh` -- mostly fixed, two gaps

The `note()` both-branches defect is confirmed fixed at line 24, and the script
fails closed for the principal SLURM/provenance/hash/schema mismatches.

1. **Medium, lines 54--58:** zero or multiple BED objects in `outputs` are not
   rejected; shell `read` silently takes the first line. Replace the Python with
   a check that `len(beds)==1`, printing one tab-separated record only then.
2. **Medium, lines 43--48 and 63--85:** a zero-call log plus an empty, correctly
   hashed BED passes count and schema checks. For these three known genomes add
   `[[ -s "$BED" ]] || note FAIL ...` and require schema `rows > 0`.

## Packaging and CI -- metadata parses, installation claim refuted

1. **Critical, `pyproject.toml` lines 1--3 and CI lines 21--30:** `pip wheel .`
   succeeds but creates `bwtandem-0.9.0-py3-none-any.whl`; no extension build is
   declared in project metadata. CI manually compiles into the checkout, masking
   the fact that users installing the wheel do not receive `_accelerators`.
   Fix with a small tracked `setup.py` (or a backend/config that declares both
   extensions), then make CI run `pip install .`, assert the native modules'
   loaded paths end in an extension suffix, and run tests from outside the repo.
2. **High, CI lines 18--38:** CI never installs the project and never invokes the
   declared `bwtandem` entry point. The current `src.main:main` is importable from
   the checkout, but that does not validate an installed artifact. Build a wheel,
   install it into a clean venv, `cd /tmp`, run `bwtandem --help` and the smoke
   input via an absolute fixture path.
3. **Medium, CI lines 33--34:** the fallback leg runs in the same checkout after
   native compilation and only two test files. Environment gating may work, but
   this is not a packaging parity test. Use a separate matrix job/clean venv and
   run the full suite in both modes.
4. **Low, `pyproject.toml` lines 5--22:** fields are syntactically valid and the
   entry target is importable, but distribution name `bwtandem` exposing package
   name `src` is collision-prone. Rename the import package to `bwtandem` before
   release, then use `bwtandem = "bwtandem.main:main"`.
5. **Confirmed:** `CITATION.cff` has valid CFF 1.2 core fields. Its omitted license
   is consistent with the explicit TODO; add the chosen SPDX identifier to both
   CFF and project metadata when LICENSE lands.

## `results/comparator_baselines.md` -- one factual error

**High, line 43:** “three maize BEDs are byte-identical copies of one run” is
false for TRASH. Manifest rows 35, 40, and 43 have distinct sizes and hashes:
`512556/c497...`, `509733/7063...`, and `507355/14ef...`. The byte-identical
statement belongs to **tantan** (manifest rows 33, 39, 56 all
`68279539/6d3f...`). Correct line 43 to:

```markdown
| TRASH-dn/tpl | n/a | published | experiment/template-specific BEDs; see manifest rows 35, 40, 43, 54, 55 |
```

and add to the tantan row: “3A/3B/3C published files are byte-identical.”
