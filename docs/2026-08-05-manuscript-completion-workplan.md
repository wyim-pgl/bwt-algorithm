# Manuscript Completion Workplan (2026-08-05)

> **For agentic workers:** Use `superpowers:executing-plans` (inline) or
> `superpowers:subagent-driven-development` (fresh subagent per task) to work through
> this task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close every remaining defect in `manuscript.md` so that no sentence
contradicts the tables, the logs, or the measurements behind them.

**Architecture:** Tasks A–E are compute-zero prose/number corrections against evidence
already on disk, and can be done in any order — they are listed in the order that
minimises re-reading. Tasks F–G are blocked on two SLURM jobs still running. Task H is
metadata. Task I is the final consistency gate and must run last.

**Tooling:**
- Manuscript: `/data/gpfs/assoc/pgl/devel/bwt-algorithm/manuscript.md` (untracked; back up before editing)
- Evidence root: `/data/gpfs/assoc/pgl/devel/exp1_human/wp0/`
- Scorers: `wp0/score_table1.py`, `wp0/score_maize_3a.py`, `wp0/score_colcen.py`
- Python: `/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python`
- Competitor benchmark scripts: `/data/gpfs/assoc/pgl/filip/old_bwt/bwtandem/run_bench_*.sbatch`

## Global Constraints

- **Honest science only — no benchmark gaming, no GT overfitting, no hiding
  unfavourable results.** Several corrections below weaken BWTandem's claims. They stay in.
- **Never edit a number without naming the file it came from.** Every task's final step
  records the evidence path in the commit message or in `resume.md`.
- **Back up before the first edit of a session:** `cp -n manuscript.md manuscript.md.<tag>-bak`.
- **Apply prose edits with assertion-checked replacement**, never a bare `sed s///g` on
  a phrase that may occur more than once. Pattern to reuse (Task A shows it in full):
  count occurrences, abort if the count is not what you expect, only then write.
- **`int()`, never `int(float())`, when reading a period from a BED.** BWTandem's own
  BED is 8 columns whose column 5 is a copy count; the period is `len(motif)`.
- Do not re-open the hardware-asymmetry disclosure (node generation / thread count /
  input FASTA). **The user decided on 2026-08-05 that it is not needed.**

## Edit Map

| File | Responsibility | Tasks |
|---|---|---|
| `manuscript.md` §4.4 L541, §4.5 L546 | unique-region period claim | A |
| `manuscript.md` §2.2.1 L82 | thread count disclosure | B |
| `manuscript.md` §3.2 L228 | TRF period cap on Arabidopsis | C |
| `manuscript.md` §3.3.1 L320 | maize TAG definition provenance | D |
| `manuscript.md` §4.3 L534–535 | TRF density-vs-size argument | E |
| `manuscript.md` §3.1 L99, Table 1d, Figure 1 | single-knob operating-point curve | F |
| `manuscript.md` §3.1 L95, §1 L42 | index reuse across operating points | G |
| `results/manifest.tsv` | memory metric provenance | H |
| all of the above | final consistency sweep | I |

---

### Task A: §4.4 and §4.5 — the unique-region period claim is inverted

The manuscript says the regions only BWTandem reports are **short**-period. After the
`period_of()` scorer fix they are **longer**-period than the corroborated calls. The
same wrong claim appears in two sections. Fixing it *strengthens* the surrounding
"not low-complexity artefacts" conclusion, which is why it must not be quietly dropped.

**Files:**
- Modify: `manuscript.md:541` (§4.4 Limitations), `manuscript.md:546` (§4.5 Future directions)
- Evidence: `exp1_human/wp0/logs/wp31_refix_5985651.out` lines 5–20

**Ground truth (BWTandem-F, 875,224 unique calls = 21.9% of output):**

| period band | unique | shared |
|---|--:|--:|
| 1–2 bp | **13.2%** | **24.8%** |
| 3–6 bp | 61.2% | 51.7% |
| 7–20 bp | 3.7% | 17.4% |
| 21–100 bp | **21.9%** | **6.1%** |

The manuscript's "a quarter of them mononucleotide or dinucleotide" is the **shared**
figure (24.8%), misattributed to the unique set (13.2%).

- [x] **Step 1: Re-read the evidence before touching prose**

```bash
sed -n '1,25p' /data/gpfs/assoc/pgl/devel/exp1_human/wp0/logs/wp31_refix_5985651.out
```

Expected: line 10 `period  3-6:61.2%  21-100:21.9%  1-2:13.2%  7-20:3.7%` (unique),
line 15 `period  3-6:51.7%  1-2:24.8%  7-20:17.4%  21-100:6.1%` (shared).
If these differ, stop and re-derive — the numbers below are copied from this log.

- [x] **Step 2: Back up**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
cp -n manuscript.md manuscript.md.pre-wpA-bak
```

- [x] **Step 3: Apply both edits with assertion-checked replacement**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
python3 - <<'PY'
import io,sys
p="manuscript.md"; s=io.open(p,encoding="utf-8").read()
edits=[
 ("They are short-period arrays, a quarter of them mononucleotide or dinucleotide, and 91% of them lie outside the adotto catalog,",
  "They skew toward longer periods than the corroborated calls, not shorter: 21.9% have a period of 21-100 bp against 6.1% of shared calls, while 13.2% are mononucleotide or dinucleotide against 24.8%. 91% of them lie outside the adotto catalog,",1),
 ("and are now characterised well enough to design it: short period, higher motif diversity than corroborated calls, and 91% outside the catalog.",
  "and are now characterised well enough to design it: longer periods than the corroborated calls, higher motif diversity, and 91% outside the catalog.",1),
]
ok=True
for old,new,want in edits:
    got=s.count(old)
    print(("  ok" if got==want else f"  !! found {got}, expected {want}")+f": {old[:60]!r}")
    if got!=want: ok=False
    else: s=s.replace(old,new)
if not ok: print("\nABORTED - no write"); sys.exit(1)
io.open(p,"w",encoding="utf-8").write(s); print("\nwritten:",p)
PY
```

- [x] **Step 4: Verify the sentence reads correctly and the claim no longer inverts**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
grep -nE "skew toward longer periods|longer periods than the corroborated" manuscript.md
grep -nE "short-period arrays|short period, higher motif diversity" manuscript.md || echo "  OK: no stale short-period claim"
```

Expected: two hits from the first grep, `OK:` from the second.

- [x] **Step 5: Check the surrounding argument still follows**

Read `manuscript.md:541` in full. The sentence after the edit must still support "This
benchmark cannot adjudicate them." It does — longer periods and higher entropy both cut
against the low-complexity reading — but confirm no connective ("therefore", "so") now
points the wrong way.

- [x] **Step 6: Record in resume.md** under 다음 작업 as `[x] Task A` with the evidence path.

---

### Task B: Methods §2.2.1 — thread count does not match the runs behind Table 1b/1d

L82 states two threads for everything. Per `results/manifest.tsv`, that is true of
Tables 1a/1c/2/3A–3C (jobs 5983792/93/94) but **false of Table 1b, Table 1d and
Figure 1**, whose four runs used four threads.

**Files:**
- Modify: `manuscript.md:82`
- Evidence: `results/manifest.tsv` column 9 (`threads`)

- [x] **Step 1: Confirm the split**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
awk -F'\t' 'NR==1||$2~/BWTandem/{print $1"\t"$2"\t"$8"\t"$9}' results/manifest.tsv
```

Expected: `1a/1c/2/3A/3B/3C → 2 threads`, `1b + all four 1d rows → 4 threads`.

- [x] **Step 2: Confirm output is thread-count independent, so only cost is affected**

```bash
grep -aE "Total repeats found" /data/gpfs/assoc/pgl/devel/exp1_human/wp0/logs/*5982391*.out \
                               /data/gpfs/assoc/pgl/devel/exp1_human/wp0/logs/*5982392*.out 2>/dev/null
```

Expected: both `4,009,261` (1-thread run 5982392 and 2-thread run 5982391 agree). This
licenses the phrase "affects cost only". If they differ, say so instead and stop.

- [x] **Step 3: Apply the edit**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
python3 - <<'PY'
import io,sys
p="manuscript.md"; s=io.open(p,encoding="utf-8").read()
old="BWTandem was configured for the full period range of 1 to 2,000 bp using two threads."
new=("BWTandem was configured for the full period range of 1 to 2,000 bp using two threads "
     "for the whole-genome comparisons of Tables 1a, 2 and 3, and four threads for the "
     "matched-range operating-point runs behind Table 1b, Table 1d and Figure 1. Detection "
     "output is identical across thread counts, which we verified at one and two threads "
     "(4,009,261 calls in both); thread count changes runtime and peak memory only.")
assert s.count(old)==1, f"found {s.count(old)}"
io.open(p,"w",encoding="utf-8").write(s.replace(old,new)); print("written")
PY
```

- [x] **Step 4: Check the next sentence does not now contradict it**

```bash
sed -n '82p' manuscript.md | fold -s -w 100
```

The later clause "multithreaded tools including BWTandem and ULTRA were constrained to
two threads" is now inconsistent with the four-thread runs. Rewrite it in place to name
the comparison it applies to (the whole-genome cost tables), or delete it if Step 3's
sentence already carries the information. Do not leave both readings.

- [x] **Step 5: Record in resume.md.**

---

### Task C: ~~§3.2 — "TRF, capped at 200 bp" contradicts the command that was run~~ — RETRACTED

> **This task was based on a wrong premise and its edit was reverted (2026-08-05 23:5x).**
> `run_bench_colcen_trf_mreps_tantan.sbatch` passes `50 500`, but that script is not what
> produced the published Col-CEN row. The run log records the command that did:
> `trf .../Col-CEN_v1.2.fasta 2 7 7 80 10 50 200 -ngs -h`, elapsed `131:57:00`. **200 bp is
> correct**, as were `manifest.tsv`'s `2-200` and the Abstract. Lesson: a run log's
> `Command being timed` outranks any script on disk, because scripts get edited after the run.
> Logged TRF MaxPeriod: Col-CEN 200 · human 500 · maize 3A 6 · 3B 500 · 3C 200.
> The useful residue: §4.3's matched comparison is Arabidopsis 200 bp against maize **3C**
> 200 bp (5:28:56), not 3B — the manuscript now says so.


L228 says TRF was capped at 200 bp on Arabidopsis. The benchmark script passes
`trf $FASTA 2 7 7 80 10 50 500 -ngs -h`, where the trailing `500` is TRF's MaxPeriod.
Methods elsewhere states 500. One of the three is wrong.

**Files:**
- Modify: `manuscript.md:228`
- Evidence: `/data/gpfs/assoc/pgl/filip/old_bwt/bwtandem/run_bench_colcen_trf_mreps_tantan.sbatch`

- [x] **Step 1: Read the actual TRF invocation for every genome**

```bash
grep -hnoE "trf [^>]*" /data/gpfs/assoc/pgl/filip/old_bwt/bwtandem/run_bench_*trf*.sbatch
```

TRF argument order is `File Match Mismatch Delta PM PI Minscore MaxPeriod`, so in
`2 7 7 80 10 50 500` the MaxPeriod is **500**.

- [x] **Step 2: Decide which number is right, then make all three agree**

If MaxPeriod is 500 on Arabidopsis, edit L228 to say 500 and check that the
"more than 250 times longer" ratio and the Table 2 TRF row are unaffected (they are —
they depend on runtime, not on the cap). If a *different* Col-CEN TRF run with cap 200
produced the 131.9 h in Table 2, then instead cite that run and record its path.

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
python3 - <<'PY'
import io,sys
p="manuscript.md"; s=io.open(p,encoding="utf-8").read()
old="TRF, capped at 200 bp, required 131.9 hours"
new="TRF, capped at 500 bp, required 131.9 hours"
assert s.count(old)==1, f"found {s.count(old)}"
io.open(p,"w",encoding="utf-8").write(s.replace(old,new)); print("written")
PY
```

- [x] **Step 3: Sweep for any other TRF period-cap statement**

```bash
grep -nE "capped at|maximum period|max period|500 bp maximum|200 bp" manuscript.md | grep -i trf
```

Expected: every TRF cap statement now reads 500 bp. Reconcile any that do not.

- [x] **Step 4: Record in resume.md.**

---

### Task D: §3.3.1 — disclose that the maize TAG rule was chosen after the first one failed

L320 presents the substring rule and the rotation rule side by side but does not say the
rotation rule was adopted **after** the substring rule was found to inflate short-motif
tools. A reviewer who spots the ordering will read it as rule-shopping; stating it
plainly removes that reading, and the justification is independent of the outcome.

**Files:**
- Modify: `manuscript.md:320`
- Evidence: `exp1_human/wp0/WP0-AC_MAIZE_FINDINGS.md`

- [x] **Step 1: Re-read the finding to get the ordering right**

```bash
grep -naE "substring|rotation|fair" /data/gpfs/assoc/pgl/devel/exp1_human/wp0/WP0-AC_MAIZE_FINDINGS.md | head -20
```

- [x] **Step 2: Append a disclosure sentence to the L320 paragraph**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
python3 - <<'PY'
import io,sys
p="manuscript.md"; s=io.open(p,encoding="utf-8").read()
old="The two definitions give opposite orderings."
new=("The two definitions give opposite orderings. We adopted the rotation rule after "
     "the substring rule was found to count any call whose consensus happens to contain "
     "TAG or CTA, which is why it is reported here alongside the original rather than in "
     "place of it: Table 3A retains the substring rule and Table 3A-b restates the "
     "comparison under the rotation rule, so a reader can see both.")
assert s.count(old)==1, f"found {s.count(old)}"
io.open(p,"w",encoding="utf-8").write(s.replace(old,new)); print("written")
PY
```

- [x] **Step 3: Verify both tables are still present and consistent**

```bash
grep -nE "Table 3A\.|Table 3A-b\." manuscript.md
grep -nE "^\| BWTandem \| 41,314,898" manuscript.md
```

Expected: both captions present; the 3A-b BWTandem row reads
`41,314,898 | **1,508,821** | **34,247** | 235.43`.

- [x] **Step 4: Record in resume.md.**

---

### Task E: §4.3 — verify the TRF density-vs-size argument is not a parameter artefact

§4.3 rests its central design argument on "TRF took 131.9 hours on the 134 Mb
Arabidopsis genome and only 5.2 hours on the 2.18 Gb maize genome". If TRF ran with a
different MaxPeriod on maize, the contrast is confounded and the argument must be
weakened or dropped.

**Files:**
- Modify (conditionally): `manuscript.md:535`
- Evidence: `run_bench_*.sbatch` for maize; `results/manifest.tsv` TRF rows

- [x] **Step 1: Find the maize TRF invocation**

```bash
grep -rhnoE "trf [^>]*" /data/gpfs/assoc/pgl/filip/old_bwt/bwtandem/*.sbatch \
     /data/gpfs/assoc/pgl/filip/bwtandem_results/ 2>/dev/null | sort -u
```

- [x] **Step 2: Branch on the result**

- **Same MaxPeriod on both genomes** → the argument stands. Add one clause naming the
  shared parameters so a reader can check it, e.g. "…at identical TRF parameters
  (`2 7 7 80 10 50 500`)".
- **Different MaxPeriod** → the comparison is confounded. Replace the empirical claim
  with the citation-only version: keep the Wlodzimierz et al. (2023) observation that
  TRF is slow on satellite-rich sequence, drop the two-genome contrast, and note in
  `resume.md` that isolating the effect needs a matched-parameter TRF re-run.

- [x] **Step 3: Confirm §4.3's remaining hedge is still accurate**

The paragraph already says "We have not isolated this effect to kmer seeding
specifically". Ensure that hedge survives whichever branch was taken.

- [x] **Step 4: Record in resume.md, including which branch was taken and why.**

---

### Task F: §3.1 — replace the four-point curve with the real single-knob sweep

**✅ DONE 2026-08-07.** Original job 5985873 was voided (src/ edited mid-run) and
resubmitted as **5994683** (COMPLETED, 5/5 points, per-point src mtimes verified
unchanged). Scored by **6016600** (`score_table1_remeas.py --adj published` +
5 `--extra`). Decision (Step 4): **keep P/B/F/H** — Table 1d is the matched-range
(≤100 bp) axis the figure shares with ULTRA/tantan, and the sweep spans a narrower
range — so §3.1 now says plainly that the four points differ in more than one
setting, and reports the true single-knob sweep alongside: recall 62.59→81.46,
precision 54.77→48.94, plateau at 0.80–0.76, adj-precision peak at 0.72 (79.49).
Figure 1 and its CSV are unchanged because Table 1d is unchanged. Manifest: 5
`sweep` rows added (43 rows, 0 missing). Evidence:
`logs/idsweep_5994683.out`, `logs/idsweep_score_6016600.out`,
`WP-F_SINGLE_KNOB_SWEEP_FINDINGS.md` (full-genome section).

**~~BLOCKED on job 5985873~~ → superseded by 5994683** (`idsweep`, cpu-49, 5 points: off / 0.80 / 0.76 / 0.72 / 0.68,
sharing one index via `BWT_INDEX_CACHE`).

Today's P→B→F→H "curve" changes more than one thing between points — P disables the
catch-all outright and F adds a copies≥3 gate — so calling it a sweep of one threshold
is not supported. The sweep job is the one that is.

**Files:**
- Modify: `manuscript.md:99` (§3.1 curve paragraph), Table 1d (`manuscript.md:202`+),
  Figure 1 caption (`manuscript.md:200`), `wp0/figure_curve.{pdf,png}`, `wp0/figure_curve_data.csv`
- Evidence: `exp1_human/wp0/logs/idsweep_*.out`

- [ ] **Step 1: Confirm the job finished cleanly**

```bash
sacct -j 5985873 --format=JobID%14,JobName%10,State,Elapsed,MaxRSS,NodeList
```

Expected: `COMPLETED`, exit 0. If it failed, stop and diagnose before scoring anything.

- [ ] **Step 2: Check each point's BED is self-consistent with its log**

```bash
cd /data/gpfs/assoc/pgl/devel/exp1_human/wp0
for b in out/idsweep_*.bed; do
  echo "$b  BED=$(wc -l < $b)  LOG=$(grep -ao 'Total repeats found: [0-9]*' ${b%.bed}.log | tail -1)"
done
```

Expected: every BED line count equals its log's `Total repeats found`. A mismatch means
the same defect this whole audit was about — do not score it, find the other run.

- [ ] **Step 3: Score all points on the Table 1 axis**

```bash
cd /data/gpfs/assoc/pgl/devel/exp1_human/wp0
sbatch run_score_table1.sbatch   # ~13 min, 4.4 GB — exceeds the login node's 4 GB cap
```

- [ ] **Step 4: Decide whether the sweep supersedes P/B/F/H**

If the sweep spans a comparable recall/precision range, replace Table 1d and Figure 1
with it and rewrite L99 as a genuine one-parameter sweep. If it spans a narrower range,
**keep P/B/F/H and rewrite L99 to say plainly that the four points differ in more than
one setting** — which is the honest description of what was run. Do not present four
multi-parameter configurations as a single-knob curve either way.

- [ ] **Step 5: Regenerate the figure and its CSV**

```bash
cd /data/gpfs/assoc/pgl/devel/exp1_human/wp0
/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python make_figure_curve.py
```

Verify `figure_curve_data.csv` row count matches the number of points in Table 1d.

- [ ] **Step 6: Update `results/manifest.tsv` with a row per sweep point** (source BED,
      lines, job, threads, period range, elapsed, node, scorer, rule, commit).

- [ ] **Step 7: Record in resume.md.**

---

### Task G: §3.1 / §1 — support or retract "one index yields an operating-point curve"

**BLOCKED on job 5985762** (`cache_bench`, cpu-51, Col-CEN, 4 runs sharing one
`BWT_INDEX_CACHE`).

L42 claims "one index yields an operating-point curve". Across the runs actually
published, each operating point rebuilt the suffix array — the claim was false when
written. `BWT_INDEX_CACHE` was implemented to make it true; this task supplies the
measurement or retracts the claim.

**Files:**
- Modify: `manuscript.md:42` (Introduction), `manuscript.md:95` (§3.1 cost paragraph)
- Evidence: `exp1_human/wp0/logs/cache_bench_*.out`

- [x] **Step 1: Confirm the job finished and read the four wall-clocks**

```bash
sacct -j 5985762 --format=JobID%14,State,Elapsed,MaxRSS,NodeList
grep -haE "Elapsed \(wall clock\)|index cache (hit|miss)|Total repeats found" \
     /data/gpfs/assoc/pgl/devel/exp1_human/wp0/logs/cache_bench_*.out
```

Expected: run 1 a cache miss (builds the index), runs 2–4 hits.

- [x] **Step 2: Verify the cache did not change the output**

```bash
cd /data/gpfs/assoc/pgl/devel/exp1_human/wp0
md5sum out/cachebench_*.bed
```

Expected: identical checksums for runs at identical settings. If they differ, the cache
is not behaviour-preserving and the feature — not the manuscript — is the problem.

- [x] **Step 3: Write the saving into §3.1 as a measured number**

State what index reuse saves as `<miss time> → <hit time>` on Col-CEN, and say the
figure is for a 134 Mb assembly so no one extrapolates it to human unmeasured.

- [x] **Step 4: Make L42's claim match what was measured**

If the cache is on by default, "one index yields an operating-point curve" is now true —
add the pointer to §3.1. If it is opt-in (it is: unset `BWT_INDEX_CACHE` = old
behaviour), say so: "one index *can* serve every operating point, with the cache
enabled". Do not leave the unqualified claim.

- [x] **Step 5: Record in resume.md.**

---

### Task H: `results/manifest.tsv` — label which metric `peak_rss_gb` is

The column holds `sacct MaxRSS` (concurrent cgroup total). `/usr/bin/time -v` on the
same runs reports roughly half that, because `--threads` spawns processes
(`ProcessPoolExecutor`) and `getrusage(RUSAGE_CHILDREN)` returns the max over children
rather than their sum. This session already mistook one for the other once.

**Files:**
- Modify: `results/manifest.tsv` (header), `exp1_human/wp0/make_manifest.py`

- [x] **Step 1: Confirm the column's provenance in the generator**

```bash
grep -nE "peak_rss|MaxRSS|sacct|time -v" /data/gpfs/assoc/pgl/devel/exp1_human/wp0/make_manifest.py
```

- [x] **Step 2: Rename the column to `peak_rss_gb_sacct`** in `make_manifest.py`, and add
      a comment recording that `/usr/bin/time -v` undercounts multi-worker runs.

- [x] **Step 3: Regenerate and confirm the rows are otherwise unchanged**

```bash
cd /data/gpfs/assoc/pgl/devel/exp1_human/wp0
cp ../../bwt-algorithm/results/manifest.tsv /tmp/manifest.before.tsv
/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python make_manifest.py
diff <(cut -f1-11 /tmp/manifest.before.tsv) <(cut -f1-11 ../../bwt-algorithm/results/manifest.tsv) \
  && echo "OK: only the memory column header changed"
```

- [x] **Step 4: Record in resume.md.**

---

### Task I: Final consistency gate — run last

**✅ DONE 2026-08-07.** Step 1: no stale cost/maize figures (grep clean). Step 2:
every table's BWTandem row traces to one execution in `results/manifest.tsv` —
1a/1c → 5983793 (12:05:48 / 21.86 GB), 2 → 5983792 (0:30:49 / 1.95 GB),
3A/3B/3C → 5983794 (15:25:17 / 22.41 GB), 1b/1d-F → 5981960, 1d-P/B/H →
5982255/56/57 — all matching the printed table values. Step 3: Tables 3A/3A-b
agree on 41,314,898 bp / 1,508,821 regions (TAG columns differ by rule, as
intended). Step 4: the 2026-08-06 four-agent audit (`docs/2026-08-06-audit-triage.md`)
was the prose-figure walk; every figure added since (§3.1 single-knob sweep)
traces to `logs/idsweep_score_6016600.out`. Step 5: **115 passed** (the "91"
below is the stale pre-audit baseline). Step 6: resume.md snapshot 2026-08-07.

- [ ] **Step 1: No stale cost figure survives anywhere**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
grep -nE "7\.4 hour|18 minutes|12\.99|14\.37|1\.31 GB|6\.6 hour|0\.31 hour|41,007,584|35,017|quarter of the cost" manuscript.md \
  || echo "  OK: no stale cost or maize figures"
```

- [ ] **Step 2: Every table's BWTandem row is one execution**

For each of Tables 1a, 1b, 1c, 1d, 2, 3A, 3A-b, 3B, 3C, confirm against
`results/manifest.tsv` that accuracy, runtime and memory carry the same `producer_job`.
This is the defect the whole audit was about; it is the one check that must not be skipped.

- [ ] **Step 3: Cross-table agreement on shared quantities**

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
grep -nE "41,314,898|1,508,821" manuscript.md
```

Expected: Tables 3A and 3A-b agree on BWTandem's banded SSR bp and region count. They
differ only in the TAG column (39,011 substring vs 34,247 rotation), which is the point
of having both.

- [ ] **Step 4: Every number quoted in prose exists in a table or a named log**

Walk the Abstract, §1, §3.1, §3.2, §3.3, §4.1–4.5 and confirm each figure traces to a
table cell or an evidence file. Anything that traces to neither gets removed or measured.

- [ ] **Step 5: Full test suite still green** (the tasks above touch no source, but
      Task H edits a script under `wp0/` and Task G may touch the cache path)

```bash
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python -m pytest tests/ -q
```

Expected: **91 passed**.

- [ ] **Step 6: Update `resume.md`** — mark the workplan closed, list what remains open
      (competitor costs were never re-verified; that is a known, accepted limitation),
      and refresh the lab wiki entry per `~/.claude/CLAUDE.md`.

---

## Execution order

```
A → B → C → D → E        (compute 0, any order, do these now)
      ↓
F (needs 5985873)  G (needs 5985762)   ← independent of each other
      ↓
      H
      ↓
      I  (gate — last)
```

## Self-review

- **Coverage:** every open item from `resume.md`'s 다음 작업 maps to a task —
  §4.4 period sentence → A; Methods two-threads → B; TAG fair-definition → D;
  §4.3/4.5 unverified assertions → E and A; curve paragraph → F; cache figure → G;
  manifest metric → H. Task C is new, found while writing this plan.
- **No placeholders:** every edit step carries the exact old and new strings, every
  verification step an exact command and its expected output.
- **Consistency:** the 3A-b row `41,314,898 / 1,508,821 / 34,247 / 235.43` is quoted
  identically in Tasks D and I; the period table in Task A matches
  `wp31_refix_5985651.out` lines 10 and 15 verbatim.
- **Known gap, deliberately not a task:** competitor runtimes and memory were never
  re-verified, and ran on a different node class, thread count and input FASTA. The user
  decided on 2026-08-05 that this asymmetry does not need disclosing. It is recorded
  here so the decision is visible, not so it gets re-litigated.
