# 2026 tool additions — longdust + AniAnn's (Section 2.2.4)

Preregistered protocol: `docs/2026-09-01-longdust-anianns-benchmark-protocol.md`
(fixed before any run was scored). These deposits feed the 2026 rows of
Tables 1a, 2, 3B and 3C.

## Versions and runs

| Tool | Version | Run job | Arms |
|---|---|---|---|
| longdust | 1.4 (git g9491215), built on an el7 node: conda gcc, `-O3 -std=c99 -I/usr/include -L/usr/lib64 -lz -lm`, binary sha256 `3469655a…` | 6146419 | default; `-k8 -w20000` × human/maize/Col-CEN |
| AniAnn's | 0.7.1 (git 0d79851), shared conda env `anianns` (bioconda pysam 0.24.0), no `--classify` DB | 6146422 (human), 6146483 (Col-CEN), 6146484 (maize) | defaults, `-j 2` |

Inputs are the chromosome-only FASTAs of the regenerated BWTandem runs, so
cost cells are comparable to the BWTandem rows and NOT to the 2024
competitor rows (which consumed the accession-flavoured FASTA, larger by 3.92% in bases:
3,209,286,105 against 3,088,269,832, or 3.80% counting unambiguous bases only).
Reproducer warning: AniAnn's exits 0 having annotated nothing when its
pysam `.fai` build fails on an unwritable FASTA directory (first attempt,
jobs 6146420/6146421) — symlink the FASTA beside a pre-built index.

## Costs (GNU time, per arm)

| Arm | Wall clock | Peak RSS |
|---|---|---|
| longdust human default / `-k8 -w20000` | 27:56 / 2:00:36 | 0.47 GiB |
| longdust maize default / `-k8 -w20000` | 18:00 / 3:53:40 | 0.58 GiB |
| longdust Col-CEN default / `-k8 -w20000` | 1:56 / 10:23 | 0.06 GiB |
| AniAnn's human / maize / Col-CEN | 2:19:13 / 1:34:26 / 7:18 | 0.53 / 0.55 / 0.48 GiB |

## Conversion and scoring

`scripts/scoring/convert_2026_tools.py` maps raw outputs onto the scorers'
period-column contract — longdust's converted BED is deliberately 3-column
(its paper: "longdust is unable to report the repeat units in case of
tandem repeats"), so every period-conditional metric is structurally empty
rather than fabricated; AniAnn's carries its inferred monomer length (the
array period) as the period column and no motif. Scoring job **6146742** (repo `7b2113d`,
0 dirty; converted-BED hashes in the job log) loaded the FROZEN scorers
(`score_table1.py`, `score_colcen.py` via its LABEL:PATH interface,
`score_maize_postmerge.py`) and overrode only their source lists
(`scripts/scoring/score_2026_tools.py`), so the new rows come from exactly
the code that produced the published rows.

### The human log was replaced, and here is why that is safe to trust

`score_2026_human.txt` is the log of the clean re-run, job **6147179** (repo
`89ada02`, 0 dirty, 2026-09-02, node cpu-14). It replaced job **6146742**'s log,
which ended in a `KeyError`. Swapping a deposited log is exactly where a reader
should not have to take the authors' word, so the crashed log is deposited beside
it as `score_2026_human_6146742_crashed.txt` and the equivalence is checkable:

| file | SHA-256 (first 16) | lines |
|---|---|--:|
| `score_2026_human_6146742_crashed.txt` (job 6146742) | `2ff43f44ff53699f…` | 47 |
| `score_2026_human.txt` (job 6147179) | `87bdf4f4fd420f10…` | 33 |

```
diff score_2026_human_6146742_crashed.txt score_2026_human.txt
2,3c2,3
<  generated  : 2026-09-02T07:57:33-0700      >  generated  : 2026-09-02T15:15:27-0700
<  host       : cpu-50                        >  host       : cpu-14
34,47d33
<  (the 14-line traceback, absent from the re-run)
```

Those are the only differences. Every scored line — the PROVENANCE block's inputs
and byte counts, and all three result blocks — is byte-identical; discounting the
two header lines that record when and where the job ran, the sole remaining
difference in the shared region is `host`. No cell moved.

Three further points a sceptical reader needs, since the claim is about what the
`--adj` flag can and cannot touch:

- **The scorer itself did not change.** `git diff 7b2113d..89ada02 --
  scripts/scoring/score_table1.py` is empty. The only source change in that
  interval is the four-line wrapper edit in `score_2026_tools.py` that adds
  `--adj ""`, plus the scratch directory being added and then removed.
- **`--adj` reaches only the corroborator loop.** In `score_table1.py`, `adj_rules`
  is used at `for rule in adj_rules:` (line 408) and nowhere earlier; the crash was
  at line 421 inside that loop, on `prepared[("ULTRA", variant)]`. The baseline,
  matched-range and stratum blocks — the ones these tables read — are computed and
  printed before that loop begins, which is why they appear in full in the crashed
  log too. Passing an empty `--adj` therefore removes a section that a
  sources-only run cannot compute at all (it needs the ULTRA and tantan baselines
  in `SOURCES`), and removes nothing else.
- **The scratch directory is a separate event.** Job 6147179 regenerated its
  outputs and, on a clean exit, `score_table1.py` removed its own
  `scripts/scoring/work/` scratch, which the crashed run had left behind. Deleting
  that directory *from the repository* was a commit, `bd32c8c`, not something a
  job did; `2b2beea` had committed those 104 MB of intermediate BEDs by mistake.
  The directory holds only scorer intermediates (`gt.bed` and per-tool
  full/p100/strat BEDs), all regenerable by re-running the scorer.

## Headline results (full outputs in the three .txt files here)

- **Human (adotto, full range, one-base rule)**: longdust 41.16% region
  recall / 57.16% precision (half of BWTandem's 80.53% recall at similar
  precision); `-k8 -w20000` 35.42/60.08; AniAnn's 0.10/83.79 — expected,
  its ≥1 kb window cannot see short arrays (its own stated limitation).
- **Col-CEN**: AniAnn's posts the table's highest centromere coverage
  (86.97%) from 60 array-level intervals; longdust matches the leaders'
  monomer recall (99.86%) with no repeat units (CEN180 count 0,
  structural).
- **Maize satellites**: longdust emits **zero calls inside any curated
  knob180/TR-1/CentC array** in either arm — these satellites are not
  low-complexity under its model. AniAnn's leads array-level coverage
  (knob180 87.61%, TR-1 67.98%, CentC 81.42% versus BWTandem's
  79.79/50.34/58.55) with boundary offsets of 14.9–83.7 kb versus
  BWTandem's 253 bp (knob180, gap 0).

The complementary profile these numbers draw — interval-only speed
(longdust), array-level coverage with coarse boundaries (AniAnn's),
per-copy structured records at base-pair boundaries (BWTandem) — is the
positioning stated in the Related Work.
