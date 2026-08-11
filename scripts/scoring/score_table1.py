#!/usr/bin/env python3
"""WP0-B — Table 1 rescored over a period range every tool could actually search.

The published Table 1 compares tools that were run over different period ranges
(BWTandem 1-2000, TRF 2-500, ULTRA `-p 100`, tantan no period model). Measured on
the BEDs themselves:

    tool        calls p>100   % of calls    bp p>100   % of bp
    bwtandem         98,624        2.47%   120,089,835   35.36%
    trf              54,941        5.43%   200,273,903   67.63%
    tantan                0        0.00%             0    0.00%
    ultra                 0        0.00%             0    0.00%

So a third of BWTandem's bp and two thirds of TRF's sit in a range the other two
tools emit nothing in. The bp-level comparison is not like-for-like, and the
adjusted-precision figure inherits the same problem: its corroborators are ULTRA
and tantan, neither of which can corroborate anything above period 100.

This script produces, from one code path and with provenance:

  1. BASELINE   — every tool as published (full range)
  2. MATCHED    — every tool restricted to period <= 100
  3. ADJUSTED   — precision counting adotto OR corroborator support, at both
                  ranges, under two corroborator rules (see --adj)
  4. STRATUM    — period 101-2000 only: BWTandem vs TRF, the only two tools there

Metric definitions are those of ../score_overlap.py (itself a reproduction of
filip's scripts/compute_metrics.py):

  region_recall    = GT regions overlapped by >=1 tool call / total GT regions
  region_precision = tool calls overlapping >=1 GT region / total tool calls
  bp_recall        = intersect_bp(GT, merged_tool) / GT_bp
  bp_precision     = intersect_bp(GT, merged_tool) / merged_tool_bp
  adj_precision    = (calls overlapping GT + calls that do not but are supported
                      by a corroborator) / total calls        [fp_check.py]

They are implemented here as streaming counts rather than by slurping bedtools
output into memory, so this runs in O(1) RAM; `--check-defs` pins the two
implementations against each other on a single chromosome.

Usage:
  score_table1.py                      # full genome, both adj rules
  score_table1.py --chroms chr21       # single chromosome (fast smoke test)
  score_table1.py --check-defs chr21   # cross-check vs score_overlap.py, then exit
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
EXP1 = os.path.dirname(HERE)
BEDTOOLS = os.environ.get(
    "BEDTOOLS", "/data/gpfs/assoc/pgl/bin/bedtools2/bin/bedtools")

GT = os.path.join(EXP1, "data/adotto_primary.bed")
FILIP = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"

# The 24 chromosomes the adotto catalog covers. Tool BEDs carry extra scaffolds
# and alt contigs; scoring them against a GT that has no rows there would count
# every such call as a false positive, so they are dropped from both sides.
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

SOURCES = [
    ("BWTandem", f"{FILIP}/bwtandem/bwt_hg38.bed"),
    ("TRF",      f"{FILIP}/trf/GCA_000001405.15_GRCh38_genomic_output.bed"),
    ("ULTRA",    f"{HERE}/beds/ultra_human_GCA.bed"),
    ("tantan",   f"{FILIP}/tantan/GCA_000001405.15_GRCh38_genomic_output.bed"),
    # Table 1 prints TRASH twice (de novo and template); only one human TRASH BED
    # exists on disk and the two rows are byte-identical, so it is scored once.
    # It belongs here and not in a footnote: 77.19 % of its bp are period >100.
    ("TRASH",    f"{FILIP}/trash/GCA_000001405.15_GRCh38_genomic_trash.bed"),
]

# Rows that are scored but never act as corroborators. A second run of a tool that
# is already in SOURCES belongs here: letting it corroborate would mean a tool
# vouching for itself, which is the failure the published adjusted-precision rule
# avoids by excluding self. Populated by --extra.
EXTRA = []

# mreps is Table 1's seventh row and is deliberately NOT scored here. Its human
# BED is 13 GB / 18.2 M calls at 1.75 % region recall — the shape of the same
# defect WP0 found on Col-CEN, where mreps had been run at `-maxperiod 6` and a
# re-run at 150-400 moved CEN180 recall from 5.76 % to 78.45 %. Rescoring a run
# whose parameters are the problem would launder the defect; the fix is a re-run.
# See MREPS_HUMAN below in the log for the measured period distribution.

MAX_MATCHED_PERIOD = 100   # the lowest cap any tool was run under (ULTRA -p 100)
STRATUM = (101, 2000)      # above the cap, up to BWTandem's own ceiling


# --------------------------------------------------------------------------
# bedtools plumbing — every helper streams, nothing is held in memory
# --------------------------------------------------------------------------

def _stream(cmd):
    """Yield stdout lines of cmd; raise on non-zero exit."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line
    p.stdout.close()
    err = p.stderr.read()
    p.stderr.close()
    if p.wait() != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed:\n{err}")


def count_lines(cmd):
    return sum(1 for _ in _stream(cmd))


def count_bp(cmd):
    total = 0
    for line in _stream(cmd):
        f = line.split("\t")
        if len(f) >= 3:
            total += int(f[2]) - int(f[1])
    return total


def file_regions_bp(path):
    regions = bp = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                try:
                    bp += int(p[2]) - int(p[1])
                except ValueError:
                    continue
                regions += 1
    return regions, bp


def period_of(parts):
    """Period from column 5, falling back to motif length.

    Required, not cosmetic: TRF and TRASH put a full consensus (up to tens of kb)
    in the motif column and the period in column 5, so a len()-derived period is
    wrong for them by orders of magnitude.
    """
    if len(parts) >= 5:
        try:
            return int(parts[4])
        except ValueError:
            pass
    return len(parts[3]) if len(parts) >= 4 else None


def prepare(src, dst, chroms, pmin=None, pmax=None):
    """Write a chrom-filtered, period-filtered, coordinate-sorted copy of src."""
    keep = set(chroms)
    tmp = dst + ".unsorted"
    filtering = pmin is not None or pmax is not None
    n = 0
    with open(src) as fin, open(tmp, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            # 3 columns, not 4: the ground truth carries no motif or period
            # column, and must survive the no-filter path unchanged.
            if len(parts) < 3 or parts[0] not in keep:
                continue
            if filtering:
                per = period_of(parts)
                if per is None:
                    continue
                if pmin is not None and per < pmin:
                    continue
                if pmax is not None and per > pmax:
                    continue
            fout.write(line)
            n += 1
    with open(dst, "w") as out:
        r = subprocess.run([BEDTOOLS, "sort", "-i", tmp], stdout=out,
                           stderr=subprocess.PIPE, text=True)
    os.remove(tmp)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools sort failed on {src}:\n{r.stderr}")
    return n


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------

def score(gt, gt_regions, gt_bp, tool, workdir, tag):
    """region/bp recall + precision for one prepared tool BED."""
    n_calls = file_regions_bp(tool)[0]
    if n_calls == 0:
        return dict(calls=0, reg_recall=0.0, reg_prec=0.0,
                    bp_recall=0.0, bp_prec=0.0, tool_bp=0)

    hit_gt = count_lines([BEDTOOLS, "intersect", "-a", gt, "-b", tool, "-wa", "-u"])
    hit_tool = count_lines([BEDTOOLS, "intersect", "-a", tool, "-b", gt, "-wa", "-u"])

    merged = os.path.join(workdir, f"{tag}.merged.bed")
    with open(merged, "w") as out:
        r = subprocess.run([BEDTOOLS, "merge", "-i", tool], stdout=out,
                           stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools merge failed on {tool}:\n{r.stderr}")
    inter_bp = count_bp([BEDTOOLS, "intersect", "-a", gt, "-b", merged])
    tool_bp = file_regions_bp(merged)[1]
    os.remove(merged)

    return dict(
        calls=n_calls,
        reg_recall=100.0 * hit_gt / gt_regions if gt_regions else 0.0,
        reg_prec=100.0 * hit_tool / n_calls,
        bp_recall=100.0 * inter_bp / gt_bp if gt_bp else 0.0,
        bp_prec=100.0 * inter_bp / tool_bp if tool_bp else 0.0,
        tool_bp=tool_bp,
    )


def adjusted_precision(gt, tool, corroborators, workdir, tag):
    """(calls in GT + calls outside GT that a corroborator also calls) / calls.

    This is fp_check.py's definition. It exists because the adotto catalog is
    incomplete: a call adotto lacks is not automatically wrong. Its weakness —
    the reason WP0-B revisits it — is that a corroborator can only corroborate
    what it was able to search.
    """
    n_calls = file_regions_bp(tool)[0]
    if n_calls == 0 or not corroborators:
        return None
    fp = os.path.join(workdir, f"{tag}.fp.bed")
    with open(fp, "w") as out:
        r = subprocess.run([BEDTOOLS, "intersect", "-a", tool, "-b", gt, "-v"],
                           stdout=out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools intersect -v failed:\n{r.stderr}")
    n_fp = file_regions_bp(fp)[0]
    n_supp = 0
    if n_fp:
        n_supp = count_lines([BEDTOOLS, "intersect", "-a", fp, "-b",
                              *corroborators, "-u"])
    os.remove(fp)
    n_tp = n_calls - n_fp
    return dict(
        calls=n_calls, tp=n_tp, fp=n_fp, supported=n_supp,
        adj_prec=100.0 * (n_tp + n_supp) / n_calls,
        supp_frac=100.0 * n_supp / n_fp if n_fp else 0.0,
    )


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

HDR = (f"{'tool':<10} {'calls':>10} {'regRecall%':>11} {'regPrec%':>9} "
       f"{'bpRecall%':>10} {'bpPrec%':>8}")


def print_block(title, rows):
    print(f"\n########## {title} ##########")
    print(HDR)
    for name, m in rows:
        print(f"{name:<10} {m['calls']:>10} {m['reg_recall']:>11.2f} "
              f"{m['reg_prec']:>9.2f} {m['bp_recall']:>10.2f} {m['bp_prec']:>8.2f}")
    sys.stdout.flush()


def print_adj(title, rows):
    print(f"\n########## {title} ##########")
    print(f"{'tool':<10} {'calls':>10} {'inGT':>10} {'outGT':>10} "
          f"{'supported':>10} {'%ofOut':>8} {'adjPrec%':>9}")
    for name, a in rows:
        if a is None:
            continue
        print(f"{name:<10} {a['calls']:>10} {a['tp']:>10} {a['fp']:>10} "
              f"{a['supported']:>10} {a['supp_frac']:>7.1f}% {a['adj_prec']:>9.2f}")
    sys.stdout.flush()


# --------------------------------------------------------------------------
# definition cross-check
# --------------------------------------------------------------------------

def check_defs(chrom):
    """Run score_overlap.py and this module on one chromosome; compare."""
    print(f"# cross-checking streaming metrics vs score_overlap.py on {chrom}")
    ref = subprocess.run(
        [sys.executable, os.path.join(EXP1, "score_overlap.py"), GT]
        + [f"{src}:{name}" for name, src in SOURCES]
        + ["--chroms", chrom],
        capture_output=True, text=True)
    if ref.returncode != 0:
        print(ref.stderr)
        raise SystemExit("score_overlap.py failed")
    print("\n--- score_overlap.py (reference) ---")
    print(ref.stdout.rstrip())

    ours = run(chroms=[chrom], adj_rules=[], quiet_header=True)
    print("\n--- score_table1.py (streaming) ---")
    print(HDR)
    for name, m in ours["baseline"]:
        print(f"{name:<10} {m['calls']:>10} {m['reg_recall']:>11.2f} "
              f"{m['reg_prec']:>9.2f} {m['bp_recall']:>10.2f} {m['bp_prec']:>8.2f}")

    # parse reference table and diff
    refvals = {}
    for line in ref.stdout.splitlines():
        f = line.split()
        if len(f) == 6 and f[0] in dict(SOURCES):
            refvals[f[0]] = [float(x) for x in f[2:6]]
    bad = []
    for name, m in ours["baseline"]:
        got = [m["reg_recall"], m["reg_prec"], m["bp_recall"], m["bp_prec"]]
        exp = refvals.get(name)
        if exp is None:
            bad.append(f"{name}: missing from reference output")
            continue
        for label, g, e in zip(("regRecall", "regPrec", "bpRecall", "bpPrec"),
                               got, exp):
            if abs(g - e) > 0.005:
                bad.append(f"{name}.{label}: {g:.2f} != {e:.2f}")
    print()
    if bad:
        for b in bad:
            print(f"MISMATCH {b}")
        raise SystemExit("definitions differ — do not trust the full run")
    print("OK — both implementations agree on all four metrics for every tool")


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def run(chroms, adj_rules, workdir=None, quiet_header=False):
    workdir = workdir or os.path.join(HERE, "work")
    os.makedirs(workdir, exist_ok=True)

    gt = os.path.join(workdir, "gt.bed")
    prepare(GT, gt, chroms)
    gt_regions, gt_bp = file_regions_bp(gt)
    if gt_regions == 0:
        raise SystemExit(
            f"ground truth is empty after filtering to {','.join(chroms)} — "
            f"every metric would silently be 0")

    if not quiet_header:
        print("########## PROVENANCE ##########")
        print(f"generated  : {time.strftime('%Y-%m-%dT%H:%M:%S%z')}")
        print(f"host       : {os.uname().nodename}")
        print(f"bedtools   : {BEDTOOLS}")
        print(f"GT         : {GT}")
        print(f"GT regions : {gt_regions}   GT bp: {gt_bp}")
        print(f"chroms     : {','.join(chroms)}")
        print(f"matched cap: period <= {MAX_MATCHED_PERIOD}")
        print(f"stratum    : period {STRATUM[0]}-{STRATUM[1]}")
        print("period src : column 5, else len(motif)")
        print("sources    :")
        for name, src in SOURCES + EXTRA:
            st = os.stat(src)
            print(f"  {name:<10} {src}")
            print(f"  {'':<10} {st.st_size:,d} bytes  mtime "
                  f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st.st_mtime))}")
        sys.stdout.flush()

    # --- prepare every variant up front -----------------------------------
    prepared = {}   # (tool, variant) -> path
    for name, src in SOURCES + EXTRA:
        for variant, (pmin, pmax) in (
                ("full", (None, None)),
                ("p100", (None, MAX_MATCHED_PERIOD)),
                ("strat", STRATUM)):
            dst = os.path.join(workdir, f"{name}.{variant}.bed")
            prepare(src, dst, chroms, pmin, pmax)
            prepared[(name, variant)] = dst

    results = {}
    for variant, title in (("full", "BASELINE (full period range, as published)"),
                           ("p100", f"MATCHED RANGE (period <= {MAX_MATCHED_PERIOD}, all tools)")):
        rows = []
        for name, _ in SOURCES + EXTRA:
            rows.append((name, score(gt, gt_regions, gt_bp,
                                     prepared[(name, variant)], workdir,
                                     f"{name}.{variant}")))
        results["baseline" if variant == "full" else "matched"] = rows
        if not quiet_header:
            print_block(title, rows)

    # --- period 101-2000 stratum ------------------------------------------
    strat_rows = []
    for name, _ in SOURCES + EXTRA:
        p = prepared[(name, "strat")]
        if file_regions_bp(p)[0] == 0:
            continue
        strat_rows.append((name, score(gt, gt_regions, gt_bp, p, workdir,
                                       f"{name}.strat")))
    results["stratum"] = strat_rows
    if not quiet_header and strat_rows:
        print_block(f"STRATUM period {STRATUM[0]}-{STRATUM[1]} "
                    f"(only tools with calls there)", strat_rows)

    # --- adjusted precision -----------------------------------------------
    # EXTRA rows can take the published rule (its corroborators are the fixed pair
    # ULTRA+tantan, and an extra row is neither of them) but not leave-one-out,
    # where the corroborator set would include the tool's own other run.
    if EXTRA and adj_rules and not quiet_header:
        print(f"\n(leave-one-out adjusted precision omits "
              f"{', '.join(n for n, _ in EXTRA)}: under that rule a tool's own "
              f"other run would corroborate it. The published rule is fine for "
              f"them — its corroborators are ULTRA and tantan.)")
    for rule in adj_rules:
        for variant, label in (("full", "full range"),
                               ("p100", f"period <= {MAX_MATCHED_PERIOD}")):
            rows = []
            scored = SOURCES + EXTRA if rule == "published" else SOURCES
            for name, _ in scored:
                if rule == "published":
                    # The manuscript's rule (Methods L73): predictions outside
                    # the catalog are corroborated by ULTRA and tantan. Applied
                    # to ULTRA or tantan itself a tool would corroborate its own
                    # calls and score 100 %, and the manuscript's ULTRA (78.44)
                    # and tantan (82.93) rows are not 100 %, so self must be
                    # excluded. Dropping self reproduces all four published rows.
                    corr = [prepared[(c, variant)] for c in ("ULTRA", "tantan")
                            if c != name]
                else:  # leave-one-out: every other tool corroborates
                    corr = [prepared[(c, variant)] for c, _ in SOURCES
                            if c != name]
                rows.append((name, adjusted_precision(
                    gt, prepared[(name, variant)], corr, workdir,
                    f"{name}.{variant}.{rule}")))
            results[f"adj_{rule}_{variant}"] = rows
            if not quiet_header:
                rule_txt = ("corroborators ULTRA+tantan" if rule == "published"
                            else "corroborators = all other tools")
                print_adj(f"ADJUSTED PRECISION — {label}, {rule_txt}", rows)

    for path in prepared.values():
        try:
            os.remove(path)
        except OSError:
            pass
    try:
        os.remove(gt)
        os.rmdir(workdir)
    except OSError:
        pass
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=None,
                    help="comma-separated subset (default: all 24)")
    ap.add_argument("--check-defs", metavar="CHROM", default=None,
                    help="cross-check against score_overlap.py on CHROM, then exit")
    ap.add_argument("--adj", default="published,loo",
                    help="adjusted-precision rules: published, loo, or both")
    ap.add_argument("--extra", action="append", default=[], metavar="NAME:PATH",
                    help="additional row, scored but never used as a corroborator "
                         "(repeatable)")
    ap.add_argument("--workdir", default=None)
    a = ap.parse_args()

    for spec in a.extra:
        name, _, path = spec.partition(":")
        if not path:
            raise SystemExit(f"--extra needs NAME:PATH, got {spec!r}")
        if not os.path.exists(path):
            raise SystemExit(f"--extra path does not exist: {path}")
        EXTRA.append((name, path))

    if a.check_defs:
        check_defs(a.check_defs)
        return

    chroms = a.chroms.split(",") if a.chroms else CHROMS
    rules = [r for r in a.adj.split(",") if r]
    run(chroms=chroms, adj_rules=rules, workdir=a.workdir)


if __name__ == "__main__":
    main()
