#!/usr/bin/env python3
"""Region/bp recall + precision for a tool BED vs a ground-truth BED.

Reproduces the exact definitions used in filip's scripts/compute_metrics.py so
results are directly comparable to the published benchmark:

  region_recall    = GT regions overlapped by >=1 tool call / total GT regions
  region_precision = tool calls overlapping >=1 GT region / total tool calls
  bp_recall        = intersect_bp(GT, merged_tool) / GT_bp
  bp_precision     = intersect_bp(GT, merged_tool) / merged_tool_bp

Uses bedtools (path via $BEDTOOLS or PATH). GT/tool BEDs may be optionally
restricted to a set of chromosomes with --chroms for single-chromosome work.

Usage:
  score_overlap.py GT.bed TOOL.bed[:name] [TOOL2.bed:name ...] [--chroms chr21,chr22]
"""
import os
import sys
import subprocess
import tempfile

BEDTOOLS = os.environ.get("BEDTOOLS", "/data/gpfs/assoc/pgl/bin/bedtools2/bin/bedtools")


def sh(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"ERR {' '.join(cmd)}\n{r.stderr}\n")
    return r.stdout


def regions_bp(path):
    regions = bp = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                try:
                    regions += 1
                    bp += int(p[2]) - int(p[1])
                except ValueError:
                    pass
    return regions, bp


def filter_chroms(path, chroms):
    """Write a chrom-filtered, sorted copy; return its path (caller deletes)."""
    keep = set(chroms)
    tf = tempfile.NamedTemporaryFile("w", delete=False, suffix=".bed")
    with open(path) as f:
        for line in f:
            p = line.split("\t")
            if len(p) >= 3 and p[0] in keep:
                tf.write(line)
    tf.close()
    srt = sh([BEDTOOLS, "sort", "-i", tf.name])
    with open(tf.name, "w") as f:
        f.write(srt)
    return tf.name


def region_recall(gt, tool):
    out = sh([BEDTOOLS, "intersect", "-a", gt, "-b", tool, "-wa", "-u"])
    hit = sum(1 for x in out.split("\n") if x)
    tot = regions_bp(gt)[0]
    return 100.0 * hit / tot if tot else 0.0


def region_precision(gt, tool):
    out = sh([BEDTOOLS, "intersect", "-a", tool, "-b", gt, "-wa", "-u"])
    hit = sum(1 for x in out.split("\n") if x)
    tot = regions_bp(tool)[0]
    return 100.0 * hit / tot if tot else 0.0


def bp_recall_precision(gt, tool):
    merged = sh([BEDTOOLS, "merge", "-i", tool])  # tool already sorted
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".bed") as tf:
        tf.write(merged)
        mp = tf.name
    inter = sh([BEDTOOLS, "intersect", "-a", gt, "-b", mp])
    ibp = 0
    for line in inter.split("\n"):
        if line.strip():
            p = line.split("\t")
            if len(p) >= 3:
                ibp += int(p[2]) - int(p[1])
    gt_bp = regions_bp(gt)[1]
    tool_bp = regions_bp(mp)[1]
    os.remove(mp)
    rec = 100.0 * ibp / gt_bp if gt_bp else 0.0
    prec = 100.0 * ibp / tool_bp if tool_bp else 0.0
    return rec, prec


def main():
    raw = sys.argv[1:]
    chroms = None
    args = []
    i = 0
    while i < len(raw):
        a = raw[i]
        if a == "--chroms":
            chroms = raw[i + 1].split(",")
            i += 2
            continue
        if a.startswith("--chroms="):
            chroms = a.split("=", 1)[1].split(",")
            i += 1
            continue
        args.append(a)
        i += 1

    gt_in = args[0]
    tools = args[1:]
    tmp = []
    gt = gt_in
    if chroms:
        gt = filter_chroms(gt_in, chroms)
        tmp.append(gt)
    else:
        srt = sh([BEDTOOLS, "sort", "-i", gt_in])
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".bed") as tf:
            tf.write(srt)
            gt = tf.name
        tmp.append(gt)

    gt_regions, gt_bp = regions_bp(gt)
    print(f"GT: {gt_in}  regions={gt_regions}  bp={gt_bp}"
          + (f"  chroms={','.join(chroms)}" if chroms else ""))
    print(f"{'tool':<16} {'regions':>10} {'regRecall%':>11} {'regPrec%':>9} "
          f"{'bpRecall%':>10} {'bpPrec%':>8}")
    for spec in tools:
        path, _, name = spec.partition(":")
        name = name or os.path.basename(path)
        t = path
        if chroms:
            t = filter_chroms(path, chroms)
            tmp.append(t)
        else:
            srt = sh([BEDTOOLS, "sort", "-i", path])
            with tempfile.NamedTemporaryFile("w", delete=False, suffix=".bed") as tf:
                tf.write(srt)
                t = tf.name
            tmp.append(t)
        n = regions_bp(t)[0]
        rr = region_recall(gt, t)
        rp = region_precision(gt, t)
        br, bp_ = bp_recall_precision(gt, t)
        print(f"{name:<16} {n:>10} {rr:>11.2f} {rp:>9.2f} {br:>10.2f} {bp_:>8.2f}")

    for p in tmp:
        try:
            os.remove(p)
        except OSError:
            pass


if __name__ == "__main__":
    main()
