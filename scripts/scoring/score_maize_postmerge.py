#!/usr/bin/env python3
"""REV-38 -- post-merge accuracy on the maize arrays, and why the built-in merge
does not fire on them.

Section 3.3.2 reports that every tool except TRASH splits each curated array into
hundreds of calls and closes with "an analysis needing one interval per array
will have to merge". The reviewer asks for the accuracy of a stated merge, and
whether period-assignment drift is what defeats the pipeline's own merge.

The built-in merge (finder.py _merge_adjacent_repeats) joins two calls only when
their canonical motifs are equal, with a 10% Hamming tolerance that applies only
when the two motifs have equal length and are at least 50 bp. Its gap allowance
is generous on satellites, 100 times the period for periods of 100 bp or more.
So if adjacent calls on one array carry motifs of different lengths, the gap is
irrelevant: the equality test fails, the tolerance path is not reachable, and no
merge happens however close the calls are. That is the hypothesis this measures.

Part 1 is the diagnostic: among each tool's calls inside a curated array, how
many distinct reported periods and motif lengths appear, and what fraction of
adjacent pairs agree closely enough for the built-in rule to fire.

Part 2 is the stated merge: coordinates only, motif ignored, calls joined when
the gap is at most G. Coverage, boundary offset and fragmentation are recomputed
through score_exp3.array_metrics, the same function behind Tables 3B and 3C, so
the post-merge rows are differenceable against the published ones. G is swept so
the reader can see the trade rather than one chosen point.

Usage: score_maize_postmerge.py [G ...]     (default sweep: 0 100 1000 10000)
"""
import copy
import importlib.util
import os
import sys
from collections import defaultdict

SCORER = "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"
GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
B = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"
WP0 = "/data/gpfs/assoc/pgl/devel/exp1_human/wp0"
MZ = "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic"

spec = importlib.util.spec_from_file_location("score_exp3", SCORER)
sx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sx)

TOOLS = [
    ("BWTandem", f"{B}/bwtandem/bwt_maize.bed"),
    ("TRF", f"{B}/trf/{MZ}_trf_exp3B_satellite.bed"),
    ("ULTRA", f"{WP0}/beds/ultra_maize_exp3B.bed"),
    ("tantan (500bp)",
     f"{WP0}/fixcampaign/tantan/tantan_maize_w500.bed"),
]


def load(path):
    """(start, end, motif, period) per chromosome."""
    by_chrom = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            motif = p[3]
            period = None
            if len(p) >= 5:
                try:
                    period = int(p[4])
                except ValueError:
                    period = None
            if period is None:
                period = len(motif)
            by_chrom[p[0]].append((s, e, motif, period))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def merge_coords(intervals, gap):
    """Coordinate-only merge: join when the gap is at most `gap`."""
    out = []
    for s, e in sorted(intervals):
        if out and s - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def coverage(merged, arrays):
    """`merged` is {chrom: [(start, end), ...]}, already merged and sorted."""
    cov = 0
    for chrom, gs, ge in arrays:
        for s, e in merged.get(chrom, []):
            if e <= gs:
                continue
            if s >= ge:
                break
            cov += min(e, ge) - max(s, gs)
    return cov


def mean_calls_per_array(merged, arrays):
    """Mean intervals overlapping an array, over the arrays that have any.

    The quantity score_exp3 reports as a rounded reciprocal; computed here at
    full precision so the post-merge sweep can show hundreds falling to units.
    """
    counts = []
    for chrom, gs, ge in arrays:
        n = sum(1 for s, e in merged.get(chrom, []) if s < ge and e > gs)
        if n:
            counts.append(n)
    return sum(counts) / len(counts) if counts else 0.0


def in_array(by_chrom, arrays):
    """Calls overlapping any curated array, per array."""
    per = []
    for chrom, gs, ge in arrays:
        hits = [c for c in by_chrom.get(chrom, []) if c[0] < ge and c[1] > gs]
        if hits:
            per.append(hits)
    return per


def builtin_would_fire(a, b):
    """Does finder.py's motif test pass for this adjacent pair?

    Canonicalisation is strand-aware in the pipeline; here the comparison is on
    the reported motif itself, which is an upper bound on agreement -- if the raw
    motifs already differ in length the canonical forms cannot be equal and the
    Hamming path is closed, which is the case the diagnostic is after.
    """
    ma, mb = a[2], b[2]
    if ma == mb:
        return True
    if len(ma) == len(mb) and len(ma) >= 50:
        ham = sum(1 for x, y in zip(ma, mb) if x != y)
        return ham / len(ma) <= 0.10
    return False


def diagnose(label, by_chrom, arrays, cls):
    per = in_array(by_chrom, arrays)
    if not per:
        print(f"  {label:16s} no calls inside any {cls} array")
        return
    ncalls = sum(len(h) for h in per)
    periods, lengths, pairs_ok, pairs = set(), set(), 0, 0
    for hits in per:
        for h in hits:
            periods.add(h[3])
            lengths.add(len(h[2]))
        for i in range(len(hits) - 1):
            pairs += 1
            if builtin_would_fire(hits[i], hits[i + 1]):
                pairs_ok += 1
    print(f"  {label:16s} {ncalls:>9,d} calls in {len(per):>2d} arrays  "
          f"distinct periods {len(periods):>5d}  distinct motif lengths "
          f"{len(lengths):>5d}  adjacent pairs the built-in merge would join "
          f"{100.0*pairs_ok/pairs if pairs else 0:>6.2f}%")


def main():
    gaps = [int(g) for g in sys.argv[1:]] or [0, 100, 1000, 10000]
    knob = sx.load_gt(f"{GT}/mo17_knob180_arrays.bed")
    tr1 = sx.load_gt(f"{GT}/mo17_tr1_arrays.bed")
    centc = sx.load_gt(f"{GT}/mo17_centc_arrays.bed")

    loaded = {}
    for label, path in TOOLS:
        if not os.path.exists(path):
            print(f"MISSING {label}: {path}")
            continue
        loaded[label] = load(path)

    for cls, arrays in (("knob180", knob), ("TR-1", tr1), ("CentC", centc)):
        gt_bp = sum(e - s for _, s, e in arrays)
        print(f"\n{'='*100}\n{cls}  --  {len(arrays)} arrays, {gt_bp:,} bp\n{'='*100}")
        print("\n [why the built-in merge does not fire]")
        for label in loaded:
            diagnose(label, loaded[label], arrays, cls)

        print(f"\n [stated post-merge, coordinates only]")
        print(f"  {'tool':16s} {'gap':>7s} {'det':>7s} {'calls/array':>12s} "
              f"{'frag':>6s} {'cov%':>8s} {'offset':>13s}")
        for label in loaded:
            for g in gaps:
                merged = {c: merge_coords([(s, e) for s, e, _, _ in v], g)
                          for c, v in loaded[label].items()}
                det, frag, off = sx.array_metrics(arrays, copy.deepcopy(merged))
                cov = coverage(merged, arrays)
                # Calls per array direct, not as 1/frag: frag is rounded to two
                # decimals, so inverting it turns any array split into more than
                # a hundred calls into a division by zero. This is the same
                # precision loss the Table 3B caption warns about.
                cpa = mean_calls_per_array(merged, arrays)
                print(f"  {label:16s} {g:>7,d} {det:>3d}/{len(arrays):<3d} "
                      f"{cpa:>12.1f} {frag:>6.2f} {100.0*cov/gt_bp:>7.2f}% "
                      f"{off:>13,.0f}")


if __name__ == "__main__":
    main()
