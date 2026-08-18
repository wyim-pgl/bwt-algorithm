#!/usr/bin/env python3
"""WP0.1/0.3 -- rescore maize 3B and 3C for every tool under BOTH scoring rules.

Established in test_filter_hypothesis.py: the published Table 3C scored TRF and
tantan with no motif-length filter (reproducing 69.62 and 9244 bp exactly) but
scored BWTandem through the band-filtered path (reproducing 3051 exactly). The
two rules give BWTandem 59.02 % vs 41.19 % coverage and 212 vs 3051 bp offset, so
the mixed table understates BWTandem by the width of that gap.

Both rules are defensible -- "did the tool cover the array with anything" vs
"did the tool call it as a monomer-period satellite" -- but the table has to pick
one and apply it to every tool. This prints both, for every tool, for 3B and 3C.
"""
import copy
import importlib.util
import os
from collections import defaultdict

SCORER = "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"
GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
B = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"
HERE = os.path.dirname(os.path.abspath(__file__))
MZ = "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic"

spec = importlib.util.spec_from_file_location("score_exp3", SCORER)
sx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sx)

# The three published maize tantan BEDs (exp3A/3B/3C) are byte-identical copies
# of one genome-wide run whose maximum period is exactly 100 -- tantan's default
# -w 100 window. Neither the 3B band (100-500) nor the 3C band (100-200) can
# therefore hold anything but a single period-100 point, so those rows measure
# the window we passed rather than the tool. TANTAN_MAIZE_BED swaps in a re-run
# at a wider window; unset, the published BED is used and the table reproduces.
TANTAN_3B = os.environ.get("TANTAN_MAIZE_BED",
                           f"{B}/tantan/{MZ}_tantan_exp3B_satellite.bed")
TANTAN_3C = os.environ.get("TANTAN_MAIZE_BED",
                           f"{B}/tantan/{MZ}_tantan_exp3C_centc.bed")

BEDS_3B = [
    ("TRF",      f"{B}/trf/{MZ}_trf_exp3B_satellite.bed"),
    ("ULTRA",    f"{HERE}/beds/ultra_maize_exp3B.bed"),
    ("BWTandem", os.environ.get("BWT_MAIZE_BED", f"{B}/bwtandem/bwt_maize.bed")),
    ("tantan",   TANTAN_3B),
]
BEDS_3C = [
    ("TRF",      f"{B}/trf/{MZ}_trf_exp3C_centc.bed"),
    ("ULTRA",    f"{HERE}/beds/ultra_maize_exp3C.bed"),
    ("BWTandem", os.environ.get("BWT_MAIZE_BED", f"{B}/bwtandem/bwt_maize.bed")),
    ("tantan",   TANTAN_3C),
]


def load(path, lo=None, hi=None):
    by_chrom = defaultdict(list)
    n = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            if lo is not None:
                period = None
                if len(p) >= 5:
                    try:
                        period = int(p[4])
                    except ValueError:
                        period = None
                if period is None:
                    period = len(p[3])
                if not (lo <= period <= hi):
                    continue
            by_chrom[p[0]].append((s, e))
            n += 1
    return by_chrom, n


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def coverage(by_chrom, arrays):
    cov = 0
    merged = {c: merge(v) for c, v in by_chrom.items()}
    for chrom, gs, ge in arrays:
        for s, e in merged.get(chrom, []):
            if e <= gs:
                continue
            if s >= ge:
                break
            cov += min(e, ge) - max(s, gs)
    return cov


def block(title, beds, arrays_list, band):
    print(f"\n{'='*78}\n{title}\n{'='*78}")
    for gt_label, arrays in arrays_list:
        gt_bp = sum(e - s for _, s, e in arrays)
        print(f"\n  [{gt_label}]  {len(arrays)} arrays, {gt_bp:,} bp")
        print(f"  {'tool':10s} {'rule':13s} {'det':>8s} {'cov%':>8s} "
              f"{'offset':>11s} {'frag':>6s}")
        for label, bed in beds:
            for tag, lo, hi in (("unfiltered", None, None),
                                (f"band {band[0]}-{band[1]}", band[0], band[1])):
                by_chrom, _ = load(bed, lo, hi)
                det, frag, off = sx.array_metrics(arrays, copy.deepcopy(by_chrom))
                cov = coverage(by_chrom, arrays)
                print(f"  {label:10s} {tag:13s} {det:>4d}/{len(arrays):<3d} "
                      f"{100.0*cov/gt_bp:>7.2f}% {off:>11.2f} {frag:>6.2f}")


def main():
    knob = sx.load_gt(f"{GT}/mo17_knob180_arrays.bed")
    tr1 = sx.load_gt(f"{GT}/mo17_tr1_arrays.bed")
    centc = sx.load_gt(f"{GT}/mo17_centc_arrays.bed")
    block("TABLE 3B -- knob180 / TR-1", BEDS_3B,
          [("knob180", knob), ("TR-1", tr1)], (100, 500))
    block("TABLE 3C -- CentC", BEDS_3C, [("CentC", centc)], (100, 200))


if __name__ == "__main__":
    main()
