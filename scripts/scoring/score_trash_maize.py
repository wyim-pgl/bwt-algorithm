#!/usr/bin/env python3
"""WP0 -- rescore TRASH on maize.

TRASH's BED differs from every other tool's in two ways that affect scoring:

1. Its intervals are fixed analysis WINDOWS (0-24999, 728000-796999, ...), not
   repeat arrays, so a window overlapping a ground-truth array counts as a
   detection regardless of what TRASH actually found in it.
2. 14-18 % of rows are `none_identified` with period `N/A` -- windows where TRASH
   identified no repeat at all. They are still intervals in the BED, so
   bedtools-based detection and coverage count them as calls.

Like TRF, TRASH puts a full consensus (up to ~kb) in the motif column with the
period in column 4, so any len(motif)-derived band filter mis-handles it.

This scores every TRASH maize variant under four rules: as-published
(unfiltered, none_identified kept), none_identified dropped, band-filtered, and
both.
"""
import copy
import importlib.util
from collections import defaultdict

SCORER = "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"
GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
T = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds/trash"
MZ = "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic"

spec = importlib.util.spec_from_file_location("score_exp3", SCORER)
sx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sx)


def load(path, lo=None, hi=None, drop_none=False):
    by_chrom = defaultdict(list)
    n = n_none = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            is_none = (p[3] == "none_identified")
            if is_none:
                n_none += 1
                if drop_none:
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
    return by_chrom, n, n_none


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


RULES = [
    ("as published", None, None, False),
    ("drop none_id", None, None, True),
    ("band",         "band", None, False),
    ("band+drop",    "band", None, True),
]


def run(title, beds, arrays, band, published):
    gt_bp = sum(e - s for _, s, e in arrays)
    print(f"\n{'='*80}\n{title}   ({len(arrays)} arrays, {gt_bp:,} bp)"
          f"   published: {published}\n{'='*80}")
    print(f"  {'variant':22s} {'rule':13s} {'calls':>8s} {'det':>8s} "
          f"{'cov%':>8s} {'offset':>11s} {'frag':>6s}")
    for label, bed in beds:
        for tag, mode, _, drop in RULES:
            lo, hi = (band if mode == "band" else (None, None))
            by_chrom, n, n_none = load(bed, lo, hi, drop)
            det, frag, off = sx.array_metrics(arrays, copy.deepcopy(by_chrom))
            cov = coverage(by_chrom, arrays)
            print(f"  {label:22s} {tag:13s} {n:>8,d} {det:>5d}/{len(arrays):<2d} "
                  f"{100.0*cov/gt_bp:>7.2f}% {off:>11.2f} {frag:>6.2f}")


def main():
    knob = sx.load_gt(f"{GT}/mo17_knob180_arrays.bed")
    tr1 = sx.load_gt(f"{GT}/mo17_tr1_arrays.bed")
    centc = sx.load_gt(f"{GT}/mo17_centc_arrays.bed")

    b3b = [("de novo", f"{T}/{MZ}_trash_denovo_exp3B_satellite.bed"),
           ("template knob180", f"{T}/{MZ}_trash_templates_knob180_exp3B_satellite.bed"),
           ("template TR1", f"{T}/{MZ}_trash_templates_TR1_exp3B_satellite.bed"),
           ("merged exp3B", f"{T}/{MZ}_trash_exp3B.bed")]
    b3c = [("template CentC", f"{T}/{MZ}_trash_templates_CentC_exp3C_centc.bed"),
           ("merged exp3C", f"{T}/{MZ}_trash_exp3C.bed")]

    run("TRASH knob180", b3b, knob, (100, 500), "25/25 frag 0.05 off 9845 (dn) / 10086 (tpl)")
    run("TRASH TR-1", b3b, tr1, (100, 500), "17/17 off 10561 (dn) / 8973 (tpl)")
    run("TRASH CentC", b3c, centc, (100, 200), "17/17 cov 58.33 frag 0.05 off 5319")


if __name__ == "__main__":
    main()
