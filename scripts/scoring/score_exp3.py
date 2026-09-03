#!/usr/bin/env python3
"""Score Experiment 3 (Maize Mo17) for bwtandem_v2, replicating Filip's
compute_metrics.py semantics in pure Python (no bedtools dependency).

3A microsatellite : total regions + bp of motif-length 1-6 bp calls
3B satellite      : knob180/TR-1 array detection + norm-frag + boundary offset
                    (motif-length 100-500 bp calls)
3C centc          : CentC array detection + norm-frag + boundary offset
                    (motif-length 100-200 bp calls)

Validates the implementation against the OLD bwtandem published row first.
"""
import os
import sys
from collections import defaultdict

GT = os.environ.get("BWT_MAIZE_GT_ROOT",
                    "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth")
OLD_B = os.environ.get("BWT_MAIZE_OLD_BEDS",
                       "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds/bwtandem")
V2 = os.environ.get("BWT_MAIZE_V2_BED",
                    "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/out/bwt_maize_v2.bed")


def load_gt(path):
    arrays = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                arrays.append((p[0], int(p[1]), int(p[2])))
    return arrays


def load_calls_by_chrom(path, motif_col, min_len=None, max_len=None):
    """Load (start,end) calls bucketed by chrom; filter by motif length window.
    motif_col: 0-based index of motif column used to derive period (=len)."""
    by_chrom = defaultdict(list)
    n = 0
    bp = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= motif_col:
                continue
            try:
                s = int(p[1]); e = int(p[2])
            except ValueError:
                continue
            L = len(p[motif_col])
            if min_len is not None and L < min_len:
                continue
            if max_len is not None and L > max_len:
                continue
            by_chrom[p[0]].append((s, e))
            n += 1
            bp += (e - s)
    return by_chrom, n, bp


def overlaps(a_s, a_e, b_s, b_e):
    return max(a_s, b_s) < min(a_e, b_e)


def array_metrics(arrays, by_chrom):
    """Replicate compute_array_detection / norm_fragmentation / boundary_offset."""
    # sort calls per chrom by start for slight speed
    for c in by_chrom:
        by_chrom[c].sort()
    detected = 0
    frags = []          # fragments per overlapped array (count of overlapping calls)
    offsets = []        # boundary offset per overlapped array
    for (chrom, gs, ge) in arrays:
        calls = by_chrom.get(chrom, [])
        hits = [(s, e) for (s, e) in calls if overlaps(gs, ge, s, e)]
        cnt = len(hits)
        if cnt > 0:
            detected += 1
            frags.append(cnt)
            min_ts = min(s for s, e in hits)
            max_te = max(e for s, e in hits)
            offsets.append((abs(gs - min_ts) + abs(ge - max_te)) / 2.0)
    mean_frag = sum(frags) / len(frags) if frags else 0
    norm_frag = round(1.0 / mean_frag, 2) if mean_frag > 0 else 0
    mean_off = round(sum(offsets) / len(offsets), 2) if offsets else 0
    return detected, norm_frag, mean_off


def score(label, bed, motif_col):
    knob = load_gt(f"{GT}/mo17_knob180_arrays.bed")
    tr1 = load_gt(f"{GT}/mo17_tr1_arrays.bed")
    centc = load_gt(f"{GT}/mo17_centc_arrays.bed")

    # 3A microsatellite: motif len 1-6
    _, n3a, bp3a = load_calls_by_chrom(bed, motif_col, 1, 6)

    # 3B satellite: motif len 100-500
    b3b, _, _ = load_calls_by_chrom(bed, motif_col, 100, 500)
    # need a fresh dict per GT set (array_metrics sorts in place; reuse ok)
    import copy
    k_det, k_frag, k_off = array_metrics(knob, copy.deepcopy(b3b))
    t_det, t_frag, t_off = array_metrics(tr1, copy.deepcopy(b3b))

    # 3C centc: motif len 100-200
    b3c, _, _ = load_calls_by_chrom(bed, motif_col, 100, 200)
    c_det, c_frag, c_off = array_metrics(centc, copy.deepcopy(b3c))

    print(f"\n===== {label} =====")
    print(f"3A microsat : regions={n3a:,}  total_bp={bp3a:,}")
    print(f"3B knob180  : {k_det}/25  norm_frag={k_frag}  offset={k_off}")
    print(f"3B TR-1     : {t_det}/17")
    print(f"3C CentC    : {c_det}/17  norm_frag={c_frag}  offset={c_off}")
    return dict(n3a=n3a, bp3a=bp3a, k_det=k_det, k_frag=k_frag, k_off=k_off,
                t_det=t_det, c_det=c_det, c_frag=c_frag, c_off=c_off)


if __name__ == "__main__":
    # --- VALIDATION: old bwtandem standardized beds (col3=motif idx) ---
    # standardized format: chrom start end motif period BWTandem  -> motif_col=3
    import glob
    old_3a = glob.glob(f"{OLD_B}/*exp3A_microsatellite.bed")[0]
    old_3b = glob.glob(f"{OLD_B}/*exp3B_satellite.bed")[0]
    old_3c = glob.glob(f"{OLD_B}/*exp3C_centc.bed")[0]
    # old beds are pre-split; score each directly with no length filter
    knob = load_gt(f"{GT}/mo17_knob180_arrays.bed")
    tr1 = load_gt(f"{GT}/mo17_tr1_arrays.bed")
    centc = load_gt(f"{GT}/mo17_centc_arrays.bed")
    import copy

    _, n3a, bp3a = load_calls_by_chrom(old_3a, 3)
    b3b, _, _ = load_calls_by_chrom(old_3b, 3)
    k = array_metrics(knob, copy.deepcopy(b3b)); t = array_metrics(tr1, copy.deepcopy(b3b))
    b3c, _, _ = load_calls_by_chrom(old_3c, 3)
    c = array_metrics(centc, copy.deepcopy(b3c))
    print("===== VALIDATION: OLD bwtandem (recompute vs published) =====")
    print(f"3A regions={n3a:,} bp={bp3a:,}   [published 54,283 / 3,044,731]")
    print(f"3B knob180={k[0]}/25 frag={k[1]} off={k[2]}   [published 25 / 0 / 985.06]")
    print(f"3B TR-1={t[0]}/17   [published 16]")
    print(f"3C CentC={c[0]}/17 frag={c[1]} off={c[2]}   [published 17 / 0.01 / 4837.53]")

    # --- v2 raw bed: motif_col=3 (chrom start end motif copies tier mm strand) ---
    score("bwtandem_v2 (maize, motif-len filtered from single 1-2000 run)", V2, 3)
