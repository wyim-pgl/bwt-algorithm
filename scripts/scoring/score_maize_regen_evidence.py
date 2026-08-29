#!/usr/bin/env python3
"""Emit deposited evidence for regenerated maize manuscript values.

This scorer covers Table 3A/3A-b, the coordinate-only post-merge sweep, and
the 100--200-bp whole-genome call totals quoted in Section 3.3.3.  It does not
change any detection output.  Period is integer BED column 5 when available,
otherwise motif length (BWTandem column 5 is a decimal copy count).
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import shlex
import sys
from collections import defaultdict
from pathlib import Path

TAG_ROT = {"TAG", "AGT", "GTA", "CTA", "TAC", "ACT"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def period_of(parts):
    if len(parts) >= 5:
        try:
            return int(parts[4])
        except ValueError:
            pass
    return len(parts[3]) if len(parts) >= 4 else None


def load_calls(path):
    by_chrom = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            by_chrom[p[0]].append((start, end, p[3], period_of(p)))
    for calls in by_chrom.values():
        calls.sort()
    return by_chrom


def merge_coords(intervals, gap=0):
    out = []
    for start, end in sorted(intervals):
        if out and start - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(start, end) for start, end in out]


def scan_3a(calls, band=None):
    n = tag_sub = tag_fair = longest_sub = longest_fair = 0
    intervals = defaultdict(list)
    for chrom, rows in calls.items():
        for start, end, motif, period in rows:
            if band and (period is None or not band[0] <= period <= band[1]):
                continue
            n += 1
            intervals[chrom].append((start, end))
            motif = motif.upper()
            span = end - start
            if "TAG" in motif or "CTA" in motif:
                tag_sub += 1
                longest_sub = max(longest_sub, span)
            if period == 3 and motif[:3] in TAG_ROT:
                tag_fair += 1
                longest_fair = max(longest_fair, span)
    bp = sum(e - s for rows in intervals.values() for s, e in merge_coords(rows))
    return {"bp": bp, "regions": n, "tag_substring": tag_sub,
            "longest_substring_bp": longest_sub, "tag_rotation": tag_fair,
            "longest_rotation_bp": longest_fair}


def coverage(merged, arrays):
    total = 0
    for chrom, gs, ge in arrays:
        for start, end in merged.get(chrom, []):
            if end <= gs:
                continue
            if start >= ge:
                break
            total += min(end, ge) - max(start, gs)
    return total


def calls_per_array(merged, arrays):
    counts = []
    for chrom, gs, ge in arrays:
        n = sum(start < ge and end > gs for start, end in merged.get(chrom, []))
        if n:
            counts.append(n)
    return sum(counts) / len(counts) if counts else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bwt", required=True)
    ap.add_argument("--trf-3a", required=True)
    ap.add_argument("--ultra-3a", required=True)
    ap.add_argument("--tantan-3a", required=True)
    ap.add_argument("--trf-3b", required=True)
    ap.add_argument("--ultra-3b", required=True)
    ap.add_argument("--tantan-500", required=True)
    ap.add_argument("--trf-3c", required=True)
    ap.add_argument("--ultra-3c", required=True)
    ap.add_argument("--tantan-200", required=True)
    ap.add_argument("--knob-gt", required=True)
    ap.add_argument("--tr1-gt", required=True)
    ap.add_argument("--centc-gt", required=True)
    ap.add_argument("--score-exp3", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    paths = {key: str(Path(value).resolve()) for key, value in vars(args).items()
             if key not in {"output"}}
    for label, path in paths.items():
        if not Path(path).is_file():
            ap.error(f"missing {label}: {path}")

    calls = {label: load_calls(path) for label, path in paths.items()
             if label not in {"knob_gt", "tr1_gt", "centc_gt", "score_exp3"}}
    table3a = {}
    for label in ("bwt", "trf_3a", "ultra_3a", "tantan_3a"):
        table3a[label] = {"published_rule": scan_3a(
            calls[label], (1, 6) if label == "bwt" else None),
            "band_1_6": scan_3a(calls[label], (1, 6))}

    spec = importlib.util.spec_from_file_location("score_exp3", paths["score_exp3"])
    score_exp3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(score_exp3)
    arrays = {"knob180": score_exp3.load_gt(paths["knob_gt"]),
              "TR-1": score_exp3.load_gt(paths["tr1_gt"]),
              "CentC": score_exp3.load_gt(paths["centc_gt"])}
    postmerge = {}
    # The regenerated coordinate-only prose reports BWTandem only.  Do not
    # spend hours sweeping inherited comparator BEDs that feed no quoted cell.
    post_sources = {"BWTandem": "bwt"}
    for cls, truth in arrays.items():
        gt_bp = sum(end - start for _, start, end in truth)
        postmerge[cls] = {}
        for tool, source in post_sources.items():
            postmerge[cls][tool] = {}
            for gap in (0, 100, 500, 1000, 2000, 5000, 10000):
                merged = {chrom: merge_coords([(s, e) for s, e, _, _ in rows], gap)
                          for chrom, rows in calls[source].items()}
                detected, frag, offset = score_exp3.array_metrics(
                    truth, copy.deepcopy(merged))
                postmerge[cls][tool][str(gap)] = {
                    "detected": detected, "total_arrays": len(truth),
                    "calls_per_array": round(calls_per_array(merged, truth), 1),
                    "fragmentation": frag,
                    "coverage_percent": round(100 * coverage(merged, truth) / gt_bp, 2),
                    "mean_offset_bp": round(offset),
                }

    band_sources = {"BWTandem": "bwt", "TRF": "trf_3c",
                    "ULTRA": "ultra_3c", "tantan-w200": "tantan_200"}
    band_totals = {}
    for tool, source in band_sources.items():
        intervals = {chrom: [(s, e) for s, e, _, p in rows if 100 <= p <= 200]
                     for chrom, rows in calls[source].items()}
        band_totals[tool] = sum(e - s for rows in intervals.values()
                                for s, e in merge_coords(rows))

    result = {
        "command": " ".join(shlex.quote(x) for x in ["python3", *sys.argv]),
        "conventions": {
            "coordinates": "BED 0-based half-open",
            "period": "integer column 5, else motif length",
            "total_bp": "overlap/touch merged within chromosome",
            "postmerge": "coordinate-only; join when gap <= threshold; motif ignored",
        },
        "inputs": {label: {"path": path, "sha256": sha256(path),
                            "bytes": Path(path).stat().st_size}
                   for label, path in paths.items()},
        "table3a": table3a,
        "coordinate_postmerge": postmerge,
        "period_100_200_merged_bp": band_totals,
    }
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
