#!/usr/bin/env python3
"""Re-score manuscript Tables 3B/3B-b/3C/3C-b under their original rule.

Coordinates are BED half-open.  Detection, fragmentation and boundary offset use
raw overlapping calls.  Coverage alone is computed after merging overlapping or
touching calls (no positive gap).  Period bands use integer BED column 5 when
available and otherwise motif length.  TRASH keeps ``none_identified`` windows,
matching the manuscript, and is parsed by the same column-5 rule.
"""
from __future__ import annotations

import copy
import datetime as dt
import hashlib
import importlib.util
import json
import os
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
SCORER = Path(os.environ.get("BWT_SCORE_EXP3",
                             str(Path(__file__).resolve().parent / "score_exp3.py")))
GT = Path(os.environ.get("BWT_MAIZE_GT_ROOT",
                         "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"))
BASE = Path(os.environ.get("BWT_COMPETITOR_BEDS",
                           "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"))
WP0 = Path(os.environ.get("BWT_WP0", "/data/gpfs/assoc/pgl/devel/exp1_human/wp0"))
REGEN = Path(os.environ.get("BWT_REGEN_MAIZE_BED",
                            "/data/gpfs/assoc/pgl/devel/exp1_human/regen/regen_maize.bed"))
OLD_BWT = BASE / "bwtandem/bwt_maize.bed"
MZ = "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic"

spec = importlib.util.spec_from_file_location("score_exp3", SCORER)
sx = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(sx)

TRF3B = BASE / f"trf/{MZ}_trf_exp3B_satellite.bed"
TRF3C = BASE / f"trf/{MZ}_trf_exp3C_centc.bed"
ULTRA3B = WP0 / "beds/ultra_maize_exp3B.bed"
ULTRA3C = WP0 / "beds/ultra_maize_exp3C.bed"
TANTAN3B = WP0 / "fixcampaign/tantan/tantan_maize_w500.bed"
TANTAN3C = WP0 / "fixcampaign/tantan/tantan_maize_w200.bed"
TRASH_DN3B = BASE / f"trash/{MZ}_trash_denovo_exp3B_satellite.bed"
# This is deliberately the knob180-template file for BOTH 3B classes: that is
# the source that reproduces the manuscript's template offsets 10086 / 8973.
TRASH_TPL3B = BASE / f"trash/{MZ}_trash_templates_knob180_exp3B_satellite.bed"
TRASH_TPL3C = BASE / f"trash/{MZ}_trash_templates_CentC_exp3C_centc.bed"

PATHS = [SCORER, GT / "mo17_knob180_arrays.bed", GT / "mo17_tr1_arrays.bed",
         GT / "mo17_centc_arrays.bed", REGEN, OLD_BWT, TRF3B, TRF3C,
         ULTRA3B, ULTRA3C, TANTAN3B, TANTAN3C, TRASH_DN3B,
         TRASH_TPL3B, TRASH_TPL3C]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path):
    rows = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            period = None
            if len(p) >= 5:
                try:
                    period = int(p[4])
                except ValueError:
                    pass
            if period is None:
                period = len(p[3])
            rows[p[0]].append((start, end, period))
    return rows


def select(rows, band):
    lo, hi = band
    return {c: [(s, e) for s, e, p in values if lo <= p <= hi]
            for c, values in rows.items()}


def unfiltered(rows):
    return {c: [(s, e) for s, e, _ in values] for c, values in rows.items()}


def merge(intervals):
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out


def coverage(by_chrom, arrays):
    merged = {c: merge(v) for c, v in by_chrom.items()}
    covered = 0
    for chrom, gs, ge in arrays:
        for start, end in merged.get(chrom, []):
            if end <= gs:
                continue
            if start >= ge:
                break
            covered += min(end, ge) - max(start, gs)
    return covered


def metrics(rows, arrays, band=None):
    calls = unfiltered(rows) if band is None else select(rows, band)
    detected, frag, offset = sx.array_metrics(arrays, copy.deepcopy(calls))
    gt_bp = sum(e - s for _, s, e in arrays)
    cov = 100.0 * coverage(calls, arrays) / gt_bp
    return {"detected": detected, "total": len(arrays), "coverage": cov,
            "frag": frag, "offset": offset}


def fmt(x):
    return f"{x:.2f}"


def evaluate(bwt_path: Path):
    knob = sx.load_gt(str(GT / "mo17_knob180_arrays.bed"))
    tr1 = sx.load_gt(str(GT / "mo17_tr1_arrays.bed"))
    centc = sx.load_gt(str(GT / "mo17_centc_arrays.bed"))
    beds3b = [("BWTandem", bwt_path), ("TRF", TRF3B), ("ULTRA", ULTRA3B),
              ("tantan-w500", TANTAN3B), ("TRASH-de-novo", TRASH_DN3B),
              ("TRASH-template", TRASH_TPL3B)]
    beds3c = [("BWTandem", bwt_path), ("TRF", TRF3C), ("ULTRA", ULTRA3C),
              ("tantan-w200", TANTAN3C), ("TRASH-template", TRASH_TPL3C)]
    result = {"3B": {}, "3C": {}}
    for tool, path in beds3b:
        rows = load(path)
        result["3B"][tool] = {
            "knob180": {"unfiltered": metrics(rows, knob),
                        "banded": metrics(rows, knob, (100, 500))},
            "TR-1": {"unfiltered": metrics(rows, tr1),
                     "banded": metrics(rows, tr1, (100, 500))},
        }
    for tool, path in beds3c:
        rows = load(path)
        result["3C"][tool] = {
            "CentC": {"unfiltered": metrics(rows, centc),
                      "banded": metrics(rows, centc, (100, 200))}}
    return result


def render_tables(result):
    lines = []
    lines += ["## Table 3B replacement (accuracy columns)", "",
              "| Tool | knob180 | TR-1 | knob180 Coverage (%) | TR-1 Coverage (%) | knob180 Frag. | knob180 Offset (bp) | TR-1 Offset (bp) |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for tool, d in result["3B"].items():
        k, t = d["knob180"]["unfiltered"], d["TR-1"]["unfiltered"]
        lines.append(f"| {tool} | {k['detected']}/{k['total']} | {t['detected']}/{t['total']} | {fmt(k['coverage'])} | {fmt(t['coverage'])} | {k['frag']:.2f} | {k['offset']:.2f} | {t['offset']:.2f} |")
    lines += ["", "## Table 3B-b replacement", "",
              "| Class | Tool | Rule | Detected | Coverage (%) | Offset (bp) | Loss under band (pp) |",
              "|---|---|---|---:|---:|---:|---:|"]
    for cls in ("knob180", "TR-1"):
        for tool, d in result["3B"].items():
            u, b = d[cls]["unfiltered"], d[cls]["banded"]
            lines.append(f"| {cls} | {tool} | unfiltered | {u['detected']}/{u['total']} | {fmt(u['coverage'])} | {u['offset']:.2f} |  |")
            lines.append(f"| {cls} | {tool} | banded | {b['detected']}/{b['total']} | {fmt(b['coverage'])} | {b['offset']:.2f} | {b['coverage']-u['coverage']:+.2f} |")
    lines += ["", "## Table 3C replacement (accuracy columns)", "",
              "| Tool | Detected | Coverage (%) | Frag. | Offset (bp) |",
              "|---|---:|---:|---:|---:|"]
    for tool, d in result["3C"].items():
        x = d["CentC"]["unfiltered"]
        lines.append(f"| {tool} | {x['detected']}/{x['total']} | {fmt(x['coverage'])} | {x['frag']:.2f} | {x['offset']:.2f} |")
    lines += ["", "## Table 3C-b replacement", "",
              "| Tool | Rule | Detected | Coverage (%) | Offset (bp) | Loss under band (pp) |",
              "|---|---|---:|---:|---:|---:|"]
    for tool, d in result["3C"].items():
        u, b = d["CentC"]["unfiltered"], d["CentC"]["banded"]
        lines.append(f"| {tool} | unfiltered | {u['detected']}/{u['total']} | {fmt(u['coverage'])} | {u['offset']:.2f} |  |")
        lines.append(f"| {tool} | banded | {b['detected']}/{b['total']} | {fmt(b['coverage'])} | {b['offset']:.2f} | {b['coverage']-u['coverage']:+.2f} |")
    return "\n".join(lines) + "\n"


def main():
    missing = [str(p) for p in PATHS if not p.exists()]
    if missing:
        raise SystemExit("Missing inputs:\n" + "\n".join(missing))
    old = evaluate(OLD_BWT)
    new = evaluate(REGEN)
    # Hard validation against manuscript values. This prevents silently adopting
    # the postmerge scorer's different convention.
    checks = [
        (old["3B"]["TRF"]["knob180"]["unfiltered"]["coverage"], 80.01),
        (old["3B"]["TRF"]["knob180"]["unfiltered"]["offset"], 651.96),
        (old["3B"]["TRASH-template"]["TR-1"]["unfiltered"]["offset"], 8973.15),
        (old["3C"]["TRF"]["CentC"]["unfiltered"]["coverage"], 58.50),
        (old["3C"]["TRASH-template"]["CentC"]["unfiltered"]["coverage"], 58.68),
    ]
    failures = [(got, expected) for got, expected in checks
                if round(got, 2) != round(expected, 2)]
    if failures:
        raise SystemExit(f"Historical-convention validation failed: {failures}")
    provenance = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "convention": {
            "coordinates": "BED half-open",
            "detection_fragmentation_offset": "raw calls; any >=1 bp overlap",
            "coverage": "pooled bp; overlap/touch merge only; no positive gap",
            "period": "integer column 5, else len(column 4 motif)",
            "trash": "fixed windows retained including none_identified",
        },
        "inputs": [{"path": str(p), "size": p.stat().st_size,
                    "mtime": dt.datetime.fromtimestamp(p.stat().st_mtime,
                        dt.timezone.utc).isoformat(), "sha256": sha256(p)}
                   for p in PATHS],
        "historical_validation": "PASS",
    }
    (OUT / "table3bc_results.json").write_text(json.dumps(
        {"old": old, "regen": new}, indent=2) + "\n")
    (OUT / "table3bc_replacement.md").write_text(render_tables(new))
    (OUT / "table3bc_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(render_tables(new), end="")


if __name__ == "__main__":
    main()
