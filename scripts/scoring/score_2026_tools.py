#!/usr/bin/env python3
"""Score the 2026 tools (longdust, AniAnn's) with the FROZEN scorers (#16).

Same discipline as the regen-safe wrappers: the published rows' code paths
are reused verbatim by loading the frozen modules and overriding only their
source lists — no second scoring implementation. Protocol:
docs/2026-09-01-longdust-anianns-benchmark-protocol.md.

Modes:
  score_2026_tools.py human  NAME:BED [NAME:BED ...]
      -> frozen score_table1.py over the adotto GT with SOURCES = the given
         tools only (existing deposited rows are untouched).
  score_2026_tools.py maize  NAME:BED [NAME:BED ...]
      -> frozen score_maize_postmerge.py with TOOLS = the given tools.
  (Col-CEN needs no wrapper: score_colcen.py already accepts LABEL:PATH.)
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def need(path, label):
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"FATAL: {label} missing or empty: {path}")
    return path


def parse_tools(args):
    tools = []
    for spec in args:
        name, _, path = spec.partition(":")
        if not path:
            raise SystemExit(f"FATAL: expected NAME:BED, got {spec!r}")
        tools.append((name, need(path, name)))
    if not tools:
        raise SystemExit("FATAL: no tools given")
    return tools


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ("human", "maize"):
        raise SystemExit(__doc__)
    mode, tools = sys.argv[1], parse_tools(sys.argv[2:])

    if mode == "human":
        gt = need(os.environ.get(
            "BWT_TABLE1_GT",
            "/data/gpfs/assoc/pgl/devel/exp1_human/data/adotto_primary.bed"),
            "ground truth")
        src = os.path.join(HERE, "score_table1.py")
        spec = importlib.util.spec_from_file_location("_frozen_table1", src)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.GT, mod.SOURCES, mod.EXTRA = gt, tools, []
        sys.argv = [sys.argv[0]]
        mod.main()
    else:
        gt_root = os.environ.get(
            "BWT_MAIZE_GT_ROOT",
            "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth")
        scorer = need(os.environ.get(
            "BWT_SCORE_EXP3",
            "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"),
            "score_exp3 module")
        for fn in ("mo17_knob180_arrays.bed", "mo17_tr1_arrays.bed",
                   "mo17_centc_arrays.bed"):
            need(os.path.join(gt_root, fn), fn)
        src = os.path.join(HERE, "score_maize_postmerge.py")
        text = open(src).read().replace(
            'SCORER = "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"',
            f"SCORER = {scorer!r}", 1)
        mod = type(sys)("_frozen_maize_postmerge")
        mod.__file__ = src
        exec(compile(text, src, "exec"), mod.__dict__)
        mod.GT, mod.TOOLS = os.path.abspath(gt_root), tools
        sys.argv = [sys.argv[0]]
        mod.main()


if __name__ == "__main__":
    main()
