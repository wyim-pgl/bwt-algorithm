#!/usr/bin/env python3
"""Untracked regen-safe wrapper around the frozen Table-1 scorer at 0363d8b."""
import argparse, importlib.util, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
def need(path, label):
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"FATAL: {label} missing or empty: {path}")
    return path
def main():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--bwt-bed", required=True)
    p.add_argument("--gt", default=os.environ.get("BWT_TABLE1_GT", "/data/gpfs/assoc/pgl/devel/exp1_human/data/adotto_primary.bed"))
    p.add_argument("--competitor-root", default=os.environ.get("BWT_COMPETITOR_ROOT", "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"))
    p.add_argument("--ultra-bed", default=os.environ.get("BWT_ULTRA_HUMAN_BED", "/data/gpfs/assoc/pgl/devel/exp1_human/wp0/beds/ultra_human_GCA.bed"))
    args, rest = p.parse_known_args()
    root = os.path.abspath(args.competitor_root)
    sources = [("BWTandem", args.bwt_bed), ("TRF", os.path.join(root, "trf/GCA_000001405.15_GRCh38_genomic_output.bed")), ("ULTRA", args.ultra_bed), ("tantan", os.path.join(root, "tantan/GCA_000001405.15_GRCh38_genomic_output.bed")), ("TRASH", os.path.join(root, "trash/GCA_000001405.15_GRCh38_genomic_trash.bed"))]
    sources = [(name, need(path, name)) for name, path in sources]
    original = os.environ.get("BWT_TABLE1_ORIGINAL", "score_table1.py")
    spec = importlib.util.spec_from_file_location("_frozen_table1", os.path.join(HERE, original))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.GT, mod.SOURCES, mod.EXTRA = need(args.gt, "ground truth"), sources, []
    sys.argv = [sys.argv[0], *rest]
    mod.main()
if __name__ == "__main__": main()
