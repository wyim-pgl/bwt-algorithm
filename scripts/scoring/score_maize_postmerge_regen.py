#!/usr/bin/env python3
"""Untracked regen-safe wrapper around score_maize_postmerge.py at 0363d8b."""
import argparse, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
def need(path, label):
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"FATAL: {label} missing or empty: {path}")
    return path
def main():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--bwt-bed", required=True)
    p.add_argument("--gt-root", default=os.environ.get("BWT_MAIZE_GT_ROOT", "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"))
    p.add_argument("--competitor-root", default=os.environ.get("BWT_COMPETITOR_ROOT", "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"))
    p.add_argument("--wp0-root", default=os.environ.get("BWT_WP0_ROOT", "/data/gpfs/assoc/pgl/devel/exp1_human/wp0"))
    p.add_argument("--score-exp3", default=os.environ.get("BWT_SCORE_EXP3", "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"))
    a, gaps = p.parse_known_args(); bwt = need(a.bwt_bed, "primary BWTandem BED"); scorer = need(a.score_exp3, "score_exp3 module")
    mz = "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic"
    tools = [("BWTandem", bwt), ("TRF", os.path.join(a.competitor_root, f"trf/{mz}_trf_exp3B_satellite.bed")), ("ULTRA", os.path.join(a.wp0_root, "beds/ultra_maize_exp3B.bed")), ("tantan (500bp)", os.path.join(a.wp0_root, "fixcampaign/tantan/tantan_maize_w500.bed"))]
    tools = [(n, need(x, n)) for n, x in tools]
    for fn in ("mo17_knob180_arrays.bed", "mo17_tr1_arrays.bed", "mo17_centc_arrays.bed"): need(os.path.join(a.gt_root, fn), fn)
    source = os.path.join(HERE, "score_maize_postmerge.py")
    text = open(source).read().replace('SCORER = "/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/score_exp3.py"', f"SCORER = {scorer!r}", 1)
    mod = type(sys)("_frozen_maize_postmerge"); mod.__file__ = source; exec(compile(text, source, "exec"), mod.__dict__)
    mod.GT, mod.TOOLS = os.path.abspath(a.gt_root), tools; sys.argv = [sys.argv[0], *gaps]; mod.main()
if __name__ == "__main__": main()
