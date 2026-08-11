#!/usr/bin/env python3
"""Characterize a tool's calls that fall OUTSIDE the adotto catalog.

Are the "false positives" (calls not overlapping adotto) genuine garbage, or
real repeats that other de novo tools (ultra/tantan) also call but adotto simply
did not catalog? The latter would mean the adotto-precision metric is unfair.

Usage: fp_check.py tool.bed adotto.bed ultra.bed tantan.bed
"""
import os
import subprocess
import sys

BEDTOOLS = os.environ.get("BEDTOOLS", "/data/gpfs/assoc/pgl/bin/bedtools2/bin/bedtools")


def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout


def count(bed):
    return sum(1 for l in open(bed) if l.strip())


def main():
    tool, adotto, ultra, tantan = sys.argv[1:5]
    total = count(tool)
    # calls NOT overlapping adotto
    fp = sh([BEDTOOLS, "intersect", "-a", tool, "-b", adotto, "-v"])
    with open("/tmp/_fp.bed", "w") as f:
        f.write(fp)
    n_fp = sum(1 for l in fp.split("\n") if l.strip())
    # of those, how many overlap ultra OR tantan
    supp = sh([BEDTOOLS, "intersect", "-a", "/tmp/_fp.bed", "-b", ultra, tantan, "-u"])
    n_supp = sum(1 for l in supp.split("\n") if l.strip())
    n_tp = total - n_fp
    print(f"tool={os.path.basename(tool)}")
    print(f"  total calls            : {total}")
    print(f"  overlap adotto (TP)    : {n_tp} ({100*n_tp/total:.1f}%)")
    print(f"  NOT in adotto (FP)     : {n_fp} ({100*n_fp/total:.1f}%)")
    if n_fp:
        print(f"    of FP, ultra/tantan-supported: {n_supp} ({100*n_supp/n_fp:.1f}% of FP)")
        print(f"    of FP, unsupported (garbage) : {n_fp-n_supp} ({100*(n_fp-n_supp)/n_fp:.1f}% of FP)")
        adj_prec = 100 * (n_tp + n_supp) / total
        print(f"  adjusted precision (adotto OR ultra/tantan): {adj_prec:.1f}%")


if __name__ == "__main__":
    main()
