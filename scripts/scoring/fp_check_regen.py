#!/usr/bin/env python3
"""Fail-closed, collision-safe variant of fp_check.py at commit 0363d8b."""
import argparse, os, subprocess, tempfile
BEDTOOLS = os.environ.get("BEDTOOLS", "/data/gpfs/assoc/pgl/bin/bedtools2/bin/bedtools")
def need(path, label):
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(path) or os.path.getsize(path) == 0: raise SystemExit(f"FATAL: {label} missing or empty: {path}")
    return path
def run(cmd, stdout=None):
    try: return subprocess.run(cmd, stdout=stdout, capture_output=stdout is None, text=True, check=True)
    except (OSError, subprocess.CalledProcessError) as e: raise SystemExit(f"FATAL: command failed: {' '.join(cmd)}\n{getattr(e, 'stderr', '') or e}")
def lines(path):
    with open(path) as f: return sum(bool(x.strip()) for x in f)
def main():
    p=argparse.ArgumentParser(); [p.add_argument(x) for x in ("tool","adotto","ultra","tantan")]; a=p.parse_args()
    tool,adotto,ultra,tantan=[need(x,n) for x,n in zip((a.tool,a.adotto,a.ultra,a.tantan),("primary tool BED","adotto BED","ULTRA BED","tantan BED"))]
    total=lines(tool)
    if not total: raise SystemExit(f"FATAL: primary tool BED has no records: {tool}")
    with tempfile.TemporaryDirectory(prefix="bwt-fp-check-") as td:
        fp=os.path.join(td,"fp.bed")
        with open(fp,"w") as out: run([BEDTOOLS,"intersect","-a",tool,"-b",adotto,"-v"],out)
        nfp=lines(fp); ns=0
        if nfp: ns=sum(bool(x.strip()) for x in run([BEDTOOLS,"intersect","-a",fp,"-b",ultra,tantan,"-u"]).stdout.splitlines())
    ntp=total-nfp
    print(f"tool={os.path.basename(tool)}\n  total calls            : {total}\n  overlap adotto (TP)    : {ntp} ({100*ntp/total:.1f}%)\n  NOT in adotto (FP)     : {nfp} ({100*nfp/total:.1f}%)")
    if nfp: print(f"    of FP, ultra/tantan-supported: {ns} ({100*ns/nfp:.1f}% of FP)\n    of FP, unsupported (garbage) : {nfp-ns} ({100*(nfp-ns)/nfp:.1f}% of FP)\n  adjusted precision (adotto OR ultra/tantan): {100*(ntp+ns)/total:.1f}%")
if __name__ == "__main__": main()
