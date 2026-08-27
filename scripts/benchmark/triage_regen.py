#!/usr/bin/env python3
"""Triage a regenerated BED before any table is rescored.

Implements the agreed triage order: (1) schema and count sanity against the run
log/provenance, (2) per-tier call and base totals old vs new, (3) which changed
calls touch ambiguous (N) spans of the input, (4) satellite (tier 4) and
catch-all (tier 5) movement, (5) Tier 3 period-band redistribution
(100-200 / 100-500 / 101-2000). Interpretation and rescoring come only after
this passes.

Usage:
  triage_regen.py NEW_BED OLD_BED[.gz] FASTA [--provenance JSON]
"""
import argparse
import gzip
import json
import re
import sys
from collections import Counter, defaultdict


def load_bed(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as f:
        for i, line in enumerate(f, 1):
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                sys.exit(f"FATAL {path}:{i}: fewer than 3 columns")
            rows.append(p)
    return rows


def n_spans(fasta):
    spans = defaultdict(list)
    name, seq = None, []
    def flush():
        if name is None:
            return
        s = "".join(seq).upper()
        for m in re.finditer(r"[^ACGT]+", s):
            spans[name].append((m.start(), m.end()))
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    flush()
    return spans


def touches(spans, chrom, s, e):
    for a, b in spans.get(chrom, []):
        if s < b and e > a:
            return True
    return False


def summarize(rows, label):
    ncol = Counter(len(p) for p in rows)
    print(f"[{label}] {len(rows):,} calls; column widths: {dict(ncol)}")
    bad = 0
    tiers = Counter()
    bases = Counter()
    per_tier_calls = Counter()
    band = Counter()
    for p in rows:
        s, e = int(p[1]), int(p[2])
        if e <= s or s < 0:
            bad += 1
        tier = p[5] if len(p) > 5 else "?"
        tiers[tier] += 1
        bases[tier] += e - s
        per_tier_calls[tier] += 1
        if tier == "3":
            plen = len(p[3]) if len(p) > 3 else 0
            for lo, hi, name in ((100, 200, "100-200"), (100, 500, "100-500"),
                                 (101, 2000, "101-2000")):
                if lo <= plen <= hi:
                    band[name] += 1
    print(f"  invalid intervals: {bad}")
    for t in sorted(tiers):
        print(f"  tier {t}: {tiers[t]:>9,} calls  {bases[t]:>13,} bp")
    if band:
        print(f"  tier-3 period bands: {dict(band)}")
    return tiers, bases


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("new_bed")
    ap.add_argument("old_bed")
    ap.add_argument("fasta")
    ap.add_argument("--provenance")
    args = ap.parse_args()

    new = load_bed(args.new_bed)
    old = load_bed(args.old_bed)

    print("=== 1) schema & count sanity ===")
    nt, nb = summarize(new, "NEW")
    ot, ob = summarize(old, "OLD")
    if args.provenance:
        prov = json.load(open(args.provenance))
        print(f"  provenance: state={prov['state']} commit={prov['commit'][:9]} "
              f"elapsed={prov.get('elapsed_wall_s')}s rss_kb={prov.get('max_rss_kb')}")
        for o in prov.get("outputs", []):
            print(f"    {o['file']}  sha256={o['sha256'][:16]}...  {o['bytes']:,} B")

    print("\n=== 2) old vs new deltas per tier ===")
    for t in sorted(set(nt) | set(ot)):
        dc = nt.get(t, 0) - ot.get(t, 0)
        db = nb.get(t, 0) - ob.get(t, 0)
        print(f"  tier {t}: calls {ot.get(t,0):,} -> {nt.get(t,0):,} ({dc:+,}); "
              f"bases {ob.get(t,0):,} -> {nb.get(t,0):,} ({db:+,})")

    print("\n=== 3) changed calls touching N spans ===")
    key = lambda p: (p[0], p[1], p[2], p[3] if len(p) > 3 else "")
    old_set = {key(p) for p in old}
    new_set = {key(p) for p in new}
    removed = [p for p in old if key(p) not in new_set]
    added = [p for p in new if key(p) not in old_set]
    print(f"  removed: {len(removed):,}   added: {len(added):,}")
    spans = n_spans(args.fasta)
    total_n = sum(b - a for v in spans.values() for a, b in v)
    print(f"  input non-ACGT spans: {sum(len(v) for v in spans.values()):,} "
          f"({total_n:,} bp)")
    for label, group in (("removed", removed), ("added", added)):
        hit = sum(1 for p in group if touches(spans, p[0], int(p[1]), int(p[2])))
        print(f"  {label} touching N spans: {hit:,} / {len(group):,}")
    for p in removed[:5]:
        print(f"    - {p[0]}:{p[1]}-{p[2]} tier={p[5] if len(p)>5 else '?'} "
              f"motif_len={len(p[3]) if len(p)>3 else 0}")

    print("\n=== 4) satellite (tier 4) / catch-all (tier 5) movement ===")
    for t in ("4", "5"):
        r = sum(1 for p in removed if len(p) > 5 and p[5] == t)
        a = sum(1 for p in added if len(p) > 5 and p[5] == t)
        print(f"  tier {t}: removed {r:,}, added {a:,}")

    print("\ntriage complete -- interpret before rescoring.")


if __name__ == "__main__":
    main()
