#!/usr/bin/env python3
"""Derive a motif/copy-annotated truth BED from the adotto TR catalog (#15).

The region-level scoring GT (`adotto_primary.bed`) carries coordinates only, so
one-to-one scoring against it cannot report motif-period or copy-number error.
The full catalog (`adotto_TRregions_v1.2.1.bed`) carries a JSON list of
per-region repeat annotations in its last column; this script emits a 5-column
BED for use as the annotated truth. Design decisions (each surfaced by the
2026-09-01 review of the first deposited measurement):

- **One annotation per region, the highest-score one** — a disclosed
  simplification for compound regions: a tool that correctly reports one of
  the region's other motifs scores as a period disagreement.
- **The ANNOTATION's own coordinates**, not the region's. Region bounds are
  merged and slop-padded (~25 bp beyond the annotation at both ends), so
  scoring boundary error against them imposes a ~25 bp floor on every tool,
  and the annotation's copy count covers only the annotation sub-interval.
- **Primitive-period reduction, with the copy count rescaled.** ~1% of
  catalog motifs are non-primitive (an 87-mer that is a 29-mer ×3); tools
  that reduce to the primitive period by design (BWTandem does) would be
  demoted from 'exact' and charged ~200% copy error on a perfect call.
  The emitted motif is the primitive unit and copies are multiplied by the
  reduction factor.

Regions whose annotation list is empty or of an unexpected shape are emitted
with the REGION coordinates and no motif (they still score for overlap);
malformed catalog lines abort loudly rather than truncating the truth.

Usage: derive_adotto_annotated_truth.py adotto_TRregions_v1.2.1.bed \
           [--restrict-to REGIONS.bed] > truth5.bed

--restrict-to keeps only catalog rows whose REGION coordinates appear in the
given 3+ column BED (the primary scoring GT), so the annotated truth covers
exactly the scored region set.
"""
import argparse
import json
import sys


def primitive(motif):
    """Smallest repeating unit of `motif` and the reduction factor."""
    n = len(motif)
    for p in range(1, n // 2 + 1):
        if n % p == 0 and motif == motif[:p] * (n // p):
            return motif[:p], n // p
    return motif, 1


def best_annotation(cell):
    try:
        anns = json.loads(cell)
    except json.JSONDecodeError:
        return None
    if not isinstance(anns, list):
        return None
    best = None
    for a in anns:
        if not isinstance(a, dict) or not a.get("motif"):
            continue
        if best is None or a.get("score", 0) > best.get("score", 0):
            best = a
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("catalog_bed")
    ap.add_argument("--restrict-to", metavar="REGIONS_BED",
                    help="keep only catalog rows whose region coordinates "
                         "appear in this BED")
    args = ap.parse_args()

    keep = None
    if args.restrict_to:
        keep = set()
        with open(args.restrict_to) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3 and not line.startswith(("#", "track")):
                    keep.add((p[0], p[1], p[2]))

    n_rows = n_annotated = n_reduced = 0
    with open(args.catalog_bed) as f:
        for i, line in enumerate(f, 1):
            if line.startswith(("#", "track")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom = p[0]
            try:
                r_start, r_end = int(p[1]), int(p[2])
            except ValueError:
                sys.exit(f"FATAL {args.catalog_bed}:{i}: bad region coordinates")
            if keep is not None and (chrom, p[1], p[2]) not in keep:
                continue
            best = best_annotation(p[-1])
            if best is None:
                print(f"{chrom}\t{r_start}\t{r_end}\t\t")
                n_rows += 1
                continue
            try:
                a_start, a_end = int(best["start"]), int(best["end"])
                copies = float(best.get("copies", 0)) or ""
            except (KeyError, TypeError, ValueError):
                sys.exit(f"FATAL {args.catalog_bed}:{i}: malformed annotation")
            motif, factor = primitive(str(best["motif"]).upper())
            if factor > 1:
                n_reduced += 1
                if copies != "":
                    copies = round(copies * factor, 1)
            print(f"{chrom}\t{a_start}\t{a_end}\t{motif}\t{copies}")
            n_rows += 1
            n_annotated += 1
    print(f"rows={n_rows} annotated={n_annotated} primitive_reduced={n_reduced}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
