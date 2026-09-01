#!/usr/bin/env python3
"""Derive a motif/copy-annotated truth BED from the adotto TR catalog (#15).

The region-level scoring GT (`adotto_primary.bed`) carries coordinates only, so
one-to-one scoring against it cannot report motif-period or copy-number error.
The full catalog (`adotto_TRregions_v1.2.1.bed`) carries a JSON list of
per-region repeat annotations in its last column; this script emits a 5-column
BED — region chrom/start/end plus the motif and copy count of the
**highest-score annotation** in the region — for use as the annotated truth.

Choosing the top-score annotation is a disclosed simplification: compound
regions carry several motifs, and the strict scorer then measures period/copy
agreement against the dominant one. Regions whose annotation list is empty or
unparseable are emitted with no motif (they still score for overlap/boundary).

Usage: derive_adotto_annotated_truth.py adotto_TRregions_v1.2.1.bed > truth5.bed
"""
import json
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    n_rows = n_annotated = 0
    with open(sys.argv[1]) as f:
        for i, line in enumerate(f, 1):
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom, start, end = p[0], int(p[1]), int(p[2])
            motif, copies = "", ""
            try:
                anns = json.loads(p[-1])
                best = max((a for a in anns if a.get("motif")),
                           key=lambda a: a.get("score", 0), default=None)
                if best:
                    motif = best["motif"]
                    copies = best.get("copies", "")
            except (json.JSONDecodeError, TypeError, ValueError):
                pass
            n_rows += 1
            if motif:
                n_annotated += 1
            print(f"{chrom}\t{start}\t{end}\t{motif}\t{copies}")
    print(f"rows={n_rows} annotated={n_annotated}", file=sys.stderr)


if __name__ == "__main__":
    main()
