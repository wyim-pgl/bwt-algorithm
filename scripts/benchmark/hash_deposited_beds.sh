#!/usr/bin/env bash
# Full-file SHA-256 for every deposited benchmark artifact (#13).
#
# results/manifest.tsv records first-megabyte hashes (sha256_1MB), which
# cannot detect truncation or tail corruption. This writes a companion file
# results/manifest.sha256 with whole-file hashes for everything under
# results/beds/ plus the manifest itself, verifiable with `sha256sum -c`.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
OUT=results/manifest.sha256
: > "$OUT"
find results/beds results/ground_truth -type f 2>/dev/null | sort | while read -r f; do
    sha256sum "$f" >> "$OUT"
done
sha256sum results/manifest.tsv >> "$OUT"
echo "$(wc -l < "$OUT") full-file hashes -> $OUT"
