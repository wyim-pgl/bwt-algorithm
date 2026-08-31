#!/bin/bash
# BWTandem smoke test on the bundled test data (~30 s total).
# Run from the repository root: bash examples/quickstart.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

echo "== 1/3 toy STR (35 bp, expect exactly 1 call: TCATCGG x5)"
python3 -m src.main arabadopsis_chrs/test_seq1.fa --format bed -o "$OUT/toy"
cat "$OUT/toy.bed"
grep -q $'^repeatTCATCGG_5\t0\t35\tTCATCGG' "$OUT/toy.bed"

echo "== 2/3 real 367 kb sequence (Arabidopsis ChrM, expect ~45 calls)"
python3 -m src.main arabadopsis_chrs/ChrM.fa --format bed -o "$OUT/chrm"
n=$(wc -l < "$OUT/chrm.bed")
echo "   $n calls"
[ "$n" -ge 30 ] && [ "$n" -le 60 ]

echo "== 3/3 output formats"
python3 -m src.main arabadopsis_chrs/test_seq1.fa --format vcf -o "$OUT/toy" >/dev/null
python3 -m src.main arabadopsis_chrs/test_seq1.fa --format trf -o "$OUT/toy" >/dev/null
ls "$OUT"/toy.*

echo "SMOKE TEST PASSED"
