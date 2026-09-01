#!/usr/bin/env bash
# Full-file SHA-256 for the deposited evidence (#13, #22).
#
# results/manifest.tsv records first-megabyte hashes (sha256_1MB), which
# cannot detect truncation or tail corruption. This writes:
#   results/manifest.sha256          — whole-file hashes for everything under
#                                      results/ (plus the evidence scripts),
#                                      verifiable on any clone with
#                                      `sha256sum -c results/manifest.sha256`
#   results/external_evidence.sha256 — whole-file hashes for every manifest
#                                      source_bed that lives OUTSIDE the repo
#                                      (benchmark workspace paths); only
#                                      verifiable on the cluster, but closes
#                                      the truncation blindness for the
#                                      largest evidence files
# then verifies the in-repo manifest it just wrote.
#
# ⚠️ Run this LAST, after every other edit under results/ — a results/ file
# edited after hashing makes `sha256sum -c` fail on a fresh clone (that
# exact failure shipped in e178707).
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
OUT=results/manifest.sha256
EXT=results/external_evidence.sha256

# Evidence scripts deposited beside the scoring wrappers; every listed file
# must exist — a silently skipped rename would shrink the verified set.
EVIDENCE_SCRIPTS=(
    scripts/scoring/score_maize_regen_evidence.py
    scripts/scoring/score_one_to_one.py
    scripts/scoring/derive_adotto_annotated_truth.py
)
for f in "${EVIDENCE_SCRIPTS[@]}"; do
    [[ -f "$f" ]] || { echo "FATAL: evidence script missing: $f" >&2; exit 1; }
done

: > "$OUT"
{ find results -type f ! -name manifest.sha256 ! -name external_evidence.sha256 \
      ! -path '*__pycache__*' ! -name '*.pyc'
  printf '%s\n' "${EVIDENCE_SCRIPTS[@]}"
} | sort | while read -r f; do
    sha256sum "$f" >> "$OUT"
done

: > "$EXT"
missing=0
while IFS= read -r f; do
    if [[ -f "$f" ]]; then
        sha256sum "$f" >> "$EXT"
    else
        echo "WARN: external evidence path unreachable, not hashed: $f" >&2
        missing=$((missing + 1))
    fi
done < <(awk -F'\t' 'NR>1 && $3 ~ /^\// {print $3}' results/manifest.tsv | sort -u)

sha256sum -c "$OUT" --quiet
echo "$(wc -l < "$OUT") in-repo hashes -> $OUT (verified)"
echo "$(wc -l < "$EXT") external evidence hashes -> $EXT (${missing} unreachable)"
