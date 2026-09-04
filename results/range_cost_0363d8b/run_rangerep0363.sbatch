#!/bin/bash
#SBATCH --job-name=rangerep0363
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=44G
#SBATCH --time=24:00:00
#SBATCH --output=/data/gpfs/assoc/pgl/devel/exp1_human/wp0/rangerep0363/rangerep0363_%j.log
set -uo pipefail
: "${REP:?set REP=1|2|3}"
W=/data/gpfs/assoc/pgl/devel/exp1_human/wp0/rangerep0363
PY=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python
FASTA=/data/gpfs/assoc/pgl/devel/exp1_human/data/hg38_primary.fa
export PATH=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin:$PATH   # build.py needs gcc
cd "$W/tree"

# C-1(b): the published range-cost pairs all ran on 07ad6fa, where Tier 3's window was
# fixed at 100 bp-100 kb regardless of --max-period. At 0363d8b the ceiling bounds Tier 3
# as well, so the narrow arm no longer performs the wide search and discards it. Both arms
# run in ONE job so the pair is guaranteed to share a node, which is what makes the ratio
# interpretable. Identical to run_rangerep.sbatch except for the code commit.
echo "START $(date) HOST $(hostname) REP=$REP CODE=0363d8b(detached worktree)"
git -C "$W/tree" rev-parse HEAD
ls -l --time-style=full-iso src/*.py | md5sum

export TIER1_FMSCAN=1 TIER1_FMSCAN_MIN_DENSITY=0.45 TIER1_FMSCAN_MIN_LLR=6.0 \
       TIER1_MIN_ARRAY_LEN=20 TIER1_MIN_SCORE=20 TIER1_MIN_COPIES=2 \
       TIER1_COPYBASE=6 TIER1_COPYADD=2 TIER1_EXT_COPIES=2 \
       TIER1_SHORT_PERIOD_MAX=9 TIER1_SHORT_MIN_ARRAY_LEN=17 TIER1_SHORT_MIN_SCORE=17 \
       TIER2_MISMATCH=0.30 TIER2_SHORT_REQ_COPIES=2 \
       CATCHALL_SCAN=1 CATCHALL_MIN_IDENTITY=0.72 CATCHALL_MIN_COPIES=3

for MAXP in 100 2000; do
  OUT="$W/rep${REP}_p${MAXP}"
  echo "--- ARM MAXP=$MAXP start $(date) ---"
  /usr/bin/time -v "$PY" -m src.main "$FASTA" --min-period 1 --max-period "$MAXP" \
    --threads 4 --format bed -o "$OUT" -v > "$OUT.stdout" 2> "$OUT.time"
  echo "exit=$? calls=$(wc -l < "$OUT.bed" 2>/dev/null)"
  grep -E "Elapsed|Maximum resident" "$OUT.time"
done
echo "DONE $(date)"
