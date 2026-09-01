#!/bin/bash
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=36:00:00
#SBATCH --output=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/regen_%x_%j.out
#SBATCH --error=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/regen_%x_%j.err

# Benchmark regeneration for issues #3/#4/#5 (GRCh38 / Mo17 / Col-CEN).
#
# Reproduces the recorded remeasurement configuration exactly (v2.2 gate,
# threads 2, period 1-2000; catch-all id 0.72 MIN_COPIES 3 on human and maize,
# catch-all OFF on colcen -- see wp0/run_remeasure.sbatch), but from a CLEAN,
# PINNED commit through scripts/benchmark/run_with_provenance.sh, which is the
# provenance discipline the deposited BEDs lacked (they trace to 0e17d1a plus
# uncommitted code). The BEDs produced here are EXPECTED to differ from the
# deposited ones: the detector gained ambiguous-base masking, an off-by-one fix
# in windowed counts, satellite segmentation, and Tier 3 period wiring since.
#
# Submit (after committing so HEAD == EXPECT_COMMIT):
#   sbatch --job-name=regen_colcen --export=ALL,GENOME=colcen scripts/benchmark/regen_345.sbatch
#   sbatch --job-name=regen_human  --export=ALL,GENOME=human  scripts/benchmark/regen_345.sbatch
#   sbatch --job-name=regen_maize  --export=ALL,GENOME=maize  scripts/benchmark/regen_345.sbatch
#
# After completion, record sacct MaxRSS (authoritative for multi-process runs;
# rusage max-over-children is not the sum):
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State

set -euo pipefail
: "${GENOME:?set GENOME=human|maize|colcen}"

REPO=/data/gpfs/assoc/pgl/devel/bwt-algorithm
PY=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python
OUTDIR=/data/gpfs/assoc/pgl/devel/exp1_human/regen
EXPECT_COMMIT_FILE=$OUTDIR/EXPECT_COMMIT
mkdir -p "$OUTDIR/logs"
cd "$REPO"

# The commit every run of this campaign must be executed from. Written once at
# submission time; a job landing on a moved HEAD refuses to run.
if [[ ! -f "$EXPECT_COMMIT_FILE" ]]; then
    echo "FATAL: $EXPECT_COMMIT_FILE not found. Write the campaign commit hash" >&2
    echo "into it before submitting (git rev-parse HEAD > $EXPECT_COMMIT_FILE)." >&2
    exit 1
fi
EXPECT=$(cat "$EXPECT_COMMIT_FILE")
HEAD_NOW=$(git rev-parse HEAD)
if [[ "$HEAD_NOW" != "$EXPECT" ]]; then
    echo "FATAL: HEAD is $HEAD_NOW but the campaign is pinned to $EXPECT." >&2
    exit 1
fi

# v2.2 gate base, common to all three genomes (from wp0/run_remeasure.sbatch).
export TIER1_FMSCAN=1 TIER1_FMSCAN_MIN_DENSITY=0.45 TIER1_FMSCAN_MIN_LLR=6.0 \
       TIER1_MIN_ARRAY_LEN=20 TIER1_MIN_SCORE=20 TIER1_MIN_COPIES=2 \
       TIER1_COPYBASE=6 TIER1_COPYADD=2 TIER1_EXT_COPIES=2 \
       TIER1_SHORT_PERIOD_MAX=9 TIER1_SHORT_MIN_ARRAY_LEN=17 TIER1_SHORT_MIN_SCORE=17 \
       TIER2_MISMATCH=0.30 TIER2_SHORT_REQ_COPIES=2

case "$GENOME" in
  human)
    FASTA=/data/gpfs/assoc/pgl/devel/exp1_human/data/hg38_primary.fa
    export CATCHALL_SCAN=1 CATCHALL_MIN_IDENTITY=0.72 CATCHALL_MIN_COPIES=3
    ;;
  maize)
    FASTA=/data/gpfs/assoc/pgl/filip/bwt/zmays/zmays.fna
    export CATCHALL_SCAN=1 CATCHALL_MIN_IDENTITY=0.72 CATCHALL_MIN_COPIES=3
    ;;
  colcen)
    FASTA=/data/gpfs/assoc/pgl/filip/bwt/colcen/colcen.fa
    # catch-all off for Arabidopsis, per Methods 2.2.1
    ;;
  *) echo "FATAL: unknown GENOME '$GENOME'"; exit 1 ;;
esac
OUT=$OUTDIR/regen_${GENOME}

echo "=== campaign provenance (job header) ==="
echo "jobid        : ${SLURM_JOB_ID:-none}   host: $(hostname)"
echo "genome       : $GENOME"
echo "fasta        : $FASTA"
echo "fasta sha256 : $(sha256sum "$FASTA" | cut -d' ' -f1)"
echo "fasta bytes  : $(stat -c %s "$FASTA")   sequences: $(grep -c '^>' "$FASTA")"
echo "commit       : $HEAD_NOW (pinned)"
echo "--- env knobs ---"
env | grep -E '^(TIER1_|TIER2_|CATCHALL_|SAT_|BWT_)' | sort
echo "========================================"

scripts/benchmark/run_with_provenance.sh --tag "regen_${GENOME}_${SLURM_JOB_ID:-manual}" --out "$OUT" -- \
    "$PY" -m bwtandem.main "$FASTA" \
    --min-period 1 --max-period 2000 --threads 2 --format bed -o "$OUT" -v
RC=$?

echo "BED lines: $(wc -l < "${OUT}.bed" 2>/dev/null || echo MISSING)"
exit $RC
