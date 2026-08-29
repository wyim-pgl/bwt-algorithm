#!/usr/bin/env bash
# External acceptance gate for one regenerated genome (UNTRACKED during the
# campaign: the pin check forbids commits and the dirty check forbids edits to
# tracked files until every job has started; a new untracked file passes both).
#
# Accepts the run only when all of the following agree (the four-way agreement
# rule): SLURM COMPLETED with ExitCode 0; provenance exit_code 0, state CLEAN,
# commit == EXPECT_COMMIT; exactly one "Total repeats found: N" in the job log
# and wc -l BED == N; recomputed SHA-256 and byte count equal the provenance
# record; and every BED row passes the 8-column schema check.
#
# Usage: acceptance_gate.sh GENOME JOBID   (e.g. acceptance_gate.sh colcen 6110900)
set -uo pipefail
G=${1:?genome}; J=${2:?jobid}
R=/data/gpfs/assoc/pgl/devel/exp1_human/regen
BED=$R/regen_$G.bed
PROV=$R/regen_$G.provenance.json
LOG=$(ls $R/logs/regen_regen_${G}_${J}.out 2>/dev/null || ls $R/logs/*_${J}.out 2>/dev/null | head -1)
PY=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python
export PATH=/cm/shared/apps/slurm/current/bin:$PATH
fail=0
# NB: must not end on a failing test -- `[ PASS = FAIL ] && fail=1` returns 1,
# which made every `note PASS ... || note FAIL ...` call site take BOTH branches.
note() { echo "  [$1] $2"; if [ "$1" = FAIL ]; then fail=1; fi; return 0; }

echo "=== acceptance gate: $G (job $J) ==="

st=$(sacct -j "$J" -X --format=State,ExitCode -Pn)
[[ "$st" == "COMPLETED|0:0" ]] && note PASS "SLURM $st" || note FAIL "SLURM state '$st' != COMPLETED|0:0"

if [[ -f "$PROV" ]]; then
  $PY -m json.tool "$PROV" >/dev/null 2>&1 && note PASS "provenance is valid JSON" || note FAIL "provenance JSON invalid"
  pe=$($PY -c "import json;p=json.load(open('$PROV'));print(p['exit_code'],p['state'],p['commit'])")
  read -r pex pst pcm <<< "$pe"
  [[ "$pex" == 0 ]] && note PASS "provenance exit_code 0" || note FAIL "provenance exit_code $pex"
  [[ "$pst" == CLEAN ]] && note PASS "provenance state CLEAN" || note FAIL "provenance state $pst"
  exp=$(cat $R/EXPECT_COMMIT)
  [[ "$pcm" == "$exp" ]] && note PASS "commit matches pin ${exp:0:9}" || note FAIL "commit $pcm != pin $exp"
else
  note FAIL "provenance missing: $PROV"
fi

if [[ -n "${LOG:-}" && -f "$LOG" ]]; then
  nl=$(grep -c "^Total repeats found:" "$LOG")
  [[ "$nl" == 1 ]] && note PASS "exactly one total-count line" || note FAIL "$nl total-count lines"
  total=$(grep "^Total repeats found:" "$LOG" | awk '{print $4}')
  wl=$(wc -l < "$BED" 2>/dev/null || echo MISSING)
  [[ "$total" == "$wl" ]] && note PASS "log total $total == wc -l $wl" || note FAIL "log total $total != wc -l $wl"
else
  note FAIL "job log not found under $R/logs for job $J"
fi

[[ -s "$BED" ]] && note PASS "BED is non-empty" || note FAIL "BED missing or empty: $BED"

if [[ -f "$BED" && -f "$PROV" ]]; then
  want=$($PY -c "
import json, sys
beds = [o for o in json.load(open('$PROV'))['outputs'] if o['file'].endswith('.bed')]
if len(beds) != 1:
    sys.exit(f'expected exactly one .bed output in provenance, found {len(beds)}')
print(beds[0]['sha256'], beds[0]['bytes'])") || { note FAIL "provenance outputs: $want"; want=""; }
  read -r wsha wbytes <<< "$want"
  gsha=$(sha256sum "$BED" | cut -d' ' -f1); gbytes=$(stat -c%s "$BED")
  [[ "$gsha" == "$wsha" && "$gbytes" == "$wbytes" ]] \
    && note PASS "sha256/bytes match provenance" \
    || note FAIL "sha256/bytes drifted: got $gsha/$gbytes want $wsha/$wbytes"
  $PY - "$BED" <<'PYEOF'
import sys
bad = rows = 0
tiers = set()
with open(sys.argv[1]) as f:
    for i, line in enumerate(f, 1):
        p = line.rstrip("\n").split("\t")
        rows += 1
        try:
            ok = (len(p) == 8 and int(p[2]) > int(p[1]) >= 0
                  and p[5] in {"1","2","3","4","5"} and p[7] in {"+","-"}
                  and float(p[4]) > 0 and 0.0 <= float(p[6]) <= 1.0 and p[3])
        except (ValueError, IndexError):
            ok = False
        if not ok:
            bad += 1
            if bad <= 3: print(f"  bad row {i}: {line.rstrip()[:100]}", file=sys.stderr)
        else:
            tiers.add(p[5])
print(f"  schema: {rows:,} rows, {bad} bad, tiers seen {sorted(tiers)}")
sys.exit(1 if (bad or rows == 0) else 0)
PYEOF
  [[ $? == 0 ]] && note PASS "8-column schema clean" || note FAIL "schema violations"
else
  note FAIL "BED missing: $BED"
fi

echo
if [[ $fail == 0 ]]; then echo ">>> ACCEPTED: $G may proceed to triage/scoring."; else echo ">>> REJECTED: fix before any scoring."; exit 1; fi
