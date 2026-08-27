#!/usr/bin/env bash
# Run one benchmark command with a complete provenance record (#7 / #13).
#
# Every manuscript result must be traceable as
#   clean commit -> exact command -> environment -> output checksum,
# which none of the historic runs were (results/README.md documents BEDs
# produced by a dirty tree). This wrapper refuses to run from a dirty tree
# unless --allow-dirty is given, and even then stamps the record DIRTY so the
# output can never silently masquerade as reproducible.
#
# Usage:
#   scripts/benchmark/run_with_provenance.sh [--allow-dirty] \
#       --tag NAME --out OUTPUT_PREFIX -- <command...>
#
# Writes OUTPUT_PREFIX.provenance.json next to the run outputs. The command is
# executed under /usr/bin/time -v, so elapsed wall time and MaxRSS come from
# the kernel, not self-reporting.
set -euo pipefail

ALLOW_DIRTY=0
TAG=""
OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --allow-dirty) ALLOW_DIRTY=1; shift ;;
        --tag) TAG="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --) shift; break ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$TAG" && -n "$OUT" && $# -gt 0 ]] || {
    echo "usage: $0 [--allow-dirty] --tag NAME --out PREFIX -- <command...>" >&2
    exit 2
}

REPO_ROOT=$(git rev-parse --show-toplevel)
COMMIT=$(git -C "$REPO_ROOT" rev-parse HEAD)
DIRTY=$(git -C "$REPO_ROOT" status --porcelain --untracked-files=no | head -c1 || true)
if [[ -n "$DIRTY" ]]; then
    if [[ $ALLOW_DIRTY -eq 1 ]]; then
        echo "WARNING: working tree is dirty; record will be stamped DIRTY." >&2
        STATE="DIRTY"
    else
        echo "FATAL: working tree has uncommitted changes. Commit them or pass" >&2
        echo "--allow-dirty (the run will then be unusable for the manuscript)." >&2
        exit 1
    fi
else
    STATE="CLEAN"
fi

PYBIN=$(command -v python3)
TIMELOG=$(mktemp)
START_ISO=$(date -u +%Y-%m-%dT%H:%M:%SZ)

# GNU time is absent on some cluster nodes, so measure with wait4() rusage
# from a tiny python driver instead: same kernel numbers, no dependency.
set +e
"$PYBIN" - "$TIMELOG" "$@" <<'PYDRV'
import json, os, sys, time
logf, argv = sys.argv[1], sys.argv[2:]
t0 = time.monotonic()
pid = os.fork()
if pid == 0:
    os.execvp(argv[0], argv)
_, status, ru = os.wait4(pid, 0)
elapsed = time.monotonic() - t0
rc = os.waitstatus_to_exitcode(status)
with open(logf, "w") as fh:
    json.dump({"elapsed_s": round(elapsed, 2),
               "max_rss_kb": ru.ru_maxrss}, fh)
sys.exit(rc if rc >= 0 else 128 - rc)
PYDRV
RC=$?
set -e

ELAPSED=$("$PYBIN" -c "import json;print(json.load(open('$TIMELOG'))['elapsed_s'])" 2>/dev/null || echo "")
MAXRSS_KB=$("$PYBIN" -c "import json;print(json.load(open('$TIMELOG'))['max_rss_kb'])" 2>/dev/null || echo "")

# checksum every file the run produced under the output prefix
CHECKSUMS=$(for f in "$OUT".*; do
    [[ -f "$f" && "$f" != *.provenance.json ]] || continue
    printf '    {"file": "%s", "sha256": "%s", "bytes": %s},\n' \
        "$(basename "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")"
done | sed '$ s/,$//')

cat > "$OUT.provenance.json" <<JSON
{
  "tag": "$TAG",
  "state": "$STATE",
  "exit_code": $RC,
  "commit": "$COMMIT",
  "branch": "$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD)",
  "command": $(printf '%s\n' "$*" | "$PYBIN" -c 'import json,sys; print(json.dumps(sys.stdin.read().strip()))'),
  "started_utc": "$START_ISO",
  "finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "elapsed_wall_s": ${ELAPSED:-null},
  "max_rss_kb": ${MAXRSS_KB:-null},
  "host": "$(hostname)",
  "cpu": "$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')",
  "python": "$("$PYBIN" -c 'import sys; print(sys.version.split()[0])')",
  "numpy": "$("$PYBIN" -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo n/a)",
  "outputs": [
$CHECKSUMS
  ]
}
JSON
rm -f "$TIMELOG"
echo "provenance -> $OUT.provenance.json (state=$STATE, exit=$RC)"
exit $RC
