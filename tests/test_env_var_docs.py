"""Every tuning knob the code reads must be documented, and vice versa.

Why this exists (2026-08-05): `CATCHALL_MAX_SEEDS` was documented in CLAUDE.md with
a stated default of 200000 and did not exist anywhere in the source. Anyone setting
it got silence — no error, no effect, and a plausible-looking knob in the docs to
explain a result it had no part in. In the other direction `TIER1_FMSCAN` was absent
from the tuning documentation while being set to 1 by *every* benchmark run in the
project, including the one behind the published figures; the manuscript's "Tier 1
used FM-index enumeration mode" is that knob, and a reader had no way to find it.

Both directions are failures of the same kind: the documented interface and the real
interface drifted apart, silently. This test pins them together.

If this fails, do not edit the test to make it pass. Either document the new knob in
CLAUDE.md's "Tuning detection sensitivity (env vars)" section with its default, or
delete the documentation for the knob that no longer exists.
"""
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "src")
CLAUDE_MD = os.path.join(REPO, "CLAUDE.md")

# Prefixes that denote a user-facing tuning knob. Anything outside these (PATH,
# HOME, SLURM_*, ...) is not part of the tool's interface and is not checked.
PREFIXES = ("TIER1_", "TIER2_", "TIER3_", "CATCHALL_", "SAT_FILL_", "SAT_ANCHOR_", "BWT_")

_READ = re.compile(r"""environ(?:\.get\(|\[)\s*["']([A-Z0-9_]+)["']""")
_DOC = re.compile(r"\b((?:TIER[123]|CATCHALL|SAT_FILL|SAT_ANCHOR|BWT)_[A-Z0-9_]+)\b")


def _knobs_read_by_code():
    """Env var names read anywhere under src/, restricted to the tuning prefixes."""
    found = set()
    for root, _dirs, files in os.walk(SRC):
        for fn in files:
            if not fn.endswith((".py", ".pyx")):
                continue
            with open(os.path.join(root, fn), encoding="utf-8", errors="replace") as fh:
                found.update(_READ.findall(fh.read()))
    return {n for n in found if n.startswith(PREFIXES)}


def _knobs_documented():
    with open(CLAUDE_MD, encoding="utf-8") as fh:
        return set(_DOC.findall(fh.read()))


def test_no_documented_knob_is_absent_from_the_code():
    """A documented knob that nothing reads is worse than an undocumented one.

    It looks settable, so a user sets it, and the run silently ignores it.
    """
    phantom = sorted(_knobs_documented() - _knobs_read_by_code())
    assert not phantom, (
        "CLAUDE.md documents environment variables that no source file reads, so "
        "setting them does nothing: " + ", ".join(phantom) + ". Remove them from the "
        "docs, or implement them."
    )


def test_no_knob_the_code_reads_is_undocumented():
    """A knob nobody can find is a knob whose effect on a result cannot be reported."""
    undocumented = sorted(_knobs_read_by_code() - _knobs_documented())
    assert not undocumented, (
        "src/ reads environment variables that CLAUDE.md does not document: "
        + ", ".join(undocumented)
        + ". Add them to the 'Tuning detection sensitivity (env vars)' section with "
        "their defaults."
    )


def test_the_scan_finds_a_plausible_number_of_knobs():
    """Guard against the regexes silently matching nothing and passing vacuously."""
    read = _knobs_read_by_code()
    assert len(read) >= 30, (
        f"only {len(read)} tuning knobs found in src/ — the scan regex has probably "
        "stopped matching, which would make the two tests above pass for the wrong "
        "reason"
    )
    assert "TIER1_FMSCAN" in read, "TIER1_FMSCAN should be read by src/tier1.py"
    assert "BWT_INDEX_CACHE" in read, "BWT_INDEX_CACHE should be read by src/finder.py"
