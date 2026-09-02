"""The deposited-evidence hash manifest must verify against the working tree.

Guard for a failure mode that shipped twice (e178707, 91f9763): a results/
file edited after `hash_deposited_beds.sh` ran makes `sha256sum -c` fail on
every fresh clone. CI checks this on push; this test catches it at the local
`pytest` run, before the commit exists. Fix is always the same: rerun
`bash scripts/benchmark/hash_deposited_beds.sh` LAST and commit the refreshed
manifest.
"""
import hashlib
import os

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MANIFEST = os.path.join(REPO, "results", "manifest.sha256")


def test_manifest_hashes_verify():
    assert os.path.isfile(MANIFEST), "results/manifest.sha256 missing"
    bad, missing = [], []
    with open(MANIFEST) as fh:
        for line in fh:
            digest, _, rel = line.strip().partition("  ")
            path = os.path.join(REPO, rel)
            if not os.path.isfile(path):
                missing.append(rel)
                continue
            h = hashlib.sha256()
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1 << 20), b""):
                    h.update(chunk)
            if h.hexdigest() != digest:
                bad.append(rel)
    if bad or missing:
        pytest.fail(
            "deposited-evidence hashes stale — rerun "
            "scripts/benchmark/hash_deposited_beds.sh LAST and commit it. "
            f"mismatched: {bad[:5]} missing: {missing[:5]}")


def test_manifest_paths_are_tracked_by_git():
    """Every hashed path must be committed, not merely present in the working
    tree. A `.gitignore` rule (`*.settings`, 2026-09-02) silently dropped an
    evidence file from `git add`; the hash test above still passed locally
    because the file existed on disk, and CI's `sha256sum -c` failed on the
    fresh clone. Skipped outside a git checkout (release tarballs)."""
    import subprocess
    try:
        out = subprocess.run(["git", "-C", REPO, "ls-files", "-z"],
                             capture_output=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout")
    tracked = set(out.decode().split("\0"))
    with open(MANIFEST) as fh:
        hashed = [line.strip().partition("  ")[2] for line in fh if line.strip()]
    untracked = [p for p in hashed if p not in tracked]
    assert not untracked, (
        "hashed evidence not tracked by git (ignored by .gitignore? use "
        f"`git add -f`): {untracked[:5]}")
