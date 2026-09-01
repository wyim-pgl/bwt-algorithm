"""Every STRfinder row must parse back to the same number of fields as the header.

`STR_genotype_structure` is specified as `motif_len[MOTIF]copies,truncated`, so it
legitimately contains a comma. It used to be interpolated into an f-string
unquoted, which gave every data row 13 fields against the 12-field header and
shifted `STR_genotype` and everything after it one column to the left. A consumer
reading `Variations` got `Full_seq`, and `Full_seq` got a truncated core sequence.
"""
import csv
import io
import os
import subprocess
import sys

import pytest

from bwtandem.models import TandemRepeat, STRFINDER_HEADER, STRFINDER_COLUMNS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(REPO_ROOT, "tests", "fixtures")


def parse_row(line: str):
    return next(csv.reader(io.StringIO(line)))


def make_repeat(**kw) -> TandemRepeat:
    fields = dict(chrom="chr1", start=100, end=118, motif="AT", copies=9.0,
                  length=18, tier=1, consensus_motif="AT", n_copies_evaluated=9,
                  percent_matches=100.0)
    fields.update(kw)
    return TandemRepeat(**fields)


def test_simple_repeat_row_has_header_field_count():
    row = parse_row(make_repeat().to_strfinder())
    assert len(row) == STRFINDER_COLUMNS


def test_genotype_structure_keeps_its_comma_in_one_field():
    """The spec wants the comma; quoting is what makes it survive the round trip."""
    # 18 bp span, 2 bp motif, 9 complete copies -> no truncated remainder.
    row = parse_row(make_repeat().to_strfinder())
    assert row[3] == "2[AT]9,0"
    assert row[4] == "9", "STR_genotype must not be the tail of the previous field"


def test_truncated_remainder_is_still_one_field():
    # 19 bp span over a 2 bp motif: 9 complete copies with 1 bp left over.
    row = parse_row(make_repeat(end=119, copies=9.5, length=19).to_strfinder())
    assert row[3] == "2[AT]9,1"
    assert len(row) == STRFINDER_COLUMNS


def test_variations_land_in_the_last_column():
    """The column that the off-by-one shift silently emptied."""
    r = make_repeat(variations=["1:2:A>G", "3:1:del(1)"])
    row = parse_row(r.to_strfinder())
    assert row[-1] == "1:2:A>G;3:1:del(1)"


def test_compound_repeat_row_has_header_field_count():
    partner = make_repeat(start=118, end=136, motif="GC", consensus_motif="GC")
    r = make_repeat(is_compound=True, compound_partner=partner)
    row = parse_row(r.to_strfinder())
    assert len(row) == STRFINDER_COLUMNS
    assert row[3] == "2[AT]9;2[GC]9,0"
    assert row[4] == "9/9"


def test_header_and_rows_agree_end_to_end(tmp_path):
    """The reported reproducer: synth_mixed.fa --format strfinder."""
    out_prefix = str(tmp_path / "mixed")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    proc = subprocess.run(
        [sys.executable, "-m", "bwtandem.main",
         os.path.join(FIXTURES, "synth_mixed.fa"), "--format", "strfinder",
         "-o", out_prefix],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=900,
    )
    assert proc.returncode == 0, proc.stderr

    with open(out_prefix + ".csv", newline="") as fh:
        rows = list(csv.reader(fh))

    assert rows[0] == STRFINDER_HEADER.split(",")
    assert len(rows) > 1, "fixture produced no records to check"
    bad = [(i, len(r)) for i, r in enumerate(rows[1:], 1) if len(r) != STRFINDER_COLUMNS]
    assert not bad, f"{len(bad)} rows do not have {STRFINDER_COLUMNS} fields: {bad[:5]}"
