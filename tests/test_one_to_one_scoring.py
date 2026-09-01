"""Regression tests for scripts/scoring/score_one_to_one.py (issue #15).

The scorer is the strict counterpart to the permissive matcher in
test_ground_truth.py: one-to-one maximum-cardinality assignment, reciprocal
overlap, boundary/period/copy-number error. The adversarial review of
2026-08-27 required a counterexample regression test showing the greedy
descending-overlap heuristic is not maximum-cardinality; that test is here.
"""
import importlib.util
import json
import os
import subprocess
import sys

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "..",
                      "scripts", "scoring", "score_one_to_one.py")


def _load_module():
    spec = importlib.util.spec_from_file_location("score_one_to_one", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m = _load_module()


def rec(chrom, start, end, motif=None, copies=None):
    return {"chrom": chrom, "start": start, "end": end,
            "motif": motif, "copies": copies, "line": 0}


class TestMaximumCardinality:
    def test_greedy_counterexample(self):
        """Greedy descending overlap finds 1 match here; max-cardinality 2.

        T2-P1 is the single largest edge (120 bp). Greedy takes it first,
        which consumes T1's only partner; the optimal assignment T1-P1,
        T2-P2 matches both truths.
        """
        truth = [rec("c", 0, 100), rec("c", 80, 200)]
        preds = [rec("c", 0, 200), rec("c", 100, 200)]
        pairs, t_used, p_used = m.one_to_one(truth, preds, 0.3)
        assert len(pairs) == 2
        assert t_used == {0, 1} and p_used == {0, 1}

    def test_simple_exact(self):
        truth = [rec("c", 0, 100)]
        preds = [rec("c", 0, 100)]
        pairs, _, _ = m.one_to_one(truth, preds, 0.5)
        assert len(pairs) == 1
        assert pairs[0][2] == 100  # overlap bp

    def test_each_record_used_once(self):
        # one long prediction spanning two truths can satisfy only one
        truth = [rec("c", 0, 100), rec("c", 100, 200)]
        preds = [rec("c", 0, 200)]
        pairs, _, _ = m.one_to_one(truth, preds, 0.5)
        assert len(pairs) == 1

    def test_chromosomes_do_not_cross(self):
        truth = [rec("c1", 0, 100)]
        preds = [rec("c2", 0, 100)]
        pairs, _, _ = m.one_to_one(truth, preds, 0.5)
        assert len(pairs) == 0


class TestReciprocalOverlap:
    def test_tiny_prediction_inside_huge_truth_rejected(self):
        # covers the truth's bases poorly: 10/1000 < 0.5 of the truth
        truth = [rec("c", 0, 1000)]
        preds = [rec("c", 0, 10)]
        edges = m._eligible_edges(truth, preds, 0.5)
        assert not edges

    def test_huge_prediction_over_tiny_truth_rejected(self):
        # reciprocal: overlap must also cover 0.5 of the PREDICTION
        truth = [rec("c", 0, 10)]
        preds = [rec("c", 0, 1000)]
        edges = m._eligible_edges(truth, preds, 0.5)
        assert not edges

    def test_reciprocal_pass(self):
        truth = [rec("c", 0, 100)]
        preds = [rec("c", 25, 125)]  # overlap 75 = 0.75 of both
        edges = m._eligible_edges(truth, preds, 0.5)
        assert edges[0][0] == 75

    def test_sweepline_matches_bruteforce(self):
        import random
        rnd = random.Random(7)
        truth, preds = [], []
        for _ in range(120):
            s = rnd.randrange(0, 5000)
            truth.append(rec("c", s, s + rnd.randrange(10, 300)))
        for _ in range(150):
            s = rnd.randrange(0, 5000)
            preds.append(rec("c", s, s + rnd.randrange(10, 300)))
        edges = m._eligible_edges(truth, preds, 0.5)
        brute = {}
        for i, t in enumerate(truth):
            for j, p in enumerate(preds):
                ov = m.overlap(t, p)
                if (ov > 0 and ov >= 0.5 * (t["end"] - t["start"])
                        and ov >= 0.5 * (p["end"] - p["start"])):
                    brute.setdefault(i, {})[j] = ov
        assert {i: dict(v) for i, v in edges.items()} == brute


def _write_bed(tmp_path, name, rows):
    p = tmp_path / name
    p.write_text("".join("\t".join(str(c) for c in r) + "\n" for r in rows))
    return str(p)


def _run(args, **kw):
    return subprocess.run([sys.executable, SCRIPT, *args],
                          capture_output=True, text=True, **kw)


class TestCLI:
    def test_json_metrics(self, tmp_path):
        truth = _write_bed(tmp_path, "t.bed", [
            ("c", 0, 100, "AT", 50),        # exact period + copies match
            ("c", 200, 300, "ACG", 33),     # pred period 6 = integer multiple
            ("c", 400, 500, "ACGTA", 20),   # pred period 6 -> within 20%
            ("c", 600, 700, "ACGTACG", 14),  # pred period 25 -> outside
            ("c", 900, 950, "AT", 25),      # unmatched truth
        ])
        preds = _write_bed(tmp_path, "p.bed", [
            ("c", 0, 100, "AT", 50),
            ("c", 200, 300, "ACGACG", 16),
            ("c", 410, 500, "ACGTAC", 15),
            ("c", 600, 700, "ACGTACGTACGTACGTACGTACGTA", 4),
        ])
        out = str(tmp_path / "r.json")
        r = _run([truth, preds, "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["matched"] == 4
        assert res["sensitivity_1to1_maxcard"] == pytest.approx(80.0)
        assert res["precision_1to1_maxcard"] == pytest.approx(100.0)
        assert res["period"]["scored_pairs"] == 4
        assert res["period"]["exact_pct"] == pytest.approx(25.0)
        assert res["period"]["integer_multiple_pct"] == pytest.approx(25.0)
        assert res["period"]["within_20pct_pct"] == pytest.approx(25.0)
        assert res["period"]["outside_pct"] == pytest.approx(25.0)
        # boundary: three exact-start pairs, one 10 bp start offset
        assert res["boundary"]["exact_start_pct"] == pytest.approx(75.0)
        assert res["boundary"]["start_offset_median"] == pytest.approx(0.0)
        # copies: |50-50|/50, |16-33|/33, |15-20|/20, |4-14|/14
        assert res["copies"]["scored_pairs"] == 4

    def test_truth_without_motif_skips_period_and_copies(self, tmp_path):
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100, "AT", 50)])
        out = str(tmp_path / "r.json")
        r = _run([truth, preds, "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["matched"] == 1
        assert res["period"]["scored_pairs"] == 0
        assert res["copies"]["scored_pairs"] == 0
        assert list(res["strata"]) == ["other"]

    def test_zero_truth_copies_not_divided(self, tmp_path):
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100, "AT", 0)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100, "AT", 50)])
        out = str(tmp_path / "r.json")
        r = _run([truth, preds, "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["copies"]["scored_pairs"] == 0

    @pytest.mark.parametrize("bad", ["0", "-0.5", "1.5", "nan"])
    def test_min_overlap_validated(self, tmp_path, bad):
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100)])
        r = _run([truth, truth, "--min-overlap", bad])
        assert r.returncode != 0

    def test_inverted_interval_fatal(self, tmp_path):
        truth = _write_bed(tmp_path, "t.bed", [("c", 100, 100)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100)])
        r = _run([truth, preds])
        assert r.returncode != 0
        assert "empty/inverted" in r.stderr + r.stdout

    def test_bad_coordinates_fatal(self, tmp_path):
        truth = _write_bed(tmp_path, "t.bed", [("c", "x", 100)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100)])
        r = _run([truth, preds])
        assert r.returncode != 0

    def test_gzip_input(self, tmp_path):
        import gzip as _gzip
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100, "AT", 50)])
        gz = tmp_path / "p.bed.gz"
        with _gzip.open(gz, "wt") as f:
            f.write("c\t0\t100\tAT\t50\n")
        out = str(tmp_path / "r.json")
        r = _run([truth, str(gz), "--json", out])
        assert r.returncode == 0, r.stderr
        assert json.load(open(out))["matched"] == 1

    def test_pred_col5_period_disables_copy_metric(self, tmp_path):
        # converted ULTRA/tantan/TRF BEDs carry the PERIOD in column 5;
        # comparing it against a truth copy count would be a unit mismatch
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100, "AT", 50)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100, "AT", 2)])
        out = str(tmp_path / "r.json")
        r = _run([truth, preds, "--pred-col5", "period", "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["copies"]["scored_pairs"] == 0
        assert res["pred_col5"] == "period"
        # period metrics still score (col4 is a real motif)
        assert res["period"]["scored_pairs"] == 1

    def test_pred_motif_is_sequence_disables_period_metric(self, tmp_path):
        # TRF-style BED: col4 is the whole array sequence, so a
        # length-derived prediction period would be meaningless
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100, "AT", 50)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100, "AT" * 50, 2)])
        out = str(tmp_path / "r.json")
        r = _run([truth, preds, "--pred-col5", "period",
                  "--pred-motif-is-sequence", "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["matched"] == 1
        assert res["period"]["scored_pairs"] == 0
        # truth-based strata still fill
        assert res["strata"]["1-6"]["matched"] == 1

    def test_deep_contested_chain_no_recursion_error(self):
        """A half-shifted pileup forces one augmenting path spanning the
        whole chain; the recursive DFS died at Python's ~1000-frame default
        (reproduced at depth ~993). The iterative form must survive and
        still find the maximum matching (N pairs for N truths)."""
        n = 3000
        truth = [rec("c", i * 100, i * 100 + 100) for i in range(n)]
        # preds half-shifted, plus one extra at the front so augmenting
        # paths must propagate down the chain
        preds = [rec("c", i * 100 - 50, i * 100 + 50) for i in range(n + 1)]
        pairs, t_used, _ = m.one_to_one(truth, preds, 0.3)
        assert len(pairs) == n
        assert t_used == set(range(n))

    def test_within20_rule_matches_permissive_matcher(self, tmp_path):
        # smaller-period-relative, like periods_compatible(): truth 10 vs
        # pred 8 is 2/8 = 25% -> OUTSIDE, not within-20%
        truth = _write_bed(tmp_path, "t.bed", [("c", 0, 100, "ACGTACGTAC", 10)])
        preds = _write_bed(tmp_path, "p.bed", [("c", 0, 100, "ACGTACGT", 12)])
        out = str(tmp_path / "r.json")
        r = _run([truth, preds, "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["period"]["within_20pct_pct"] == pytest.approx(0.0)
        assert res["period"]["outside_pct"] == pytest.approx(100.0)

    def test_truth_chroms_only_drops_scaffold_preds(self, tmp_path):
        truth = _write_bed(tmp_path, "t.bed", [("chr1", 0, 100, "AT", 50)])
        preds = _write_bed(tmp_path, "p.bed", [
            ("chr1", 0, 100, "AT", 50),
            ("scaffold_7", 0, 100, "AT", 50),
        ])
        out = str(tmp_path / "r.json")
        # without the flag the scaffold call deflates precision
        r = _run([truth, preds, "--json", out])
        assert json.load(open(out))["precision_1to1_maxcard"] == pytest.approx(50.0)
        r = _run([truth, preds, "--truth-chroms-only", "--json", out])
        assert r.returncode == 0, r.stderr
        res = json.load(open(out))
        assert res["precision_1to1_maxcard"] == pytest.approx(100.0)
        assert res["dropped_off_truth_chroms"] == 1
        assert res["pred_records_loaded"] == 2


DERIVE = os.path.join(os.path.dirname(__file__), "..",
                      "scripts", "scoring", "derive_adotto_annotated_truth.py")


class TestDeriveAnnotatedTruth:
    def _catalog_line(self, chrom, r_start, r_end, anns):
        return "\t".join([chrom, str(r_start), str(r_end)] +
                         ["."] * 14 + [json.dumps(anns)]) + "\n"

    def test_annotation_coords_and_primitive_reduction(self, tmp_path):
        cat = tmp_path / "cat.bed"
        cat.write_text(self._catalog_line("chr1", 9975, 10498, [
            {"start": 10000, "end": 10467, "score": 1243,
             "motif": "ATGATGATG", "copies": 3.7},   # 9-mer = ATG x3
            {"start": 10481, "end": 10497, "score": 41,
             "motif": "GCCC", "copies": 4.2},
        ]))
        r = subprocess.run([sys.executable, DERIVE, str(cat)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        fields = r.stdout.strip().split("\t")
        # annotation coords, not the padded region's
        assert fields[:3] == ["chr1", "10000", "10467"]
        # primitive unit, copies rescaled by the reduction factor
        assert fields[3] == "ATG"
        assert float(fields[4]) == pytest.approx(11.1)
        assert "primitive_reduced=1" in r.stderr

    def test_restrict_to_region_set(self, tmp_path):
        cat = tmp_path / "cat.bed"
        cat.write_text(
            self._catalog_line("chr1", 0, 100, [
                {"start": 10, "end": 90, "score": 5, "motif": "AT", "copies": 40}])
            + self._catalog_line("chr2", 0, 100, [
                {"start": 10, "end": 90, "score": 5, "motif": "AT", "copies": 40}]))
        regions = tmp_path / "primary.bed"
        regions.write_text("chr1\t0\t100\n")
        r = subprocess.run([sys.executable, DERIVE, str(cat),
                            "--restrict-to", str(regions)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        lines = r.stdout.strip().splitlines()
        assert len(lines) == 1 and lines[0].startswith("chr1\t")

    def test_unexpected_annotation_shape_is_motifless_not_crash(self, tmp_path):
        cat = tmp_path / "cat.bed"
        cat.write_text(self._catalog_line("chr1", 0, 100, [None, "bare"]))
        r = subprocess.run([sys.executable, DERIVE, str(cat)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        fields = r.stdout.rstrip("\n").split("\t")
        assert fields[:3] == ["chr1", "0", "100"]  # region coords, no motif
        assert fields[3] == ""
