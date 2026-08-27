"""A crashed stress-test case must fail the run, not vanish from the totals.

tests/test_random_stress.py is a script with a main(), not pytest test
functions, so the ordinary suite never exercises its pass/fail logic. This
drives it directly with an injected detector fault.
"""
import importlib.util
import pathlib

import pytest

STRESS = pathlib.Path(__file__).resolve().parent / "test_random_stress.py"


def _load():
    spec = importlib.util.spec_from_file_location("stress_mod", STRESS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_a_crashed_case_fails_the_run_even_when_the_rest_score_perfectly(capsys):
    mod = _load()
    mod.NUM_SEQUENCES = 3
    real = mod.run_finder
    seen = {"n": 0}

    def flaky(*a, **k):
        seen["n"] += 1
        if seen["n"] == 2:
            raise RuntimeError("injected fault")
        return real(*a, **k)

    mod.run_finder = flaky
    rc = mod.main()
    out = capsys.readouterr().out

    assert rc == 1, "a crashed case must fail the run"
    assert "1 crashed" in out
    assert "injected fault" in out
    # The surviving cases scored perfectly; the failure must not come from the
    # metric thresholds.
    assert "Sensitivity=100.0%" in out


def test_a_clean_run_still_passes():
    mod = _load()
    mod.NUM_SEQUENCES = 2
    assert mod.main() == 0
