#!/usr/bin/env python3
"""Untracked regen-safe entry point for score_table1_reciprocal.py."""
import os, runpy
if __name__ == "__main__":
    os.environ["BWT_TABLE1_ORIGINAL"] = "score_table1_reciprocal.py"
    runpy.run_path(os.path.join(os.path.dirname(__file__), "score_table1_regen.py"), run_name="__main__")
