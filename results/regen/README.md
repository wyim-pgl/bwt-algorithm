# Regenerated benchmark evidence

These files support the manuscript rows regenerated from clean source commit
`0363d8b`.

- `regen_{colcen,human,maize}.provenance.json`: run configuration, inputs, and
  execution provenance for SLURM jobs 6110900, 6110901, and 6124640.
- `score_table1_regen_v2.txt`: human Table 1 scoring, including the separate
  tantan `-w 2000` extra row.
- `score_maize_regen.txt`: maize scores against the three satellite catalogs.
  This output uses the post-merge scorer and is not the convention used for
  manuscript Tables 3B and 3C.
- `table3bc_replacement.md` and `table3bc_provenance.json`: complete Table
  3B/3C replacement results under the legacy manuscript convention. The scorer
  first reproduced every old-BED manuscript cell before scoring the regenerated
  BED.
- `maize_extra_evidence.json`: machine-readable output for all regenerated
  maize values not covered by the reports above: the Table 3A/3A-b BWTandem
  row, the coordinate-only post-merge sweep quoted in Section 3.3.2, and the
  whole-genome 100--200-bp merged totals quoted in Section 3.3.3. It records
  the exact reproduction command, full SHA-256 and byte size of every input,
  and the scoring conventions. Table 3A was re-scored from the regenerated BED
  under both the historical mixed rule and the consistent 1--6-bp rule; it is
  not copied from the old table. The producer is
  `../../scripts/scoring/score_maize_regen_evidence.py`.

To reproduce the extra maize report on a system with the recorded external
inputs mounted, run the exact command stored in its top-level `command` field.
The output is deterministic: it contains no timestamp, and rerunning it at the
same output path must reproduce the deposited file byte for byte.

Full SHA-256 hashes are recorded in `../manifest.sha256`; row-level provenance
and scorer hashes are recorded in `../manifest.tsv`.
