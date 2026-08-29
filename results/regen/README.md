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

Full SHA-256 hashes are recorded in `../manifest.sha256`; row-level provenance
and scorer hashes are recorded in `../manifest.tsv`.
