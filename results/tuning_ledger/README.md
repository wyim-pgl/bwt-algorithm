# Tuning ledger

Deposited 2026-09-03 as the second half of the author's decision on
`quarantine.md` §3.9: disclose the configuration search rather than re-run it.

`ledger.tsv` is the append-only record of the coordinate-descent campaign that
chose BWTandem's configuration — 44 scored configurations on adotto chromosomes
21 and 22, with accept and reject decisions taken after recall and precision were
observed. `best.json` holds the selection state, including the objective in the
authors' own words: reaching ULTRA's recall at approximately ULTRA's precision.

Methods 2.2.3 states what this implies and the manuscript does not soften it:
every reported operating point is a post-selection, in-sample choice, and the
chromosome 22 result is post-selection validation rather than an independent test.

The catalog scored against is the same adotto v1.2.1 set the paper reports
against. Competitors received no equivalent search.
