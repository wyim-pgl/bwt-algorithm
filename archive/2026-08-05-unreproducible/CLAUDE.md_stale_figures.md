# Removed from CLAUDE.md on 2026-08-05 — unreproducible

See ./README.md for why. Kept verbatim for provenance only.

> **The recall/precision figures below are from the pre-`706fb76` build** (which read
> uninitialised heap in its alignment traceback and was not reproducible). All three
> species were re-measured on the fixed build `d52a4ff` (2026-07-10, in
> `docs/2026-06-25-catch-all-benchmark-for-filip.md`): the fix trades ~0.1 pp region
> recall for ~0.4 pp region precision and leaves every satellite/centromere number
> unchanged, so the op-points and their ordering hold. Refreshed genome figures:
> human catchF **80.54 % / 50.52 %** (adj 79.5 %), catchH **82.23 % / 48.35 %**.
> See `docs/2026-07-09-nondeterminism-uninitialised-ptr-table.md`.
