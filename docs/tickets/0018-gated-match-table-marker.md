# 0018 — SCORE must reject an ungated match table

**Refs:** ADR 0006, 0004, 0017

## Why

Two artefacts are both "the match table", and they are structurally interchangeable:

- `matches.parquet`, written per chromosome by MATCH — `COPY allele_match_table`, every
  scoring variant it looked at
- `match_table.parquet`, written once by MATCH_COMBINE — the same table with
  `EXCLUDE (target_ref, target_alts)` and `WHERE accession IN $accessions`, so only the
  accessions that passed `min_overlap` survive

SCORE loads whichever it is given with `INSERT INTO allele_match_table BY NAME SELECT * FROM
read_parquet(...)`. `BY NAME` matches columns by name and ignores count and order, and the two
columns the gated table drops are nullable, so **both files load without error**. The docstring
on `import_match_tables` documents this as intended, because MATCH_COMBINE genuinely does
consume the ungated partials — the function is shared by both stages.

Hand SCORE `matches.parquet` and it scores every accession, including the ones that failed the
match-rate gate. No crash, no warning, no NaN: a complete-looking `scores.txt.gz` full of
plausible numbers. Those numbers are exactly what `min_overlap` exists to prevent — a score
computed from 20% of a scoring file's variants is not a reproduction of that published score,
it is a different number wearing its name.

This is not a live bug. The wiring in `pgsc_calc`'s `main.nf` is correct today. It is a
mistake nothing would catch, and it bites whoever edits the Nextflow channels, writes a new
test harness, or re-runs SCORE by hand with whatever parquet is nearest.

Only the match table has this problem. `samples.parquet` and `partial_scores.parquet` have
schemas that nothing else shares, so feeding a stage the wrong one already errors.

## Scope

Let the database refuse the wrong file, rather than adding a mechanism to detect it.

- SCORE's `allele_match_table` declares only the columns SCORE actually reads, which means
  dropping `target_ref` and `target_alts`. Verified: only `_matchlog.py` reads them, and that
  runs in MATCH_COMBINE. `INSERT ... BY NAME` then rejects `matches.parquet` on its own,
  because the source carries columns the target does not have — on duckdb 1.4.1,
  `Binder Error: Table "allele_match_table" does not have a column with name "target_ref"`.
- `create_allele_match_table` is shared, and MATCH_COMBINE genuinely needs those two columns
  (`add_match_log` writes them into `score_log_table`). So the DDL has to differ by stage:
  give the shared function an explicit argument, or write two. That difference is not
  incidental bookkeeping — it is precisely "what this stage consumes", which is why it is
  load-bearing enough to hang the rejection on.
- The binder error is loud but uninformative: it reads like a corrupt parquet rather than a
  wiring mistake. Catch it and re-raise naming both artefacts and saying that accessions which
  failed `min_overlap` would otherwise have been scored.
- MATCH_COMBINE keeps accepting MATCH's partial tables unchanged — that is its input.

This is deliberately less machinery than marking the file. Leave a `ponytail:` comment naming
the ceiling and the upgrade path:

> `ponytail:` the rejection is structural — SCORE declares only what it reads, so the binder
> refuses the ungated table. That holds only while the two tables' column sets differ. If SCORE
> ever needs `target_ref`, mark the gate explicitly instead: `COPY ... (FORMAT parquet,
> KV_METADATA {gated: 'true'})` in MATCH_COMBINE, read back with `parquet_kv_metadata()`
> (available in duckdb 1.4.1, and it annotates the file without touching any row schema).

**Do not** implement this as a check that inspects the file and decides — e.g. querying the
parquet schema and raising if `target_ref` is present. That reads a deliberate claim ("this
table has been gated") off an incidental fact ("these columns aren't needed downstream"), and
being a check rather than a constraint, it can be skipped, reordered or quietly deleted with
nothing failing. Declaring what SCORE consumes and letting the binder refuse the rest is not
an inference and cannot be bypassed.

## Out of scope

Marking or schema-checking the other stage artefacts. The same mechanism generalises to all
four, but the other three have distinct schemas and no demonstrated silent-failure mode, so
they can wait for a reason.

## Acceptance

- Feeding SCORE a `matches.parquet` fails before any scoring work, with an error naming both
  artefacts and stating that accessions which failed `min_overlap` would otherwise be scored —
  not a bare `BinderException`.
- SCORE still accepts `match_table.parquet` from MATCH_COMBINE, and the golden scores are
  unchanged for all four fixture shapes.
- MATCH_COMBINE still accepts every chromosome's `matches.parquet`.
- The test for the rejection drives the CLI, not just the loader, so it covers the wiring a
  real mistake would go through.
- SCORE's match table declares no column SCORE does not read. If that stops being true, the
  `ponytail:` ceiling has been reached and the explicit marker is the fix.
