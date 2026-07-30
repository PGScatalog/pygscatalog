# 0016 — One connection seam for the score database

**Blocks:** 0017

## Why

Thirteen sites open a duckdb connection, under two competing conventions: take a `conn`
(`update_match_table`, `add_match_log`, the `check_*` functions) or take a `db_path` and open
your own (`create_wide_weights_table`, `get_variant_metadata`, `get_weight_matrix`,
`get_n_score_variants`, `create_score_table`, `insert_scores`, `get_ok_accessions`,
`calculate_score_statistics`, and the accession query in `calculate_scores`).

Nine of the thirteen open with no configuration, so `--max_memory_gb` and `--threads` never
reach them. Only `ScorePipeline._connect`, `load_scoring_files`, `write_scores` and the two
log exports apply the limits the CLI advertises.

This matters more in a TRE than locally. duckdb's default is a share of *visible* RAM, which
inside a scheduled task is usually the host's rather than the task's limit, so exceeding it
gets the task killed by the scheduler rather than spilled to disk. And the query it matters
most for is the `PIVOT` in `create_wide_weights_table` — the memory-hungriest step in SCORE,
at the site of the alpha1 out-of-memory failure (pgsc_calc#459) that alpha2 exists to fix.
It is currently unconfigured.

## Scope

- A module owning the database: it holds the path and the resource policy and hands out
  read-only or writable connections.
- Every caller takes a connection. One convention, not two.

## Out of scope

Changing the defaults or tuning duckdb — this ticket makes the declared policy reach every
query, it does not decide what the policy should be. Splitting `ScorePipeline` (0017).

## Acceptance

- Every connection to `pgs.db` carries `max_memory` and `threads`, asserted by a test that
  reads the setting back from inside a running stage rather than by grep alone.
- No function in `lib/score` takes a `db_path` and opens its own connection.
- Unit tests can pass an in-memory connection instead of a temp file, and at least one test
  that currently needs a database on disk stops needing one.
- Scores unchanged.
