# 0004 — MATCH_COMBINE: pool match counts and gate on min_overlap

**Depends on:** 0003
**Blocks:** 0007
**Refs:** ADR 0003

## Why

No single chromosome's MATCH task can decide `min_overlap` — the rate is genome-wide. Pooling has to happen once, centrally, and its output is the accession list that SCORE fans out over.

## Scope

- New CLI: `pgscatalog-calc match-combine`, consuming every chromosome's partial match table.
- Pool per `(sampleset, accession)` and apply `min_overlap`. **Not** per shard, and **not** per sample batch: a technical batch is a slice of one callset, so a match rate that differs between batches is a data-quality alarm, not grounds for scoring some samples and not others.
- Compute the per-batch rate anyway and log a warning when batches diverge beyond a tolerance. Alarm, not gate — a gate here yields an accession scored for 300k samples and absent for 100k with nothing in the score file explaining why.
- Verify the sample union is identical across chromosomes. Scope this to chromosomes **present in the cache**, not chromosomes named in the samplesheet — `cleanup_empty_zarr_group` (`load.py:222`) legitimately deletes a group whose index query returned nothing, which happens whenever the scoring files cover no variants on a chromosome. Phrased against the samplesheet, this check fails spuriously on any chromosome-sparse scoring file.
- Emit the passing accession list plus the existing match logs.

## Acceptance

- An accession at 74% genome-wide is rejected even where one chromosome is at 90%.
- Divergent per-batch rates warn and do not alter which accessions pass.
- A scoring file covering only chr1–chr5 does not trip the cross-chromosome union check.
- Match logs match today's output for a single-shard sampleset.

## Deferred in implementation

**The per-batch match rate divergence warning is not implemented.** It would be dead
code by construction: 0005 rejects ragged shards at MATCH, so by the time
MATCH_COMBINE runs, every batch under a chromosome has been proven to hold the same
variants and therefore the same match rate. The warning becomes meaningful only when
`homref` semantics make ragged batches legal, and it belongs with that work.
