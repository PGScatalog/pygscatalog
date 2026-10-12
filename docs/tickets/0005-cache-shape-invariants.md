# 0005 — Enforce sample disjointness and reject ragged batches

**Depends on:** 0001, 0003
**Refs:** ADR 0003, ADR 0004

## Why

Sample batching makes `GROUP BY (sampleset, accession, sample_id)` + `SUM` valid only while every `(sample, variant)` pair appears in exactly one shard. If two shards under one chromosome share a sample, that sample's dosage is double-counted — no error, no NaN, just an inflated score. This is a silent wrong answer, so it needs a real check rather than an assumption.

Checks live at MATCH because LOAD is fanned out per shard and can observe nothing cross-shard. They are free there: sample lists are already zarr group attrs (`targetvariants.py:190`), and nothing reads the 3D arrays.

## Scope

- Per-chromosome sample disjointness across leaves; their union is the sampleset's samples.
- Ragged detection: leaves under one chromosome disagreeing on which variants they contain. Under `unobserved` semantics (the only mode in alpha2 — see 0006) this is a hard error.
- **The ragged error must not tell the user to declare `homref`.** Ragged shards are ambiguous in cause: dropped monomorphic sites (where homozygous reference is true) and a re-imputed sample subset (where a variant may have failed QC, and homozygous reference is false) produce the same shape. Advise rectangularity — re-export without dropping monomorphic sites, or split on the variant axis. `homref` is a claim about how a callset was generated, not a remedy for an error message.
- Checks compare against the **cache**, which holds only positions the scoring files asked for. Narrower than the source files, and deliberately so: that is exactly the subset that can affect a score.
- Scope to chromosomes present in the cache (see 0004 for why).

## Out of scope

`get_sample_ids`' `NotImplementedError` (`_dosage.py:495`). It stops being reachable on its own — a SCORE task is staged with only its own leaf's subtree, so there is nothing left to disagree with. Leave it.

## Acceptance

- Two shards under one chromosome sharing a sample fail at MATCH, naming both.
- Ragged shards fail with a message that never mentions `homref`.
- Rectangular batches (bgen, or VCFs that kept monomorphic sites) pass.
- A chromosome absent from the cache entirely does not trip either check.
