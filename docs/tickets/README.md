# alpha2 tickets

Derived from ADRs 0001, 0003, 0004. Read `CONTEXT.md` and those ADRs before starting — they carry the reasoning, these carry the work.

| # | Ticket | Depends on |
|---|---|---|
| 0001 | [Nested cache hierarchy](0001-nested-cache-hierarchy.md) | — |
| 0002 | [Leaf uniqueness validation](0002-leaf-uniqueness-validation.md) | — |
| 0003 | [MATCH per chromosome](0003-match-stage-per-chromosome.md) | 0001 |
| 0004 | [MATCH_COMBINE and gate](0004-match-combine-and-gate.md) | 0003 |
| 0005 | [Cache shape invariants](0005-cache-shape-invariants.md) | 0001, 0003 |
| 0006 | [Absence semantics declaration](0006-absence-semantics-declaration.md) | — |
| 0007 | [SCORE per shard](0007-score-stage-per-shard.md) | 0001, 0004 |
| 0008 | [SCORE_COMBINE](0008-score-combine.md) | 0007 |
| 0009 | [Imputation floor guard](0009-impute-floor-guard.md) | — |
| 0010 | [Match rate quantisation](0010-match-rate-quantisation.md) | — |
| 0011 | [`match_flipped` always zero](0011-match-flipped-always-zero.md) | — |
| 0012 | [Give the Cache a module](0012-genotype-cache-module.md) | — |
| 0013 | [Fold the cache traversals in](0013-fold-cache-traversals.md) | 0012 |
| 0014 | [Build the leaf path once](0014-build-the-leaf-path-once.md) | 0012 |
| 0015 | [Delete the multi-shard dosage machinery](0015-delete-multi-shard-dosage-machinery.md) | — |
| 0016 | [One connection seam for the score database](0016-score-db-connection-seam.md) | — |
| 0017 | [Split `ScorePipeline` per stage](0017-split-scorepipeline-per-stage.md) | 0012, 0016 |
| 0018 | [SCORE must reject an ungated match table](0018-gated-match-table-marker.md) | — |

Grab first: 0001 unblocks most of the pipeline. 0002, 0006 and 0009 are independent and can go in parallel.

0010 and 0011 are pre-existing match log bugs found while implementing 0004. Both were left alone there because fixing them changes match log output, which would have invalidated 0004's "logs match today's output" acceptance criterion.

0012 – 0017 came out of an architecture review of the alpha2 implementation, after 0001 – 0009 landed. 0012 unblocks the rest of that group; 0015 and 0016 are independent and can go in parallel. 0017 needed a decision before it started — where scratch goes relative to published artefacts — now recorded as [ADR 0005](../adr/0005-scratch-lives-under-out-dir.md), along with [ADR 0006](../adr/0006-one-module-per-stage-no-base-class.md) for the shape of the split it produced.

## Not in this repo

Nextflow module wiring lives in `pgsc_calc` (`workflows/pgsc_calc/main.nf`, `modules/local/`): four processes replacing the single `PGSC_CALC_SCORE` task, chromosome-grouped channels for MATCH, shard-level channels for SCORE.

A pre-flight validation process reading VCF headers and bgen `.sample` files would fail faster than 0005 does — an invariant breach at MATCH_COMBINE means every LOAD task has already run, which over a biobank is hours. Purely a compute optimisation; the authoritative check still belongs at MATCH.

## Deferred

`homref` absence semantics (sparse-called WGS, ragged batches) — design recorded in ADR 0004, strictly additive on top of alpha2.

Reference-database backfill for variants absent from the whole sampleset — needs its own ADR, because it decouples "matched" from "has a dosage" and would otherwise drive the match rate to 100% by construction and kill `min_overlap`.

One finding from the same architecture review is still unticketed: deciding whether a run has one sampleset. It is asserted three ways and iterated once, with a `_matchlog.py` TODO recording that the multi-sampleset path is untested. The interchangeable match tables from that same review are now 0018.
