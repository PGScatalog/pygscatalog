# alpha2 tickets

Open work only. A ticket is deleted when it lands — the commit that closes it is its record, and `git log -- docs/tickets/` finds the ones that are gone. Decisions that outlive the work belong in `docs/adr/`, not here.

Derived from ADRs 0001, 0003, 0004. Read `CONTEXT.md` and those ADRs before starting — they carry the reasoning, these carry the work.

| # | Ticket | Depends on |
|---|---|---|
| 0010 | [Match rate quantisation](0010-match-rate-quantisation.md) | — |
| 0011 | [`match_flipped` always zero](0011-match-flipped-always-zero.md) | — |

Both are pre-existing match log bugs, found while splitting MATCH_COMBINE out. Both were left alone there because fixing them changes match log output, which would have invalidated that work's "logs match today's output" acceptance criterion. They are independent of each other.

## Not in this repo

Nextflow module wiring lives in `pgsc_calc` (`workflows/pgsc_calc/main.nf`, `modules/local/`): four processes replacing the single `PGSC_CALC_SCORE` task, chromosome-grouped channels for MATCH, shard-level channels for SCORE.

A pre-flight validation process reading VCF headers and bgen `.sample` files would fail faster than the cache shape invariants do — an invariant breach at MATCH_COMBINE means every LOAD task has already run, which over a biobank is hours. Purely a compute optimisation; the authoritative check still belongs at MATCH.

## Deferred

`homref` absence semantics (sparse-called WGS, ragged batches) — design recorded in ADR 0004, strictly additive on top of alpha2.

Reference-database backfill for variants absent from the whole sampleset — needs its own ADR, because it decouples "matched" from "has a dosage" and would otherwise drive the match rate to 100% by construction and kill `min_overlap`.

One finding from the alpha2 architecture review is still unticketed: deciding whether a run has one sampleset. It is asserted three ways and iterated once, with a `_matchlog.py` TODO recording that the multi-sampleset path is untested.
