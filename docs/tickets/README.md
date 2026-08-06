# alpha2 tickets

Open work only. A ticket is deleted when it lands — the commit that closes it is its record, and `git log -- docs/tickets/` finds the ones that are gone. Decisions that outlive the work belong in `docs/adr/`, not here.

Derived from ADRs 0001, 0003, 0004, and 0007–0010. Read `CONTEXT.md` and those ADRs before starting — they carry the reasoning, these carry the work.

## Open

X, Y and MT support (0019–0025) has landed, along with 0026, the ploidy bug it uncovered. Multiallelic handling (ADR 0010) adds 0028–0032; 0027, 0028, 0029, 0032 and 0033 have landed.

| # | Ticket | Depends on |
|---|---|---|
| 0030 | [Only the first alternative allele can match](0030-match-every-alternative-allele.md) | 0029 ✅ |
| 0031 | [`EXCLUDED` means three different things](0031-every-exclusion-names-its-reason.md) | 0027 ✅, 0030 |

0029 and 0030 are the pair that has to ship together to be observable: 0029 caches records nothing can match yet, 0030 matches records that aren't cached yet.

0033 was why 0028 was not finished in practice: 0028 made the matcher skip `*` and symbolic alleles, but the cache writer refused those records first, so the criterion was only ever met by driving the matcher directly. It was a crash on `main` independently of multiallelic support.

## Not in this repo

Nextflow module wiring lives in `pgsc_calc` (`workflows/pgsc_calc/main.nf`, `modules/local/`): four processes replacing the single `PGSC_CALC_SCORE` task, chromosome-grouped channels for MATCH, shard-level channels for SCORE.

A pre-flight validation process reading VCF headers and bgen `.sample` files would fail faster than the cache shape invariants do — an invariant breach at MATCH_COMBINE means every LOAD task has already run, which over a biobank is hours. Purely a compute optimisation; the authoritative check still belongs at MATCH.

## Deferred

`homref` absence semantics (sparse-called WGS, ragged batches) — design recorded in ADR 0004, strictly additive on top of alpha2.

Reference-database backfill for variants absent from the whole sampleset — needs its own ADR, because it decouples "matched" from "has a dosage" and would otherwise drive the match rate to 100% by construction and kill `min_overlap`.

One finding from the alpha2 architecture review is still unticketed: deciding whether a run has one sampleset. It is asserted three ways and iterated once, with a `_matchlog.py` TODO recording that the multi-sampleset path is untested.

Validating chrX scores against plink2 is deliberately unticketed — integration tests for this are planned manually. It is worth doing eventually: ADR 0007's decision that a hemizygous effect allele contributes a dosage of 1 rests on matching plink2 `--score`'s default, and nothing in 0020–0025 checks that claim against the reference implementation (0020 checks the dosage against the target VCF's own calls, which is a different thing).
