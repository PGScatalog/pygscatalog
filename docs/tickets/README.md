# alpha2 tickets

This directory contains open work only.
Delete a ticket when its implementation is complete.
The closing commit preserves the ticket's record.
Use `git log -- docs/tickets/` to find deleted tickets.
Record decisions that remain relevant after implementation in `docs/adr/`.

The tickets derive from ADRs 0001, 0003, 0004 and 0007–0014.
Read `CONTEXT.md` and those ADRs before work starts.
The glossary and ADRs explain the decisions. Tickets describe the required work.

## Open

No tickets are open.

## Completed work

The latest completed work was the Sex check in
[ADR 0014](../adr/0014-declared-sex-is-checked-never-used.md).
It amends ADR 0007 without reversing its ploidy decision.
The pipeline checks Declared sex but never uses it to determine ploidy.
Ploidy remains observed per call, and the repository contains no PAR coordinates.

- Ticket 0039 added the declaration and its requirement, which tests can check without genotypes.
- Ticket 0040 added the comparison. A mismatch raises an error. Declaring `0` exempts a sample.

The preceding work added Serial run in ADR 0013.
It provides an entry point for environments where Nextflow cannot run.
It does not add another scaling mechanism. `pgsc_calc` remains the supported path.

- Ticket 0036 added `pgsc_calc run`.
- Ticket 0037 changed integration tests to use that command. This tests the public entry point directly.
- Ticket 0038 measured the memory property that makes serial execution usable with real data.

Earlier work proceeded in this order:

1. Tickets 0019–0025 added X, Y and MT support.
2. That work exposed two defects, recorded as tickets 0026 and 0027.
3. ADR 0010 introduced multiallelic handling and tickets 0028–0032.
4. That implementation exposed tickets 0033 and 0034. Both described crashes or incorrect results on `main`, independently of multiallelic support.
5. Ticket 0035 added v2's `--ignore_strand_flips` option.

Two decisions from that work remain documented:

- A wholly missing VCF call retains the ploidy declared by its GT. A guard applies when no call at the variant is hemizygous.
  See [ADR 0007](../adr/0007-ploidy-is-observed-per-call.md).
- A strand flip uses a reverse complement. Indels do not support that inference.
  See [ADR 0011](../adr/0011-strand-flips-are-reverse-complements-and-optional.md).

## Work outside this repository

Nextflow module configuration belongs in `pgsc_calc`, under `workflows/pgsc_calc/main.nf` and `modules/local/`.
Four processes replace the single `PGSC_CALC_SCORE` task.
MATCH channels group inputs by chromosome. SCORE channels group inputs by shard.

A validation process could read VCF headers and bgen `.sample` files before LOAD.
It would detect invalid layouts earlier than the cache checks.
An error at MATCH_COMBINE occurs after every LOAD task completes, which can take hours for a biobank.
This proposed process would reduce wasted computation.
The authoritative check would remain at MATCH.

## Deferred work

`homref` absence semantics supports sparse-called WGS and ragged batches.
ADR 0004 records its design as an addition to alpha2.

Using a reference database for variants absent from the whole sampleset needs a separate ADR.
It separates a match from the availability of a dosage.
Without a separate decision, such matches could force the match rate to 100% and make `min_overlap` ineffective.

`_score_shard` in `tests/integration/test_bgen.py` defines the stage sequence for a third time.
Ticket 0037 deliberately retained it.
It receives a cache zip instead of target genomes.
Replacing it with `pgsc_calc run` would require that command to accept an existing cache.
That is a subcommand decision, not only a test cleanup.

The alpha2 architecture review left one finding without a ticket: whether a run contains one sampleset.
The code asserts this in three places and iterates over samplesets in another.
A TODO in `_matchlog.py` records that the path for multiple samplesets lacks tests.

The alpha2 plan also left chrX comparisons with plink2 without tickets and planned integration tests manually.
ADR 0007 based dosage 1 for a hemizygous effect allele on an assumption about the default for plink2 `--score`.
Tickets 0020–0025 did not check that assumption against plink2.
Ticket 0020 checked dosage against the target VCF calls, which is a different check.
