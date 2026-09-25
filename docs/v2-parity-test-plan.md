# Extend semantic parity coverage with pgsc_calc v2.3.0

Status: implemented. The work allocation below records the implementation plan.
The fixture's `EXTENDED.md` records the resulting compatibility contract.

## Implementation findings

- The unchanged v2 release enables mean imputation at 50 or more samples. The
  missingness fixture therefore has both a default reference and an explicitly
  configured no-mean-imputation reference. Integer imputation values 0, 1 and 2
  test dependence on the called donors without fractional dosage examples.
- v2 multi-score statistics and aggregate denominators are not consistently
  accession-specific. The tests retain this evidence and independently verify
  v3 statistics instead of requiring equal averages under different denominators.
- The strand-exclusion reference requires the supported MATCH process argument
  override: the apparent top-level Nextflow parameter is not wired through v2.
- Zero-weight regression tests exposed coefficient-based membership filtering in
  v3. ADR 0015 defines the corrected observed-allele statistics. Existing golden
  statistics also require the corresponding correction; weighted sums do not.
- Cross-agent review added exact v2 union-statistic assertions and grouped MATCH
  inputs by chromosome for sample-batch tests. Runtime mutations of accession
  routing, sample identity, denominators and constant-one imputation are rejected.
- Validation: 164 relevant tests passed, with Ruff and formatting checks. The
  final source review identified an additional pre-existing boundary: different
  genetic models on the same target/effect allele across accessions can conflict
  with the weight-matrix primary key. This combination remains outside coverage.

## Objective and boundaries

Demonstrate agreement in selected alleles and calculated scores wherever v2 and
v3 share a scoring policy. Test documented v3 behaviour independently where they
intentionally differ. Byte-for-byte replication and preservation of v2 defects
are not goals.

Use deterministic, committed VCF hard calls and scoring files. Exclude fractional
dosage scenarios, including DS/GP inputs. Missing-call examples must use donor
genotypes whose imputed dosages are integers, so missingness can be tested without
introducing fractional dosage coverage. Decimal effect weights remain in scope.
PAR, heteroplasmy and ancestry normalisation remain excluded. Keep the existing
25-chromosome fixture and its references unchanged; new cases only need the
chromosomes relevant to their behaviour. MT retains standalone PLINK coverage.

## Test cases

| Workstream | Cases | Required assertions |
| --- | --- | --- |
| Arithmetic | Three additive accessions with overlapping and disjoint variants; the same target allele with different weights across accessions; different effect alleles of one target variant across accessions | Correct accession/sample keys, selected alleles, chromosome sums, dosage sums, denominators and final averages; no cross-score leakage |
| Arithmetic | Signed non-binary-exact decimal weights, mixed magnitudes and cancellation; deterministic variant ordering and enough rows to cross an actual scoring chunk boundary | Agreement within an explicit reference-precision budget; independent Decimal calculation for hard-call arithmetic; invariance to input order |
| Arithmetic | Additive, dominant and recessive effects on fully called diploid genotypes with dosages 0, 1 and 2 | Correct transformed contributions and aggregation of v2 effect-type shards; use hand-calculated expectations to interpret any denominator differences |
| Missingness | Different missing GT patterns between samples; missing REF- and ALT-effect calls; one sample missing every scored call but other samples providing called alleles | Observed allele counts differ between samples; sums follow the selected imputation policy; averages use the documented denominator; zero observed denominator produces the documented null result |
| Missingness | Entire variant uncalled; too few samples for v3 imputation | Explicit v3 failure consistent with its contract, not a fabricated zero score or a broad expected failure |
| Matching | Direct, swapped and complemented SNP matches; absent and incompatible alleles; missing other allele; ambiguous SNPs with keep/exclude settings; ignore-strand-flips setting | Selected target/effect identity and contributions, expected exclusions, score acceptance and absence of scores rejected by overlap policy |
| Matching | Repeated scoring rows resolving to the same target allele; distinct ALT effects at one multiallelic site; multiple target candidates for one score row | Collision exclusion and candidate selection follow ADR 0010; explicitly separate v2 agreements from intentional differences |
| Matching | Multi-base reverse complements and indels that would require strand inference | v3 follows ADR 0011, with hand-calculated selected contributions; v2 output records the intentional difference |
| Partitioning | Independently reordered sample columns between chromosomes; reversed chromosome and scoring-file input order; complete sample batches | Results keyed by sample/accession agree with the same committed unsplit reference; no positional joins, lost samples or duplicate output rows |
| Edge semantics | A zero-weight row alongside nonzero rows, and an all-zero accession | Sum is unchanged/zero; document and test whether matched zero-weight rows contribute to counts, dosage statistics and averages before declaring parity |

Do not combine every axis into one parameter product. Use small named cases that
isolate a failure. Partitioning starts with fully called data: splitting an
imputed cohort can change the allele-frequency estimation population and is not
automatically a valid invariance assertion.

## Policy work before implementation

1. Read CONTEXT.md and ADRs 0007, 0009, 0010, 0011 and 0014 alongside current code.
   Some ADR consequence sections describe historical defects; current tests and
   implementation determine which defects remain.
2. Record a case table distinguishing shared-policy parity, intentional divergence
   and expected input rejection. Existing unit tests already cover many matching
   rules and null averages; extend end-to-end coverage rather than copying them.
3. Missingness requires an explicit comparison: the recorded v2 command uses
   `no-mean-imputation`, while v3 imputes missing calls. Keep a diagnostic v2 run
   showing that difference. Only add a shared-policy reference if the unchanged
   release exposes a supported way to align it; do not patch v2 or silently
   change flags to obtain agreement. Otherwise assert v3's documented arithmetic
   independently and label this a tested difference, not a parity pass.
4. Use at least the actual imputation sample floor in the positive missingness
   fixture; do not bypass it with a monkeypatch. Keep missing Y calls out because
   VCF cannot express the required chromosome-absence distinction here.
5. Resolve zero-weight statistics explicitly. Do not adopt either implementation
   merely because it currently behaves that way. Propose matched observed alleles
   as the denominator basis, and document the decision before any production fix.
6. Retain explicit `--xchr-model 1` for shared sex-chromosome comparisons. Do not
   imply compatibility with all v2 defaults.

## Subagent execution and ownership

The coordinator owns shared files, the comparison contract, reference generation
and integration. Use at most three worker subagents concurrently. Every worker
must be told: you are not alone in the repository; do not revert others' changes;
edit only your owned files and adapt to the agreed shared interface.

All paths below are relative to
`pgscatalog.utils/packages/pgscatalog.calc/tests/`.

| Owner | Exclusive ownership | Deliverable |
| --- | --- | --- |
| Coordinator | Existing `integration/test_v2_parity.py`, shared helper if extraction is needed, integration `conftest.py`, fixture root README, nox session, this plan, ADR updates, every new reference-output directory | Freeze the minimal case/comparison interface; generate and verify reference outputs; review divergences; integrate and validate |
| Worker A: arithmetic | `data/v2_parity/arithmetic/` excluding reference-output directories; `integration/test_v2_parity_arithmetic.py` | Deterministic generator and committed inputs for multi-score, decimal-weight and diploid effect-type cases; behavioural tests |
| Worker B: missingness | `data/v2_parity/missingness/` excluding reference-output directories; `integration/test_v2_parity_missingness.py` | Integer-valued imputation examples, unequal denominators, zero-denominator sample and rejection cases; independent expected arithmetic |
| Worker C: matching | `data/v2_parity/matching/` excluding reference-output directories; `integration/test_v2_parity_matching.py` | Orientation, exclusion, collision and documented-divergence cases; expected selected target/effect identities |
| Worker A, second assignment | `integration/test_v2_parity_partitioning.py` | Reuse committed arithmetic inputs/reference for sample-order, input-order and complete-batch invariance |
| Worker B, second assignment | `data/v2_parity/zero_weight/` excluding reference-output directories; `integration/test_v2_parity_zero_weight.py` | Zero-weight cases after the coordinator records the policy decision |

Workers must not edit production code to force a parity pass. Return a minimal
reproducer for any new defect; the coordinator assigns a bounded fix separately
after checking the contract. No blanket xfails or weakened assertions.

## Dependency order

1. Coordinator freezes case names, accession allocation, comparison fields and
   ownership. Extract only the existing stage driver/reference reader that the
   new tests actually need; do not build a fixture framework.
2. Workers A/B/C prepare generators, deterministic inputs and explicit expected
   behaviour in parallel. Each supplies v2 invocation requirements and identifies
   any policy divergence before requesting reference generation.
3. Coordinator validates each input bundle and runs the unchanged v2.3.0 release
   at commit `72ee54f2119d2fe65ac671327342dfc412bd0a3b` using its pinned images.
   Use isolated work/output directories and one reference-generation owner.
4. Commit original per-chromosome/effect-type `.sscore` files, aggregate output,
   matching reports and provenance (commands, flags, versions, input hashes).
   Never derive references from v3. Workers finish their assertions against these
   concrete outputs; absent reference files must fail rather than skip.
5. Reassign available workers to partitioning and zero weights. Coordinator
   reviews all numerical differences and updates the coverage contract.
6. Run an adversarial review of the finished suite, with a concrete mutation check
   for weight/accession routing, sample identity and denominator handling.

## Comparison and completion criteria

- Enumerate expected files, sample/accession keys and matched rows; reject missing
  files and duplicate rows. Verify fixture hashes and semantic equivalence of
  v2/v3 scoring inputs in the tests.
- Compare chromosome scores, dosage statistics and counts before final sums and
  averages. Read v2 DENOM explicitly where it can differ from ALLELE_CT. Aggregate
  effect-type shards without losing or double-counting their contributions.
- Keep counts and selected allele identities exact. For decimal scores, derive
  and document the absolute/relative error budget from PLINK's saved precision,
  aggregation and v3's public rounding. Do not choose a tolerance from observed
  discrepancies. Check near-zero sums with an absolute bound.
- For documented differences, assert a specific v3 result or rejection with an
  independent biological/arithmetic explanation. Retain v2 evidence, but do not
  assert numerical equality or exact log wording where policies differ.
- Reuse existing committed references for partitioning tests. Reordered temporary
  inputs must be deterministic transformations of committed files.
- Normal integration runs launch no PLINK, Nextflow or Docker process. The
  generator cannot overwrite reference outputs.
- All new cases and existing parity tests pass; relevant matching, scoring and
  partitioning tests plus Ruff pass. Deliberately incorrect routing/counting is
  caught. The final coverage table lists agreements, intentional differences and
  exclusions, including fractional dosages.
