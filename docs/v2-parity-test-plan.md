# Extend semantic parity coverage with pgsc_calc v2.3.0

Status: implemented. The work allocation below records the implementation plan.
The fixture's `EXTENDED.md` records the resulting compatibility contract.

## Implementation findings

- The unchanged v2 release enables mean imputation at 50 or more samples. The
  missingness fixture therefore has a default reference and an explicit
  no-mean-imputation reference. Integer dosages 0, 1 and 2 test imputation from
  called donors without fractional dosage examples.
- v2 multi-score statistics and aggregate denominators are not consistently
  accession-specific. Tests retain this evidence and check v3 statistics
  independently. They do not require equal averages when denominators differ.
- The strand-exclusion reference requires a supported argument override for the
  MATCH process. v2 does not forward the apparent top-level Nextflow parameter.
- Zero-weight regression tests exposed v3 membership filtering based on
  coefficients. ADR 0015 defines the corrected observed-allele statistics. Existing
  golden statistics require the same correction. Weighted sums do not change.
- Cross-agent review added exact assertions for v2 union statistics. It also
  grouped MATCH inputs by chromosome for sample-batch tests. Tests reject runtime
  mutations affecting accession routing, sample identity, denominators and
  constant-one imputation.
- Validation completed with 164 relevant tests passing, plus Ruff and formatting
  checks. Final source review identified another pre-existing limit. Different
  genetic models for the same target/effect allele across accessions can conflict
  with the weight-matrix primary key. This combination remains outside coverage.

## Objective and boundaries

Demonstrate agreement in selected alleles and calculated scores where v2 and v3
share a scoring policy. Test documented v3 behaviour independently where policies
intentionally differ. Byte-for-byte replication and preservation of v2 defects
are not goals.

Use deterministic, committed VCF hard calls and scoring files. Exclude fractional
dosage scenarios, including DS/GP inputs. Missing-call examples must use donors
whose imputed dosages are integers. This tests missingness without fractional
dosage coverage. Decimal effect weights remain in scope.

PAR, heteroplasmy and ancestry normalisation remain excluded. Keep the existing
25-chromosome fixture and its references unchanged. New cases need only
chromosomes relevant to their behaviour. MT retains standalone PLINK coverage.

## Test cases

| Workstream | Cases | Required assertions |
| --- | --- | --- |
| Arithmetic | Three additive accessions have overlapping and disjoint variants. Accessions assign different weights to the same target allele. They also use different effect alleles at one target variant. | Check accession/sample keys and selected alleles. Check chromosome sums, dosage sums, denominators and final averages. Exclude cross-score leakage. |
| Arithmetic | Signed decimal weights are not binary-exact. Weights vary in magnitude and include cancellation. Deterministic variant order crosses an actual scoring chunk boundary. | Apply an explicit reference-precision budget. Calculate hard-call arithmetic independently with Decimal. Check invariance to input order. |
| Arithmetic | Additive, dominant and recessive effects use fully called diploid genotypes with dosages 0, 1 and 2. | Check transformed contributions and combined v2 effect-type shards. Interpret denominator differences with hand-calculated expectations. |
| Missingness | Samples have different missing GT patterns. Missing calls include REF and ALT effects. One sample lacks every scored call, while other samples provide called alleles. | Observed allele counts must differ between samples. Sums must follow the chosen imputation policy. Averages must use the documented denominator. Zero observed denominators must produce the documented null result. |
| Missingness | An entire variant is uncalled. Another case has too few samples for v3 imputation. | Require the explicit v3 failure defined by its contract. Do not accept fabricated zero scores or broad expected failures. |
| Matching | SNP matches are direct, swapped or complemented. Cases include absent alleles, incompatible alleles and missing other alleles. Ambiguous SNPs use keep/exclude settings. Another setting ignores strand flips. | Check selected target/effect identities, contributions and exclusions. Check score acceptance. Scores rejected by overlap policy must be absent. |
| Matching | Repeated score rows resolve to the same target allele. Other rows select distinct ALT effects at one multiallelic site. One score row has multiple target candidates. | Require collision exclusion and candidate selection under ADR 0010. Separate v2 agreements from intentional differences. |
| Matching | Cases include multi-base reverse complements and indels that would require strand inference. | Require v3 behaviour under ADR 0011. Check selected contributions by hand. Retain v2 output as evidence of the intentional difference. |
| Partitioning | Chromosomes have independently reordered sample columns. Other runs reverse chromosome and scoring-file order or use complete sample batches. | Compare results by sample/accession with the same committed unsplit reference. Exclude positional joins, lost samples and duplicate rows. |
| Edge semantics | One accession combines zero and nonzero weights. Another accession has only zero weights. | Require an unchanged or zero sum. Document and test zero-weight membership in counts, dosage statistics and averages before declaring parity. |

Do not combine every axis into one parameter product. Use small named cases
that isolate a failure.

Partitioning starts with fully called data. Splitting an imputed cohort can
change the population used to estimate allele frequencies. Such a split does
not automatically support an invariance assertion.

## Policy work before implementation

1. Read CONTEXT.md and ADRs 0007, 0009, 0010, 0011 and 0014 with the current code.
   Some ADR consequences describe historical defects. Current tests and
   implementation determine which defects remain.
2. Record shared-policy parity, intentional divergence and expected input
   rejection in a case table. Existing unit tests cover many matching rules and
   null averages. Extend end-to-end coverage instead of copying those tests.
3. Compare missingness policies explicitly. The recorded v2 command uses
   `no-mean-imputation`, while v3 imputes missing calls. Keep a diagnostic v2 run
   showing this difference.

   Add a shared-policy reference only if the unchanged release provides a
   supported way to align policies. Do not patch v2. Do not silently change
   flags to obtain agreement. Otherwise, assert v3's documented arithmetic
   independently. Label this a tested difference, not a parity pass.
4. Use at least the actual imputation sample floor in the positive missingness
   fixture. Do not bypass it with a monkeypatch. Exclude missing Y calls. VCF
   cannot express the required chromosome-absence distinction here.
5. Resolve zero-weight statistics explicitly. Do not adopt either
   implementation merely because it behaves that way. Propose matched observed
   alleles as the denominator basis. Document the decision before changing
   production code.
6. Retain explicit `--xchr-model 1` for shared sex-chromosome comparisons. Do not
   imply compatibility with all v2 defaults.

## Subagent execution and ownership

The coordinator owns shared files, the comparison contract, reference generation
and integration. Use at most three worker subagents concurrently. Every worker
must receive these instructions:

- You are not alone in the repository.
- Do not revert other workers' changes.
- Edit only your owned files.
- Adapt to the agreed shared interface.

All paths below are relative to
`pgscatalog.utils/packages/pgscatalog.calc/tests/`.

| Owner | Exclusive ownership | Deliverable |
| --- | --- | --- |
| Coordinator | Existing `integration/test_v2_parity.py` and any necessary shared helper. Integration `conftest.py`, fixture root README, nox session, this plan and ADR updates. Every new reference-output directory. | Freeze the minimal case/comparison interface. Generate and check reference outputs. Review differences. Integrate and check the work. |
| Worker A: arithmetic | `data/v2_parity/arithmetic/`, excluding reference-output directories. `integration/test_v2_parity_arithmetic.py`. | Supply a deterministic generator and committed multi-score, decimal-weight and diploid effect-type inputs. Add behavioural tests. |
| Worker B: missingness | `data/v2_parity/missingness/`, excluding reference-output directories. `integration/test_v2_parity_missingness.py`. | Supply integer imputation examples, unequal denominators, a zero-denominator sample and rejection cases. Add independent arithmetic expectations. |
| Worker C: matching | `data/v2_parity/matching/`, excluding reference-output directories. `integration/test_v2_parity_matching.py`. | Supply orientation, exclusion, collision and documented-divergence cases. Define expected selected target/effect identities. |
| Worker A, second assignment | `integration/test_v2_parity_partitioning.py`. | Reuse arithmetic inputs and references for sample-order, input-order and complete-batch invariance tests. |
| Worker B, second assignment | `data/v2_parity/zero_weight/`, excluding reference-output directories. `integration/test_v2_parity_zero_weight.py`. | Add zero-weight cases after the coordinator records the policy decision. |

Workers must not change production code to force a parity pass. Return a minimal
reproducer for each new defect. The coordinator checks the contract before
assigning a separate, bounded correction. Do not add blanket xfails or weaken
assertions.

## Dependency order

1. The coordinator freezes case names, accession allocation, comparison fields
   and ownership. Extract only stage-driver or reference-reader code that the
   new tests need. Do not build a fixture framework.
2. Workers A/B/C prepare generators, deterministic inputs and explicit expected
   behaviour in parallel. Each worker supplies v2 invocation requirements.
   Each identifies policy differences before requesting reference generation.
3. The coordinator checks each input bundle. They run the unchanged v2.3.0
   release at commit `72ee54f2119d2fe65ac671327342dfc412bd0a3b` with its pinned
   images. Use isolated work/output directories. Assign reference generation to
   one owner.
4. Commit the original `.sscore` files for each chromosome and effect type.
   Commit aggregate output, matching reports and provenance. Provenance includes
   commands, flags, versions and input hashes. Never derive references from v3.
   Workers finish assertions against these outputs. Missing reference files must
   fail tests, not skip them.
5. Reassign available workers to partitioning and zero weights. The coordinator
   reviews all numerical differences and updates the coverage contract.
6. Review the completed suite adversarially. Test mutations of weight/accession
   routing, sample identity and denominator handling.

## Comparison and completion criteria

- Enumerate expected files, sample/accession keys and matched rows. Reject
  missing files and duplicate rows. Check fixture hashes and semantic
  equivalence of v2/v3 scoring inputs in tests.
- Compare chromosome scores, dosage statistics and counts before final sums and
  averages. Read v2 DENOM explicitly wherever it can differ from ALLELE_CT.
  Combine effect-type shards without losing or counting contributions twice.
- Keep counts and selected allele identities exact. Derive the absolute/relative
  decimal-score error budget from PLINK's saved precision, aggregation and v3's
  public rounding. Document that budget. Do not choose tolerances from observed
  discrepancies. Check near-zero sums with an absolute bound.
- For documented differences, assert a specific v3 result or rejection. Explain
  it independently through biology or arithmetic. Retain v2 evidence. Do not
  require numerical equality or identical log wording when policies differ.
- Reuse committed references for partitioning tests. Produce reordered temporary
  inputs through deterministic transformations of committed files.
- Normal integration runs start no PLINK, Nextflow or Docker process. The
  generator cannot overwrite reference outputs.
- All new cases and existing parity tests must pass. Relevant matching, scoring
  and partitioning tests must pass, together with Ruff. Tests must reject
  deliberately incorrect routing or counting. The final coverage table must
  list agreements, intentional differences and exclusions, including fractional
  dosages.
