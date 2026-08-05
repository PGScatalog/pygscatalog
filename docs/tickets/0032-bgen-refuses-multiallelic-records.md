# 0032 — A multiallelic bgen record is skipped silently, and `ref_first` corrupts dosage

**Refs:** ADR 0010.

## Why

Two independent problems in the bgen reader, both about allele identity, both cheap to
close together.

**Skipped, not reported.** `bgen_read_variants` skips multiallelic records
(`lib/cache/_bgen.py:283-285`) and the position is then padded as absent, filtered out by
`update_match_table`, and logged as an ordinary unmatched variant — the same
misreporting ticket 0029 fixes for VCF. Once 0029 and 0030 land, identical records behave
differently depending on input format, with nothing in the match log saying which
happened.

ADR 0010 rules bgen out of multiallelic support deliberately: `bgen_to_hardcall` is
hard-wired biallelic (`probs[haploid, :2]` and a three-column unphased layout, lines
164-170), a triallelic unphased diploid call carries six probabilities in the spec's
colex ordering, and bgen input here is population-scale biobank data from standard
imputation pipelines that split multiallelics upstream. The decision is not to build the
decoder. It is not to skip quietly.

**`ref_first` swaps labels without swapping indices.** `ref` and `alt` are chosen by the
`ref_first` parameter (lines 254, 297-299), but `bgen_to_hardcall` returns argmax indices
in the bgen's own allele order. Set `ref_first=False` and index 0 still means
`alleles[0]`, while the cache now labels `alleles[1]` as REF — so a `REFALT` match with
`effect_allele_idx=0` counts the alt. Nothing passes the parameter, no caller can reach
it, no test covers it. It is dead speculative flexibility whose only reachable behaviour
is a silently wrong score.

## Scope

Let the existing `NotImplementedError` (lines 286-287) fire on the skip path too, naming
the variant, and delete `skip_multiallelic`. A bgen holding a multiallelic record at a
scoring-file position refuses to load rather than quietly scoring a subset. The message
should say what the user can do: convert to VCF, or split the record upstream.

Delete `ref_first` and read `ref` as `alleles[0]` unconditionally.

## Out of scope

Multiallelic bgen support — the colex decoder, per ADR 0010. If it is ever wanted, it
starts here, but it is not wanted now.

Non-SNP bgen records, which `is_snp` already skips and which this ticket does not change.

## Acceptance

- A bgen holding a multiallelic record at a scoring-file position fails at LOAD with a
  message naming the variant and the position, not a log line.
- A bgen with no multiallelic records at scoring-file positions loads exactly as it does
  today.
- `skip_multiallelic` and `ref_first` are gone.
- The format support table in ADR 0010 still describes what the code does.
