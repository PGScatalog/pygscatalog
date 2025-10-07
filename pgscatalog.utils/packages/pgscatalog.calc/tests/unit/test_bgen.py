import shutil

import numpy as np
import pytest

from pgscatalog.calc import TargetVariants
from pgscatalog.calc.lib._bgen import (
    bgen_buffer_variants,
    bgen_get_sample_list,
)
from pgscatalog.calc.lib._bgen_probabilities import (
    phased_probabilities_to_hard_calls,
    unphased_probabilities_to_hard_calls,
)


def test_phased_bgen(
    test_positions,
    phased_bgen_path,
    bgen_sample,
    n_samples,
    tmp_path_factory,
    sample_ids,
):
    variants = list(
        bgen_buffer_variants(
            position_batch=test_positions,
            target_path=phased_bgen_path,
            sample_path=bgen_sample,
        )
    )

    samples = bgen_get_sample_list(
        target_path=phased_bgen_path, sample_path=bgen_sample
    )

    assert len(variants) == len(test_positions)
    assert len(samples) == n_samples

    targets = TargetVariants(
        variants, samples=sample_ids, filename=str(phased_bgen_path), sampleset="test"
    )
    gts = targets.genotypes
    assert gts.shape == (3, 3202, 2)
    assert np.ptp(gts) == np.uint8(1)


def test_unphased_bgen(
    test_positions,
    unphased_bgen_path,
    bgen_sample,
    n_samples,
    tmp_path_factory,
    sample_ids,
):
    variants = list(
        bgen_buffer_variants(
            position_batch=test_positions,
            target_path=unphased_bgen_path,
            sample_path=bgen_sample,
            is_chromosome_padded=False,
        )
    )

    assert len(variants) == len(test_positions)

    targets = TargetVariants(
        variants, samples=sample_ids, filename=str(unphased_bgen_path), sampleset="test"
    )
    gts = targets.genotypes
    assert gts.shape == (3, 3202, 2)
    assert np.ptp(gts) == np.uint8(1)


def test_dosage_equality(
    test_positions,
    unphased_bgen_path,
    phased_bgen_path,
    bgen_sample,
    n_samples,
    tmp_path_factory,
):
    unphased_variants = list(
        bgen_buffer_variants(
            position_batch=test_positions,
            target_path=unphased_bgen_path,
            sample_path=bgen_sample,
            is_chromosome_padded=False,
        )
    )

    phased_variants = list(
        bgen_buffer_variants(
            position_batch=test_positions,
            target_path=phased_bgen_path,
            sample_path=bgen_sample,
            is_chromosome_padded=False,
        )
    )

    for unphased, phased in zip(unphased_variants, phased_variants, strict=True):
        # sample dosage must be the same
        # but precise gt order for heterozygous (0, 1) or (1, 0) can differ
        # (because it's inferred in unphased samples)
        assert np.array_equal(np.sum(unphased.gts, axis=1), np.sum(phased.gts, axis=1))


def test_bgen_with_missing(
    test_positions, unphased_bgen_path, bgen_sample, tmp_path_factory
):
    tmp_path_factory.mktemp("cache")
    missing_chrom, missing_pos = "1", 9067156  # definitely missing
    positions = test_positions + [(missing_chrom, missing_pos)]

    variants = list(
        bgen_buffer_variants(
            position_batch=positions,
            target_path=unphased_bgen_path,
            sample_path=bgen_sample,
            is_chromosome_padded=False,
        )
    )

    # test that missing variants are yielded properly
    assert len(variants) == len(positions)
    missing_variant = variants[-1]
    assert missing_variant.gts is None
    assert missing_variant.ref is None
    assert missing_variant.alts is None
    assert missing_variant.chr_name == missing_chrom
    assert missing_variant.chr_pos == missing_pos


def test_bgen_no_index(test_positions, phased_bgen_path, bgen_sample, tmp_path_factory):
    tmp_path_factory.mktemp("cache")
    outdir = tmp_path_factory.mktemp("bgen_no_index")
    shutil.copy(phased_bgen_path, outdir)

    with pytest.raises(ValueError) as e:
        _ = list(
            bgen_buffer_variants(
                target_path=outdir / "tiny1000G.bgen",
                sample_path=bgen_sample,
                position_batch=test_positions,
                is_chromosome_padded=False,
            )
        )
        assert "index" in str(e.value)


@pytest.fixture
def unphased_probabilities():
    # homozygous first allele (0, 0)
    # heterozygous (0, 1)
    # homozygous second allele (1, 1)
    probs = np.array([[0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.0, 0.3, 0.7]])
    gts = np.array([[0, 0], [0, 1], [1, 1]], dtype=np.uint8)
    return probs, gts


@pytest.fixture
def phased_probabilities():
    # hap 1 - allele 1, hap 2 - allele 2 (0, 0)
    # hap 1 - allale 1, hap 2 - allele 2 (0, 1)
    # hap 1 - allele 2, hap 1 - allele 1 (1, 0)
    # hap 1 - allele 2, hap 2 - allele 2 (1, 1)
    probs = np.array(
        [
            [0.7, 0.3, 0.7, 0.3],
            [0.8, 0.2, 0.1, 0.9],
            [0.25, 0.75, 0.8, 0.2],
            [0.0, 1.0, 0.1, 0.9],
        ]
    )

    gts = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=np.uint8)
    return probs, gts


def test_phased_probs(unphased_probabilities):
    probs, gts = unphased_probabilities
    assert np.array_equal(unphased_probabilities_to_hard_calls(probs), gts)


def test_unphased_probs(phased_probabilities):
    probs, gts = phased_probabilities
    assert np.array_equal(phased_probabilities_to_hard_calls(probs), gts)


def test_bad_input(unphased_probabilities, phased_probabilities):
    unphased_probs, _ = unphased_probabilities
    phased_probs, _ = phased_probabilities

    with pytest.raises(ValueError):
        phased_probabilities_to_hard_calls(unphased_probs)

    with pytest.raises(ValueError):
        unphased_probabilities_to_hard_calls(phased_probs)
