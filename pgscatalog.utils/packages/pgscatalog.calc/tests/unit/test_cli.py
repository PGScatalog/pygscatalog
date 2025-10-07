from pgscatalog.calc.cli.load import AnnotatedGenomePath


def test_annotated_genome_path():
    genome = AnnotatedGenomePath.from_string(
        "s3://1000genomes/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz:12"
    )
    assert (
        genome.path
        == "s3://1000genomes/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    )
    assert genome.chrom == "12"

    genome = AnnotatedGenomePath.from_string("path/to/file_with_merged_chroms.vcf.gz")
    assert genome.chrom is None
    assert genome.path == "path/to/file_with_merged_chroms.vcf.gz"
