How to match variants across reference panels and target genomes
=================================================================

``pgscatalog-intersect`` is a CLI application that makes it easy to match variants (same strand) between a set of reference and target data .pvar/bim files. The application:

* Uses allelic frequency (afreq) and variant missingness (vmiss) to evaluate whether the variants in the TARGET are suitable for inclusion in a PCA analysis
* Filters matches on strand ambiguity and multi-allelic/INDEL status
* Outputs a tab delimited report describing each intersected variant

Installation
-------------

::

    $ pip install pgscatalog-match

Usage
-----

Matching HGDP vs 1000 Genomes 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To intersect variants in HGDP against 1000 Genomes:

::

    $ pgscatalog-intersect --ref GRCh38_1000G_ALL.pvar.zst \
        --target GRCh38_hgdp_5.pvar.zst \
        --chrom 5 \
        --maf_target 0.1 \
        --geno_miss 0.1 \
        --outdir . \
        -v

You'll also need gzipped `afreq <https://www.cog-genomics.org/plink/2.0/formats#afreq>`_ and `vmiss <https://www.cog-genomics.org/plink/2.0/formats#vmiss>`_ files for the target genome in the same directory. 

The output is a tab delimited text file structured to contain one variant per row, with the following columns:

.. list-table:: 
   :widths: 25 75 
   :header-rows: 1

   * - Column name 
     - Description
   * - CHR:POS:A0:A1
     - A colon delimited variant ID, consisting of chromosome, position, a0 and a1 alleles (like REF / ALT)
   * - ID_REF 
     - The ID of the reference variant
   * - REF_REF 
     - The REF allele of the reference variant
   * - IS_INDEL
     - Is the reference variant an indel?
   * - STRANDAMB
     - Is the reference variant strand ambiguous?
   * - IS_MA_REF
     - Is the reference variant multi-allelic?
   * - ID_TARGET
     - The ID of the matched variant in the target genome
   * - REF_TARGET
     - The REF allele of the target variant
   * - IS_MA_TARGET
     - Is the target variant multi-allelic?
   * - AAF
     - Allele frequency
   * - F_MISS_DOSAGE
     - Missing dosage rate
   * - SAME_REF
     - Do the reference and target variant have the same reference allele?
   * - PCA_ELIGIBLE
     - Is the variant eligible for PCA inclusion?

You can use the table to extract IDs and then use something like `plink2 <https://www.cog-genomics.org/plink/2.0/>`_ to extract a subset of variants. 

Help
----

::

    $ pgscatalog-intersect --help