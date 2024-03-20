How to match scoring file variants against target genomes
=========================================================

``pgscatalog-match`` is a CLI application that makes it easy to match genetic variants in a normalised scoring file against target variant information files. The application:

* identifies match candidates using genomic coordinates and allele information
* creates a summary log and calculates match rates using the best match candidates
* writes scoring files in plink2 --score format using the best match candidates
* creates a log that describes all possible match candidate and match status

The application will error if not enough variants in the scoring file are present in the target variant information files. This is because it's important to match variants well to faithfully reproduce published scoring files.


Installation
-------------

::

    $ pip install pgscatalog-match

Usage
-----

Match variants
~~~~~~~~~~~~~~

::

    $ mkdir matchout
    $ pgscatalog-match --dataset test --scorefiles normalised_scorefile.txt --target variants.pvar --outdir matchout --min_overlap 0.75
    
.. note::

    Variant information files in both plink2 pvar and plink1 bim format are supported

Matching very large datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're matching many millions of variants it can be a good idea to split work by matching each chromosome separately:

::

    $ mkdir matchout
    $ pgscatalog-match --dataset test --scorefiles normalised_scorefile.txt --target variants_chrom1.pvar --outdir matchout --chrom 1 --only_match

This will generate match files in Arrow IPC format in matchout/matchtmp. This process can be run in parallel to distribute work.

You can then merge these matches with a different CLI tool:

::
   
   $ pgscatalog-matchmerge --dataset test --scorefiles normalised_scorefile.txt --matches path/to/matchfile --outdir matchmergeout --min_overlap 0.75

The PGS Catalog Calculator does this automatically when working with target genomes split to have one file per chromosome.   
    
Help
----

::

    $ pgscatalog-match --help
