How to adjust PGS in the context of genetic ancestry
====================================================

``pgscatalog-ancestry-ajust`` is a CLI application that `adjusts calculated PGS using genetic ancestry similarity data`_.

.. _`adjusts calculated PGS using genetic ancestry similarity data`: https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html 

The process requires some complex inputs. If you'd like to automatically adjust PGS in the context of genetic ancestry automatically, it's probably simplest to use the `PGS Catalog Calculator`_, which uses this application internally.

.. _`PGS Catalog Calculator`: https://github.com/PGScatalog/pgsc_calc

Installation
-------------

::

    $ pip install pgscatalog-calcapp

Usage
-----

The following inputs are needed:

* PGS calculated on a target and reference population in long (melted) format (``--agg_scores``)
    * PGS must be calculated using the same set of variants (the intersection)
    * Input file structure is the same as the output from ``pgscatalog-aggregate``
* PCA derived from a reference population using `fraposa-pgsc`_ (``--ref_pcs``)
* PCA for target population projected onto reference PCA using ``fraposa-pgsc`` (``--target_pcs``)
* A relatedness cutoff file for the reference population (``--reference_related``)
* A sample information file for the reference population containing genetic ancestry group labels (``--psam``, ``SuperPop`` column by default)
* Text labels of the reference and target samples (``--reference`` and ``--dataset``), which must match the aggregated scoring file contents

.. _`fraposa-pgsc`: https://pypi.org/project/fraposa-pgsc/

::

    $ pgscatalog-ancestry-adjust --agg_scores aggregated_scores.txt.gz \
        --psam reference.psam \
        --ref_pcs ref.pcs \
        --target_pcs target.pcs \
        --reference_related ref.king.cutoff.id \
        --dataset hgdp \
        --reference reference \
        --outdir ancestry_results 

Help
----

::

    $ pgscatalog-ancestry-adjust --help