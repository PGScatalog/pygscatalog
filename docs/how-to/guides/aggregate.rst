How to aggregate PGS split across multiple files
=================================================

``pgscatalog-aggregate`` is a CLI application that makes it easy to aggregates calculated PGS that have been split across multiple files (e.g. by chromosome or effect type). 

If you'd like to automatically aggregate PGS for complicated data, it's probably simplest to use the `PGS Catalog Calculator`_, which uses this application internally.

.. _`PGS Catalog Calculator`: https://github.com/PGScatalog/pgsc_calc

Installation
-------------

::

    $ pip install pgscatalog-calc

Usage
-----

Input data (``-s``) are expected to be in ``plink2 --score`` output format.

::

    $ pgscatalog-ancestry-adjust -s hgdp_14_additive_0.sscore.zst hgdp_17_additive_0.sscore.zst ... \
        -o aggregated \
        --verbose

The application supports input data that are optionally compressed with gzip or zstandard. 

Help
----

::

    $ pgscatalog-ancestry-adjust --help