How to format scoring files from the PGS Catalog
=================================================

``pgscatalog-format`` is a CLI application that makes it easy to combine scoring files into a standardised output.

.. note:: ``pgscatalog-combine`` was recently renamed to ``pgscatalog-format``

The process involves:

* extracting important fields from scoring files
* doing some quality control checks
* optionally lifting over variants to a consistent genome build
* writing to a consistent schema

Input scoring files must follow PGS Catalog standards. The output file is useful for
doing data science tasks, like matching variants across a scoring file and target
genome.

Installation
------------

::

    $ pipx install pgscatalog-core

Usage
-----

Combining PGS Catalog scoring files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tip:: It's easiest to get started by downloading scoring files in the same genome build: :doc:`download`

::

    $ mkdir output
    $ pgscatalog-combine -s PGS000001_hmPOS_GRCh38.txt.gz PGS0001229_hmPOS_GRCh38.txt.gz -t GRCh38 -o output

For each input scoring file a formatted scoring file will be written to the output directory.

.. tip:: If you're formatting lots of scoring files in parallel the ``--threads`` parameter can help speed up the process

Lifting over scoring files
~~~~~~~~~~~~~~~~~~~~~~~~~~

It's possible to combine scoring files with different genome builds using liftover.

.. danger:: You should only do this when combining PGS Catalog and custom scoring files, because the PGS Catalog provides harmonised data

First, download chain files from UCSC:

* `hg19ToHg38.over.chain.gz`_
* `hg38ToHg19.over.chain.gz`_

.. _hg19ToHg38.over.chain.gz: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
.. _hg38ToHg19.over.chain.gz: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/

And copy them into a directory (e.g. ``my_chain_dir/``).

Assuming you have a custom scoring file in GRCh37 (``my_scorefile_grch37.txt.gz``), and you want to combine it with a PGS Catalog scoring file in GRCh38.

::

    $ mkdir output
    $ pgscatalog-format -s PGS000001_hmPOS_GRCh38.txt.gz my_scorefile_grch37.txt.gz \
        --chain_dir my_chain_dir/ \
        -t GRCh38 \
        -o output

Help
----

::

    $ pgscatalog-format --help
