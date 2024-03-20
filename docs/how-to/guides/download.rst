How to download scoring files from the PGS Catalog
==================================================

``pgscatalog-download`` is a CLI application that makes it easy to download scoring files from the
PGS Catalog with a mixture of PGS, publication, or trait accessions. The application:

* automatically retries downloads if they fail
* validates the checksum of downloaded scoring files
* automatically selects scoring files aligned to a requested genome build

Installation
-------------

::

    $ pip install pgscatalog-core

Usage
-----

Downloading PGS IDs scoring files aligned to GRCh38
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    $ mkdir downloads
    $ pgscatalog-download --pgs PGS000822 PGS001229 --build GRCh38 -o downloads

.. note::

    Setting ``--build`` will download scoring files harmonised by the PGS Catalog. This means scoring fields have consistent fields, like genomic coordinates.

Downloading all scores associated with a trait
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To download all scores associated with Alzheimer's disease:

::

    $ mkdir downloads
    $ pgscatalog-download --efo MONDO_0004975 -b GRCh38 -o downloads

By default scores associated with child traits, like late-onset Alzheimer's disease, are included.
To exclude them use:

::

    $ mkdir downloads
    $ pgscatalog-download --efo MONDO_0004975 -b GRCh38 -o downloads --efo_direct

Downloading all scores associated with a publication
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're interested in scores from a specific publication:

::

    $ mkdir downloads
    $ pgscatalog-download --pgp PGP000517 -b GRCh38 -o downloads

Help
----

::

    $ pgscatalog-download --help