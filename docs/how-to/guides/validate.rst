How to validate PGS Catalog scoring files
==================================================

``pgscatalog-validate`` is a CLI application that validates a set of scoring files
to match the PGS Catalog scoring file formats.
It can validate:

* The formatted scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring)
* The harmonized (Position) scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring_hm_pos)

Installation
-------------

::

    $ pip install pgscatalog-validate

Usage
-----

Validating one scoring file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    $ mkdir log
    $ pgscatalog-validate -t formatted -f PGS000001.txt.gz --log_dir log/


Validating all the scoring files in a directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    $ mkdir log
    $ pgscatalog-validate -t formatted --dir scores_directory/ --log_dir log/


Help
----

::

    $ pgscatalog-validate --help