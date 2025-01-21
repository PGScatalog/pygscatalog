How to validate PGS Catalog scoring files
==================================================

``pgscatalog-validate`` is a CLI application that validates a set of scoring files
to match the PGS Catalog scoring file formats.
It can validate:

* The formatted scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring)
* The harmonized (Position) scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring_hm_pos)

Installation
-------------

Download the source code:

::
    $ git clone https://github.com/PGScatalog/pygscatalog.git
    $ cd pygscatalog/pygscatalog.validate

Install the dependencies:

::
    $ poetry install

Run with a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start a new local virtual envrironment:

::
    $ poetry shell
    $ pgscatalog-validate --help

Or install the package in the current environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Build the wheel package:

::
    $ poetry build

Install the built package in the current environment via pip:

::
    $ pip install dist/pgscatalog_validate-0.1-py3-none-any.whl

Run

::
    pgscatalog-validate --help

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