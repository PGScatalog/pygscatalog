.. pygscatalog documentation master file, created by
   sphinx-quickstart on Tue Jan 23 13:49:00 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pygscatalog
===========

.. rtfd-shield::
    :project: pygscatalog

.. github-shield::
    :last-commit:
    :branch: main
    
.. actions-shield::
    :workflow: Run pytest

.. pre-commit-shield::



``pygscatalog`` provides a set of Python CLI applications and developer libraries for working with polygenic scores (PGS),
including integration with the PGS Catalog.

These applications and libraries are used internally by the `PGS Catalog Calculator`_, which is an automated workflow for calculating PGS, including adjustment of scores in the context of `genetic ancestry similarity`_.

If you're interested in PGS but aren't sure where to begin, the calculator is the best place.

If you're working with PGS data and want to do some kinds of bespoke analysis not supported by the calculator, these tools might be helpful.

.. _`PGS Catalog Calculator`: https://github.com/PGScatalog/pgsc_calc
.. _`genetic ancestry similarity`: https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   how-to/index

.. sidebar-links::
   :github: