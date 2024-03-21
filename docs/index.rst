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

.. sidebar-links::
   :github:

``pygscatalog`` provides a set of Python CLI applications and developer libraries for working with polygenic scores (PGS),
including integration with the `PGS Catalog`_.

These applications and libraries are used internally by the `PGS Catalog Calculator`_, which is an automated workflow for calculating PGS, including adjustment of scores in the context of `genetic ancestry similarity`_.

If you're interested in PGS but aren't sure where to begin, the calculator is the best place.

If you're working with PGS data and want to do some kinds of bespoke analysis not supported by the calculator, these tools might be helpful.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   how-to/index

Credits
-------

``pygscatalog`` (*aka* ``pgscatalog_utils``) is developed as part of the `PGS Catalog`_ project, a
collaboration between the University of Cambridge's Department of Public Health and Primary Care (Michael Inouye, 
Samuel Lambert) and the European Bioinformatics Institute (Helen Parkinson, Laura Harris).

This package contains code libraries and apps for working with PGS Catalog data and calculating PGS within the 
`PGS Catalog Calculator`_ (``pgsc_calc``) workflow, and is based on an earlier 
codebase (`pgscatalog_utils`_) with contributions and input from members 
of the PGS Catalog team (Samuel Lambert, Benjamin Wingfield, Aoife McMahon Laurent Gil) and Inouye lab 
(Rodrigo Canovas, Scott Ritchie, Jingqin Wu).

A manuscript describing this package and ``pgsc_calc`` pipeline is *in preparation*. In the meantime if you use the
tool we ask you to cite the repo and the paper describing the PGS Catalog resource:

- `PGS Catalog Calculator`_ *in preparation*. PGS Catalog Team. https://github.com/PGScatalog/pgsc_calc
- Lambert *et al.* (2021) The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation.  Nature Genetics. 53:420-425 doi: `10.1038/s41588-021-00783-5`_

All of our code is open source and permissively licensed with `Apache 2`_.

This work has received funding from EMBL-EBI core funds, the Baker Institute, the University of Cambridge, 
Health Data Research UK (HDRUK), and the European  Union’s Horizon 2020 research and innovation programme under grant 
agreement No 101016775 INTERVENE.

.. _`PGS Catalog`: https://www.pgscatalog.org/
.. _`PGS Catalog Calculator`: https://github.com/PGScatalog/pgsc_calc
.. _`genetic ancestry similarity`: https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html
.. _`pgscatalog_utils`: https://github.com/PGScatalog/pgscatalog_utils
.. _`10.1038/s41588-021-00783-5`: https://doi.org/10.1038/s41588-021-00783-5
.. _`Apache 2`: https://github.com/PGScatalog/pygscatalog/blob/main/README.md

