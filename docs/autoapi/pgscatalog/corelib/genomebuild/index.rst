:orphan:

:py:mod:`pgscatalog.corelib.genomebuild`
========================================

.. py:module:: pgscatalog.corelib.genomebuild


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.genomebuild.GenomeBuild




.. py:class:: GenomeBuild


   Bases: :py:obj:`enum.Enum`

   Enumeration of genome build: the reference genome release that a scoring file
   is aligned to.

   >>> GenomeBuild.GRCh38
   GenomeBuild.GRCh38

   .. py:attribute:: GRCh37
      :value: 'GRCh37'

      

   .. py:attribute:: GRCh38
      :value: 'GRCh38'

      

   .. py:attribute:: NCBI36
      :value: 'NCBI36'

      

   .. py:method:: from_string(build)
      :classmethod:

      :param build: genome build string
      :return: :class:`GenomeBuild`
      :raises ValueError: From an unsupported build string

      >>> GenomeBuild.from_string("GRCh38")
      GenomeBuild.GRCh38
      >>> str(GenomeBuild.from_string("GRCh37"))
      'GRCh37'
      >>> GenomeBuild.from_string("NR") is None
      True
      >>> GenomeBuild.from_string("pangenome")
      Traceback (most recent call last):
      ...
      ValueError: Can't match build='pangenome'



