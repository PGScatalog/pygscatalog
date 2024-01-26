:orphan:

:py:mod:`pgscatalog.corelib._read`
==================================

.. py:module:: pgscatalog.corelib._read

.. autoapi-nested-parse::

   This module contains functions for reading data from PGS Catalog files.
   These functions aren't really meant to be imported outside corelib 



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib._read.auto_open
   pgscatalog.corelib._read.detect_wide
   pgscatalog.corelib._read.generate_header_lines
   pgscatalog.corelib._read.get_columns
   pgscatalog.corelib._read.read_header
   pgscatalog.corelib._read.read_rows_lazy



Attributes
~~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib._read.logger


.. py:function:: auto_open(filepath, mode='rt')

   Automatically open a gzipped text file or an uncompressed text file


.. py:function:: detect_wide(cols: list[str]) -> bool

   Check columns to see if multiple effect weights are present. Multiple effect weights must be present in the form:
   effect_weight_suffix1
   effect_weight_suffix2


.. py:function:: generate_header_lines(f)

   Header lines in a PGS Catalog scoring file are structured like:
       #pgs_id=PGS000348
       #pgs_name=PRS_PrCa
   Files can be big, so we want to only read header lines and stop immediately


.. py:function:: get_columns(path)

   Grab column labels from a PGS Catalog scoring file. line_no is useful to skip the header


.. py:function:: read_header(path: pathlib.Path)

   Parses the header of a PGS Catalog format scoring file into a dictionary


.. py:function:: read_rows_lazy(*, csv_reader, fields: list[str], name: str, wide: bool, row_nr: int)

   Read rows from an open scoring file and instantiate them as ScoreVariants


.. py:data:: logger

   

