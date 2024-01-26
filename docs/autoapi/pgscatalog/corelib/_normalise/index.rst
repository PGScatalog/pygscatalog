:orphan:

:py:mod:`pgscatalog.corelib._normalise`
=======================================

.. py:module:: pgscatalog.corelib._normalise

.. autoapi-nested-parse::

   This module contains a data processing pipeline to format ScoreVariants in a
   standard way. Each step in the data processing pipeline is a generator that operates
   on a list of ScoreVariants and yields updated ScoreVariants. This makes it easy to
   plug in extra steps where needed, and lazily works on millions of objects.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib._normalise.assign_effect_type
   pgscatalog.corelib._normalise.assign_other_allele
   pgscatalog.corelib._normalise.check_bad_variant
   pgscatalog.corelib._normalise.check_duplicates
   pgscatalog.corelib._normalise.check_effect_allele
   pgscatalog.corelib._normalise.check_effect_weight
   pgscatalog.corelib._normalise.detect_complex
   pgscatalog.corelib._normalise.drop_hla
   pgscatalog.corelib._normalise.lift
   pgscatalog.corelib._normalise.load_chain
   pgscatalog.corelib._normalise.normalise
   pgscatalog.corelib._normalise.remap_harmonised



Attributes
~~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib._normalise.logger


.. py:function:: assign_effect_type(variants)

   Convert PGS Catalog effect type columns to EffectType enums

   The most common type of effect type is additive:
   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "is_recessive": "FALSE", "is_dominant": "FALSE"})
   >>> list(assign_effect_type([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(...,effect_type=EffectType.ADDITIVE,...)]

   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "is_recessive": "TRUE", "is_dominant": "FALSE"})
   >>> list(assign_effect_type([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(...,effect_type=EffectType.RECESSIVE,...)]

   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "is_recessive": "FALSE", "is_dominant": "TRUE"})
   >>> list(assign_effect_type([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(...,effect_type=EffectType.DOMINANT,...)]


.. py:function:: assign_other_allele(variants)

   Check if there's more than one possible other allele, remove if true

   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "other_allele": "A"})
   >>> list(assign_other_allele([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(effect_allele='A',...,other_allele='A',...)]
   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "other_allele": "A/C"})
   >>> list(assign_other_allele([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(effect_allele='A',...,other_allele=None,...)]


.. py:function:: check_bad_variant(variants, drop_missing=False)

   Missing effect allele:
   >>> variant = ScoreVariant(**{"effect_allele": None, "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(check_bad_variant([variant], drop_missing=True)) # doctest: +ELLIPSIS
   []

   Missing chromosome name and position:
   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(check_bad_variant([variant], drop_missing=True)) # doctest: +ELLIPSIS
   []


.. py:function:: check_duplicates(variants)

   Check if a scoring file contains multiple variants with the same ID
   ID = chr:pos:effect_allele:other_allele


.. py:function:: check_effect_allele(variants, drop_missing=False)

   Odd effect allele:
   >>> variant = ScoreVariant(**{"effect_allele": "Z", "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(check_effect_allele([variant], drop_missing=True)) # doctest: +ELLIPSIS
   []


   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(check_effect_allele([variant], drop_missing=True)) # doctest: +ELLIPSIS
   [ScoreVariant(effect_allele='A'...)]


.. py:function:: check_effect_weight(variants)

   Check that effect weights are valid floats. Effect weights are intentionally
   left as strings during processing.

   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(check_effect_weight([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(effect_allele='A',effect_weight=5,...)]

   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": "potato", "accession": "test", "row_nr": 0})
   >>> list(check_effect_weight([variant])) # doctest: +ELLIPSIS
   Traceback (most recent call last):
   ...
   ValueError


.. py:function:: detect_complex(variants)

   Some older scoring files in the PGS Catalog are complicated.
   They often require bespoke set up to support interaction terms, etc
   This function only exists to provide loud warnings to end users.


.. py:function:: drop_hla(variants)

   Drop HLA alleles from a list of ScoreVariants

   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(drop_hla([variant])) # doctest: +ELLIPSIS
   [ScoreVariant(effect_allele='A',...)]

   >>> variant = ScoreVariant(**{"effect_allele": "P", "effect_weight": 5, "accession": "test", "row_nr": 0})
   >>> list(drop_hla([variant]))
   []


.. py:function:: lift(*, scoring_file, harmonised, current_build, target_build, chain_dir, min_lift=0.95)


.. py:function:: load_chain(*, current_build, target_build, chain_dir)

   Only supports loading GRCh37 and GRCh38 chain files

   >>> from ._config import Config
   >>> chain_dir = Config.ROOT_DIR / "tests" / "chain"
   >>> load_chain(current_build=GenomeBuild.GRCh37, target_build=GenomeBuild.GRCh38, chain_dir=chain_dir) # doctest: +ELLIPSIS
   <pyliftover.liftover.LiftOver object at...

   >>> load_chain(current_build=GenomeBuild.GRCh38, target_build=GenomeBuild.GRCh37, chain_dir=chain_dir) # doctest: +ELLIPSIS
   <pyliftover.liftover.LiftOver object at...

   >>> load_chain(current_build=GenomeBuild.NCBI36, target_build=GenomeBuild.GRCh38, chain_dir=chain_dir)
   Traceback (most recent call last):
   ...
   ValueError: Unsupported liftover current_build=GenomeBuild.NCBI36, target_build=GenomeBuild.GRCh38


.. py:function:: normalise(scoring_file, drop_missing=False, liftover=False, chain_dir=None, target_build=None)

   Order of steps is important:

   1. liftover non-harmonised data (quite rare), failed lifts get None'd
   2. remap harmonised data, failed harmonisations get None'd
   3. log and optionally drop bad variants


.. py:function:: remap_harmonised(variants, harmonised)

   Overwrite key attributes with harmonised data, if available.

   In this case chr_name, chr_position, and other allele are missing.
   Perhaps authors submitted rsID and effect allele originally:
   >>> variant = ScoreVariant(**{"effect_allele": "A", "effect_weight": 5, "accession": "test", "row_nr": 0, "hm_chr": 1, "hm_pos": 100, "hm_inferOtherAllele": "A"})
   >>> list(remap_harmonised([variant], harmonised=True)) # doctest: +ELLIPSIS
   [ScoreVariant(...,chr_name=1,chr_position=100,...other_allele='A'...)]


.. py:data:: logger

   

