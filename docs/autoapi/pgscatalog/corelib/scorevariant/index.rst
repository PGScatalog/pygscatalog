:orphan:

:py:mod:`pgscatalog.corelib.scorevariant`
=========================================

.. py:module:: pgscatalog.corelib.scorevariant

.. autoapi-nested-parse::

   This module contains classes that compose a ScoreVariant: a variant in a PGS
   Catalog Scoring File.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.scorevariant.EffectAllele
   pgscatalog.corelib.scorevariant.EffectType
   pgscatalog.corelib.scorevariant.ScoreVariant




.. py:class:: EffectAllele(allele)


   A class that represents an effect allele found in PGS Catalog scoring files

   The allele that's dosage is counted (e.g. {0, 1, 2}) and multiplied by the variant's
   weight (effect_weight) when calculating score. The effect allele is also known as
   the 'risk allele'.

   >>> simple_ea = EffectAllele("A")
   >>> simple_ea
   EffectAllele("A")
   >>> simple_ea.is_snp
   True
   >>> str(simple_ea)
   'A'
   >>> EffectAllele("AG")
   EffectAllele("AG")
   >>> hla_example = EffectAllele("+")
   >>> hla_example
   EffectAllele("+")
   >>> hla_example.is_snp
   False

   .. py:property:: allele


   .. py:property:: is_snp
      :type: bool

      SNPs are the most common type of effect allele in PGS Catalog scoring
      files. More complex effect alleles, like HLAs or APOE genes, often require
      extra work to represent in genomes. Users should be warned about complex
      effect alleles.

      >>> ea = EffectAllele("+")
      >>> ea.is_snp
      False
      >>> ea.allele = "A"
      >>> ea.is_snp
      True



.. py:class:: EffectType


   Bases: :py:obj:`enum.Enum`

   Enums that represent inheritance models. The vast majority of variants have
   an additive effect type.

   This changes downstream PGS calculation:

   * ScoreVariants with an additive effect type will always be added to the PGS sum.
   * ScoreVariants with a dominant effect type are only added to the PGS sum if there is at least one copy of the effect allele.
   * ScoreVariants with a recessive effect type are only added to the PGS sum if there are two copies of the effect allele.

   >>> EffectType.ADDITIVE
   EffectType.ADDITIVE
   >>> str(EffectType.ADDITIVE)
   'additive'

   .. py:attribute:: ADDITIVE
      :value: 'additive'

      

   .. py:attribute:: DOMINANT
      :value: 'dominant'

      

   .. py:attribute:: RECESSIVE
      :value: 'recessive'

      


.. py:class:: ScoreVariant(*, effect_allele: str, effect_weight: str, accession: str, row_nr: int, chr_name: str = None, chr_position: int = None, rsID: str = None, other_allele: str = None, hm_chr: str = None, hm_pos: int = None, hm_inferOtherAllele: str = None, hm_source: str = None, is_dominant: str = None, is_recessive: str = None, hm_rsID: str = None, hm_match_chr: str = None, hm_match_pos: str = None, is_duplicated: bool = False, effect_type: EffectType = EffectType.ADDITIVE, is_complex: bool = False, **kwargs)


   Represents a single row in a PGS Catalog scoring file.

   It's rare to instantiate this class directly. Instead, create a
   class:`ScoringFile`  from a path and you can lazily iterate over variants.

   .. py:attribute:: complex_fields
      :type: tuple[str]
      :value: ('is_haplotype', 'is_diplotype', 'is_interaction')

      

   .. py:attribute:: mandatory_fields
      :type: tuple[str]
      :value: ('effect_allele', 'effect_weight', 'accession', 'row_nr')

      

   .. py:attribute:: optional_fields
      :type: tuple[str]
      :value: ('chr_name', 'chr_position', 'rsID', 'other_allele', 'hm_chr', 'hm_pos', 'hm_inferOtherAllele',...

      

   .. py:attribute:: output_fields
      :type: tuple[str]
      :value: ('chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',...

      


