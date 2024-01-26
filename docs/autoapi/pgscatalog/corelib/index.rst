:py:mod:`pgscatalog.corelib`
============================

.. py:module:: pgscatalog.corelib


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.CatalogCategory
   pgscatalog.corelib.CatalogQuery
   pgscatalog.corelib.Config
   pgscatalog.corelib.GenomeBuild
   pgscatalog.corelib.ScoreQueryResult
   pgscatalog.corelib.ScoreVariant
   pgscatalog.corelib.ScoringFile
   pgscatalog.corelib.ScoringFiles




.. py:exception:: BasePGSException(message)


   Bases: :py:obj:`Exception`

   The base class from which all PGS errors must inherit.
   The purpose of this class is to simplify finding PGS exceptions and exiting python
   with a matching custom exit code.


.. py:exception:: BuildError(message)


   Bases: :py:obj:`CombineError`

   Raised when there's a problem with a scoring file genome build.


.. py:exception:: CatalogError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors when querying or downloading from the PGS Catalog


.. py:exception:: CombineError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors that are raised when combining scorefiles


.. py:exception:: DuplicateMatchError(message)


   Bases: :py:obj:`MatchError`

   Raised when a matched variant has been duplicated, so that a variant with the same ID
   would be split across two rows in an output scoring file.


.. py:exception:: GenomesNotFound(message)


   Bases: :py:obj:`SamplesheetError`

   Raised when FileNotFound


.. py:exception:: InvalidAccessionError(message)


   Bases: :py:obj:`CatalogError`

   Raised when an invalid term is used to query the Catalog


.. py:exception:: MatchError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors that are raised during variant matching


.. py:exception:: MatchRateError(message)


   Bases: :py:obj:`MatchError`

   Raised when match rate is below match threshold for one or more scoring files


.. py:exception:: MatchValueError(message)


   Bases: :py:obj:`MatchError`

   Raised when a match function receives inappropriate values.

   e.g., Multiple chromosomes detected in variant data but data is split per-chromosome


.. py:exception:: QueryError(message)


   Bases: :py:obj:`CatalogError`

   Raised when the Catalog API doesn't return a valid response


.. py:exception:: SamplesheetError(message)


   Bases: :py:obj:`BasePGSException`

   The base class for errors related to samplesheet parsing


.. py:exception:: SamplesheetFormatError(message)


   Bases: :py:obj:`SamplesheetError`

   Raised when a samplesheet is badly formatted


.. py:exception:: ScoreChecksumError(message)


   Bases: :py:obj:`CatalogError`

   Raised when a scoring file fails checksum validation


.. py:exception:: ScoreDownloadError(message)


   Bases: :py:obj:`CatalogError`

   Raised when a scoring file can't be downloaded


.. py:exception:: ScoreFormatError(message)


   Bases: :py:obj:`CombineError`

   Raised when there's a problem with a scoring file.


.. py:exception:: ZeroMatchesError(message)


   Bases: :py:obj:`MatchError`

   Raised when zero matches are found for one or more scoring files.

   Distinct from MatchRateError because it's very common, and caused by bad input data or parameters.


.. py:class:: CatalogCategory


   Bases: :py:obj:`enum.Enum`

   The three main categories in the PGS Catalog

   Enumeration values don't mean anything and are automatically generated:

   >>> CatalogCategory.SCORE
   <CatalogCategory.SCORE: 1>

   .. py:attribute:: PUBLICATION

      

   .. py:attribute:: SCORE

      

   .. py:attribute:: TRAIT

      


.. py:class:: CatalogQuery(*, accession, include_children=False, **kwargs)


   Efficiently query the PGS Catalog API using accessions

   Supports trait (EFO), score (PGS ID), or publication identifier (PGP ID)

   >>> CatalogQuery(accession="PGS000001")
   CatalogQuery(accession='PGS000001', category=CatalogCategory.SCORE, include_children=None)

   Supports multiple PGS ID input in a list:

   >>> CatalogQuery(accession=["PGS000001", "PGS000002"])
   CatalogQuery(accession=['PGS000001', 'PGS000002'], category=CatalogCategory.SCORE, include_children=None)

   Duplicates are automatically dropped:

   >>> CatalogQuery(accession=["PGS000001", "PGS000001"])
   CatalogQuery(accession=['PGS000001'], category=CatalogCategory.SCORE, include_children=None)

   Publications and trait accessions are supported too:

   >>> CatalogQuery(accession="PGP000001")
   CatalogQuery(accession='PGP000001', category=CatalogCategory.PUBLICATION, include_children=None)

   >>> CatalogQuery(accession="EFO_0001645")
   CatalogQuery(accession='EFO_0001645', category=CatalogCategory.TRAIT, include_children=False)

   .. py:method:: get_query_url()

      Automatically resolve a query URL for a PGS Catalog accession (or multiple
      score accessions).

      A list is returned because when querying multiple score accessions batches
      are created:

      >>> CatalogQuery(accession=["PGS000001","PGS000002"]).get_query_url()
      ['https://www.pgscatalog.org/rest/score/search?pgs_ids=PGS000001,PGS000002']

      (each element in this list contains up to 50 score IDs)

      Multiple score accessions are automatically deduplicated:

      >>> CatalogQuery(accession = ["PGS000001"] * 100).get_query_url()
      ['https://www.pgscatalog.org/rest/score/search?pgs_ids=PGS000001']

      Publications don't batch because they natively support many scores:

      >>> CatalogQuery(accession="PGP000001").get_query_url()
      'https://www.pgscatalog.org/rest/publication/PGP000001'

      Traits don't batch for the same reason as publications:

      >>> CatalogQuery(accession="EFO_0001645").get_query_url()
      'https://www.pgscatalog.org/rest/trait/EFO_0001645?include_children=0'

      Child traits terms aren't included by default. Only traits can have children.


   .. py:method:: infer_category()

      Inspect an accession and guess the Catalog category

      >>> CatalogQuery(accession="PGS000001").infer_category()
      <CatalogCategory.SCORE: 1>

      >>> CatalogQuery(accession="EFO_0004346").infer_category()
      <CatalogCategory.TRAIT: 2>

      >>> CatalogQuery(accession="MONDO_0005041").infer_category()
      <CatalogCategory.TRAIT: 2>

      >>> CatalogQuery(accession="PGP000001").infer_category()
      <CatalogCategory.PUBLICATION: 3>

      Be careful, assume lists of accessions only contain PGS IDs:

      >>> CatalogQuery(accession=["PGS000001", "PGS000002"]).infer_category()
      <CatalogCategory.SCORE: 1>


   .. py:method:: score_query()

      Query the PGS Catalog API and return :class:`ScoreQueryResult`

      Information about a single score is returned as a dict:

      >>> CatalogQuery(accession="PGS000001").score_query() # doctest: +ELLIPSIS
      ScoreQueryResult(pgs_id='PGS000001', ftp_url=...

      If information about multiple scores is found, it's returned as a list:

      >>> CatalogQuery(accession=["PGS000001", "PGS000002"]).score_query() # doctest: +ELLIPSIS
      [ScoreQueryResult(pgs_id='PGS000001', ftp_url=...

      Publications and traits always return a list of score information:

      >>> CatalogQuery(accession="PGP000001").score_query() # doctest: +ELLIPSIS
      [ScoreQueryResult(pgs_id='PGS000001', ftp_url=...



.. py:class:: Config


   This class stores global package configuration as class attributes

   Most of the time you won't need to change anything.

   .. py:attribute:: API_HEADER

      

   .. py:attribute:: BATCH_SIZE
      :value: 20000

      

   .. py:attribute:: FTP_EXCLUSIVE
      :value: False

      

   .. py:attribute:: MAX_RETRIES
      :value: 5

      

   .. py:attribute:: ROOT_DIR

      


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



.. py:class:: ScoreQueryResult(*, pgs_id, ftp_url, ftp_grch37_url, ftp_grch38_url, license)


   Class that holds score metadata with methods to extract important fields

   .. py:method:: from_query(result_response)
      :classmethod:

      Parses PGS Catalog API JSON response

      :param result_response: PGS Catalog API JSON response
      :returns: :class:`ScoreQueryResult`

      >>> fake_response = {"id": "fake", "ftp_harmonized_scoring_files":
      ... {"GRCh37": {"positions": "fake.txt.gz"}, "GRCh38": {"positions": "fake.txt.gz"}},
      ... "license": "fake", "ftp_scoring_file": "fake.txt.gz"}
      >>> ScoreQueryResult.from_query(fake_response) # doctest: +ELLIPSIS
      ScoreQueryResult(pgs_id='fake', ftp_url='fake.txt.gz',...


   .. py:method:: get_download_url(genome_build=None)

      Returns scoring file download URL, with support for specifying harmonised data in a specific genome build

      >>> query = CatalogQuery(accession="PGS000001").score_query()
      >>> build = GenomeBuild.GRCh38
      >>> query.get_download_url()
      'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/PGS000001.txt.gz'
      >>> query.get_download_url(build)
      'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/Harmonized/PGS000001_hmPOS_GRCh38.txt.gz'



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

      


.. py:class:: ScoringFile(identifier, target_build=None, query_result=None, **kwargs)


   Represents a single scoring file in the PGS Catalog.

   :param identifier: A PGS Catalog score accession in the format ``PGS123456`` or a path to a local scoring file
   :param target_build: An optional :class:`GenomeBuild`, which represents the build you want the scoring file to align to
   :param query_result: An optional :class:`ScoreQueryResult`, if provided with an accession identifier it prevents hitting the PGS Catalog API
   :raises InvalidAccessionError: If the PGS Catalog API can't find the provided accession
   :raises ScoreFormatError: If you try to iterate over a ``ScoringFile`` without a local path (before downloading it)

   You can make ``ScoringFiles`` with a path to a scoring file:

   >>> sf = ScoringFile(Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt.gz")
   >>> sf # doctest: +ELLIPSIS
   ScoringFile('.../PGS000001_hmPOS_GRCh38.txt.gz', target_build=None)

   >>> sf.genome_build
   GenomeBuild.GRCh38

   >>> sf.pgs_id
   'PGS000001'

   >>> for variant in sf.variants: # doctest: +ELLIPSIS
   ...     variant
   ...     break
   ScoreVariant(effect_allele='T',effect_weight='0.16220387987485377',accession='PGS000001',...

   You can also make a ``ScoringFile`` by using PGS Catalog score accessions:

   >>> sf = ScoringFile("PGS000001", target_build=GenomeBuild.GRCh38)
   >>> sf
   ScoringFile('PGS000001', target_build=GenomeBuild.GRCh38)

   It's important to use the ``.download()`` method when you're not working with local files,
   or many attributes and methods will be missing or won't work:

   >>> for variant in sf.variants:
   ...     variant
   ...     break
   Traceback (most recent call last):
   ...
   corelib.pgsexceptions.ScoreFormatError: Local file is missing. Did you .download()?

   A ``ScoringFile`` can also be constructed with a ``ScoreQueryResult`` if you want
   to be polite to the PGS Catalog API. Just add the ``query_result`` parameter:

   >>> score_query_result = sf.catalog_response  # extract score query from old query
   >>> ScoringFile(identifier=sf.pgs_id, query_result=sf.catalog_response)  # doesn't hit the PGS Catalog API again
   ScoringFile('PGS000001', target_build=None)

   :class:`InvalidAccessionError` is raised if you provide bad identifiers:

   >>> import tempfile
   >>> with tempfile.TemporaryDirectory() as tmp_dir:
   ...     ScoringFile("potato", GenomeBuild.GRCh38).download(tmp_dir)
   Traceback (most recent call last):
   ...
   corelib.pgsexceptions.InvalidAccessionError: Invalid accession: 'potato'

   The same exception is raised if you provide a well formatted identifier that doesn't exist:

   >>> with tempfile.TemporaryDirectory() as tmp_dir:
   ...     ScoringFile("PGS000000", GenomeBuild.GRCh38).download(tmp_dir)
   Traceback (most recent call last):
   ...
   corelib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGS000000'

   .. py:property:: target_build

      The ``GenomeBuild`` you want a ``ScoringFile`` to align to. Useful when using PGS
      Catalog accessions to instantiate this class.


   .. py:property:: variants

      A generator that yields rows from the scoring file as ``ScoreVariants``,
      if a local file is available (i.e. after downloading). Always available for
      class instances that have a valid local path.


   .. py:method:: download(directory, overwrite=False)

      Download a ScoringFile to a specified directory with checksum validation

      :param directory: Directory to write file to
      :param overwrite: Overwrite existing file if present

      :raises ScoreDownloadError: If there's an unrecoverable problem downloading the file
      :raises ScoreChecksumError: If md5 validation consistently fails

      :returns: None

      >>> import tempfile, os
      >>> with tempfile.TemporaryDirectory() as tmp_dir:
      ...     ScoringFile("PGS000001").download(tmp_dir)
      ...     print(os.listdir(tmp_dir))
      ['PGS000001.txt.gz']

      It's possible to request a scoring file in a specific genome build:

      >>> import tempfile, os
      >>> with tempfile.TemporaryDirectory() as tmp_dir:
      ...     ScoringFile("PGS000001", GenomeBuild.GRCh38).download(tmp_dir)
      ...     print(os.listdir(tmp_dir))
      ['PGS000001_hmPOS_GRCh38.txt.gz']



   .. py:method:: get_log(drop_missing=False, variant_log=None)

      Create a JSON log from a ScoringFile's header and variant rows. 


   .. py:method:: normalise(liftover=False, drop_missing=False, chain_dir=None, target_build=None)

      Extracts key fields from a scoring file in a normalised format.

      Takes care of quality control.

      >>> testpath = Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt.gz"
      >>> variants = ScoringFile(testpath).normalise()
      >>> for x in variants: # doctest: +ELLIPSIS
      ...     x
      ...     break
      ScoreVariant(effect_allele='T',effect_weight='0.16220387987485377',...

      Supports lifting over scoring files from GRCh37 to GRCh38:

      >>> testpath = Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh37.txt"
      >>> chaindir = Config.ROOT_DIR / "tests" / "chain"
      >>> sf = ScoringFile(testpath)
      >>> sf.harmonised = False  # lying, or liftover will be skipped
      >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh38)
      >>> for x in variants:
      ...     (x.rsID, x.chr_name, x.chr_position)
      ...     break
      ('rs78540526', '11', 69516650)

      Example of lifting down (GRCh38 to GRCh37):

      >>> testpath = Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt"
      >>> chaindir = Config.ROOT_DIR / "tests" / "chain"
      >>> sf = ScoringFile(testpath)
      >>> sf.harmonised = False  # lying, or liftover will be skipped
      >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh37)
      >>> for x in variants:
      ...     (x.rsID, x.chr_name, x.chr_position)
      ...     break
      ('rs78540526', '11', 69331418)

      Liftover support is only really useful for custom scoring files that aren't
      in the PGS Catalog. It's always best to use harmonised data when it's
      available from the PGS Catalog. Harmonised data goes through a lot of validation
      and error checking.

      For example, if you set the wrong genome build, you can get odd
      results returned without any errors, warnings, or exceptions:

      >>> testpath = Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt"
      >>> chaindir = Config.ROOT_DIR / "tests" / "chain"
      >>> sf = ScoringFile(testpath)
      >>> sf.harmonised = False  # lying, or liftover will be skipped
      >>> sf.genome_build = GenomeBuild.GRCh37  # wrong build ! it's GRCh38
      >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh38)
      >>> for x in variants:
      ...     (x.rsID, x.chr_name, x.chr_position)
      ...     break
      ('rs78540526', '11', 69701882)

      A :class:`LiftoverError` is only raised when many converted coordinates are missing.



.. py:class:: ScoringFiles(*args, target_build=None, **kwargs)


   This container class provides methods to work with multiple ScoringFile objects.

   You can use publications or trait accessions to instantiate:

   >>> ScoringFiles("PGP000001", target_build=GenomeBuild.GRCh37)
   ScoringFiles('PGS000001', 'PGS000002', 'PGS000003', target_build=GenomeBuild.GRCh37)

   Or multiple PGS IDs:

   >>> ScoringFiles("PGS000001", "PGS000002")
   ScoringFiles('PGS000001', 'PGS000002', target_build=None)

   List input is OK too:

   >>> ScoringFiles(["PGS000001", "PGS000002"])
   ScoringFiles('PGS000001', 'PGS000002', target_build=None)

   Or any mixture of publications, traits, and scores:

   >>> ScoringFiles("PGP000001", "PGS000001", "PGS000002")
   ScoringFiles('PGS000001', 'PGS000002', 'PGS000003', target_build=None)

   Scoring files with duplicate PGS IDs (accessions) are automatically dropped.
   In the example above ``PGP000001`` contains ``PGS000001``, ``PGS000002``, and ``PGS000003``.

   Traits can have children. To include these traits, use the ``include_children`` parameter:

   >>> score_with_children = ScoringFiles("MONDO_0004975", include_children=True)
   >>> score_wo_children = ScoringFiles("MONDO_0004975", include_children=False)
   >>> len(score_with_children) > len(score_wo_children)
   True

   For example, Alzheimer's disease (``MONDO_0004975``) includes Late-onset Alzheier's disease (``EFO_1001870``) as a child trait.

   Concatenation works as expected:

   >>> ScoringFiles('PGS000001') + ScoringFiles('PGS000002', 'PGS000003')
   ScoringFiles('PGS000001', 'PGS000002', 'PGS000003', target_build=None)

   But only :class:`ScoringFiles` with the same genome build can be concatenated:

   >>> ScoringFiles('PGS000001') + ScoringFiles('PGS000002', 'PGS000003', target_build=GenomeBuild.GRCh38)
   Traceback (most recent call last):
   ...
   TypeError: unsupported operand type(s) for +: 'ScoringFiles' and 'ScoringFiles'

   Multiplication doesn't make sense, because :class:`ScoringFile` elements must be unique,
   so isn't supported.

   >>> ScoringFiles('PGS000001') * 3
   Traceback (most recent call last):
   ...
   TypeError: unsupported operand type(s) for *: 'ScoringFiles' and 'int'

   You can slice and iterate over :class:`ScoringFiles`:

   >>> score = ScoringFiles("PGP000001", target_build=GenomeBuild.GRCh38)
   >>> score[0]
   ScoringFile('PGS000001', target_build=GenomeBuild.GRCh38)
   >>> for x in score:
   ...     x
   ScoringFile('PGS000001', target_build=GenomeBuild.GRCh38)
   ScoringFile('PGS000002', target_build=GenomeBuild.GRCh38)
   ScoringFile('PGS000003', target_build=GenomeBuild.GRCh38)
   >>> score[0] in score
   True

   The accession validation rules apply from :class:`ScoringFile`:

   >>> ScoringFiles("PGPpotato")
   Traceback (most recent call last):
   ...
   corelib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGPpotato'

   Local files can also be used to instantiate :class:`ScoringFiles`:

   >>> import tempfile
   >>> with tempfile.TemporaryDirectory() as d:
   ...     x = ScoringFile("PGS000001", target_build=GenomeBuild.GRCh38)
   ...     x.download(directory=d)
   ...     ScoringFiles(x.local_path) # doctest: +ELLIPSIS
   ScoringFiles('.../PGS000001_hmPOS_GRCh38.txt.gz', target_build=None)

   But the ``target_build`` parameter doesn't work with local files:

   >>> with tempfile.TemporaryDirectory() as d:
   ...     x = ScoringFile("PGS000002", target_build=GenomeBuild.GRCh38)
   ...     x.download(directory=d)
   ...     ScoringFiles(x.local_path, target_build=GenomeBuild.GRCh37)
   Traceback (most recent call last):
   ...
   ValueError: Can't load local scoring file when target_build is setTry .normalise() method to do liftover, or load harmonised scoring files from PGS Catalog

   If you have a local scoring file that needs to change genome build, and using PGS
   Catalog harmonised data isn't an option, you should make a :class:`ScoringFile` from a path, then
   use the ``normalise()`` method with liftover enabled.

   .. py:property:: elements

      Returns a list of :class:`ScoringFile` objects contained inside :class:`ScoringFiles`



