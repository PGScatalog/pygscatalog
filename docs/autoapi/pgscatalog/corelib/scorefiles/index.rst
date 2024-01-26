:orphan:

:py:mod:`pgscatalog.corelib.scorefiles`
=======================================

.. py:module:: pgscatalog.corelib.scorefiles

.. autoapi-nested-parse::

   This module contains classes to compose and contain a ``ScoringFile``: a file
   in the PGS Catalog that contains a list of genetic variants and their effect weights.
   Scoring files are used to calculate PGS for new target genomes.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.scorefiles.ScoringFile
   pgscatalog.corelib.scorefiles.ScoringFileHeader
   pgscatalog.corelib.scorefiles.ScoringFiles




Attributes
~~~~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.scorefiles.logger


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



.. py:class:: ScoringFileHeader(*, pgs_name, genome_build, pgs_id=None, pgp_id=None, variants_number=None, trait_reported=None, trait_efo=None, trait_mapped=None, weight_type=None, citation=None, HmPOS_build=None, HmPOS_date=None, format_version=None, license=None)


   Headers store useful metadata about a scoring file.

   This class provides convenient functions for reading and extracting information
   from the header. The header must follow PGS Catalog standards. It's always best
   to build headers with ``from_path()``:

   >>> testpath = Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt.gz"
   >>> ScoringFileHeader.from_path(testpath) # doctest: +ELLIPSIS
   ScoringFileHeader(pgs_id='PGS000001', pgp_id='PGP000001', pgs_name='PRS77_BC', ...

   But you can construct an instance with some minimum data:

   >>> header = ScoringFileHeader(pgs_name="PGS0000001", genome_build="hg19")
   >>> header # doctest: +ELLIPSIS
   ScoringFileHeader(pgs_id='None', pgp_id='None', pgs_name='PGS0000001', ...

   Strings are always used to construct (e.g. genome_build='GRCh37') but the header
   contains some objects:

   >>> header.genome_build
   GenomeBuild.GRCh37

   .. py:attribute:: fields
      :value: ('pgs_id', 'pgp_id', 'pgs_name', 'genome_build', 'variants_number', 'trait_reported',...

      

   .. py:method:: from_path(path)
      :classmethod:



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



.. py:data:: logger

   

