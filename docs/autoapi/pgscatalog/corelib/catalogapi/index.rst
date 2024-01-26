:orphan:

:py:mod:`pgscatalog.corelib.catalogapi`
=======================================

.. py:module:: pgscatalog.corelib.catalogapi

.. autoapi-nested-parse::

   Classes and functions related to the PGS Catalog API 



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pgscatalog.corelib.catalogapi.CatalogCategory
   pgscatalog.corelib.catalogapi.CatalogQuery
   pgscatalog.corelib.catalogapi.ScoreQueryResult




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



