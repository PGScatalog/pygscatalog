"""Classes and functions related to the PGS Catalog API"""

import enum
import logging

import httpx
import tenacity

from pgscatalog.core.lib.pgsexceptions import QueryError, InvalidAccessionError
from pgscatalog.core.lib.genomebuild import GenomeBuild
from pgscatalog.core.lib._config import Config


logger = logging.getLogger(__name__)


class CatalogCategory(enum.Enum):
    """The three main categories in the PGS Catalog

    Enumeration values don't mean anything and are automatically generated:

    >>> CatalogCategory.SCORE
    <CatalogCategory.SCORE: 1>
    """

    SCORE = enum.auto()
    TRAIT = enum.auto()
    PUBLICATION = enum.auto()


def _query_error(retry_state):
    """Couldn't query the PGS Catalog API after retrying and waiting a bunch"""
    try:
        retry_state.outcome.result()
    except Exception as e:
        raise QueryError("Can't query PGS Catalog API") from e


class CatalogQuery:
    """Efficiently query the PGS Catalog API using accessions

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
    """

    _rest_url_root = "https://www.pgscatalog.org/rest"

    def __init__(self, *, accession, include_children=False, **kwargs):
        if not isinstance(accession, str):
            # deduplicate lists, preserving input order
            self.accession = list(dict.fromkeys(accession))
        else:
            self.accession = accession

        self.category = kwargs.setdefault("category", self.infer_category())

        match self.category:
            case CatalogCategory.TRAIT:
                # only traits can have children
                self.include_children = include_children
            case CatalogCategory.SCORE | CatalogCategory.PUBLICATION:
                self.include_children = None
            case _:
                raise ValueError(f"Bad category: {self.category!r}")

    def __repr__(self):
        return (
            f"{type(self).__name__}(accession={repr(self.accession)}, category="
            f"{self.category}, include_children={self.include_children})"
        )

    def infer_category(self):
        """Inspect an accession and guess the Catalog category

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
        """
        match accession := self.accession:
            case list() if all([x.startswith("PGS") for x in accession]):
                category = CatalogCategory.SCORE
            case list():
                raise ValueError(
                    f"Invalid accession in list: {accession!r}. Lists must only contain PGS IDs."
                )
            case str() if accession.startswith("PGS"):
                category = CatalogCategory.SCORE
            case str() if accession.startswith("PGP"):
                category = CatalogCategory.PUBLICATION
            case str() if "_" in accession:
                # simple check for structured text like EFO_ACCESSION, HP_ACCESSION, etc
                category = CatalogCategory.TRAIT
            case _:
                raise InvalidAccessionError(f"Invalid accession: {accession!r}")

        logger.debug(f"{accession=} {category=}")
        return category

    def get_query_url(self):
        """
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
        """
        urls: list[str] | str = []
        match (self.category, self.accession):
            case CatalogCategory.TRAIT, str():
                child_flag: int = int(self.include_children)
                urls = f"{self._rest_url_root}/trait/{self.accession}?include_children={child_flag}"
            case CatalogCategory.SCORE, str():
                urls = [f"{self._rest_url_root}/score/search?pgs_ids={self.accession}"]
            case CatalogCategory.SCORE, list():
                for chunk in self._chunk_accessions():
                    chunked_accession = ",".join(chunk)
                    urls.append(
                        f"{self._rest_url_root}/score/search?pgs_ids="
                        f"{chunked_accession}"
                    )
                return urls
            case CatalogCategory.PUBLICATION, str():
                urls = f"{self._rest_url_root}/publication/{self.accession}"
            case _:
                raise ValueError(
                    f"Invalid CatalogCategory and accession type: {self.category!r}, "
                    f"type({self.accession!r})"
                )
        logger.debug(f"Resolved API query URL: {urls}")
        return urls

    def _chunk_accessions(self):
        size = 50  # /rest/score/{pgs_id} limit when searching multiple IDs
        # using a dict to get unique elements instead of a set to preserve order
        accessions = self.accession
        return (accessions[pos : pos + size] for pos in range(0, len(accessions), size))

    @tenacity.retry(
        stop=tenacity.stop_after_attempt(Config.MAX_RETRIES),
        retry=tenacity.retry_if_exception_type(httpx.RequestError),
        retry_error_callback=_query_error,
        wait=tenacity.wait_fixed(3) + tenacity.wait_random(0, 2),
    )
    def score_query(self):
        """Query the PGS Catalog API and return :class:`ScoreQueryResult`

        Information about a single score is returned as a dict:

        >>> CatalogQuery(accession="PGS000001").score_query() # doctest: +ELLIPSIS
        ScoreQueryResult(pgs_id='PGS000001', ftp_url=...

        If information about multiple scores is found, it's returned as a list:

        >>> CatalogQuery(accession=["PGS000001", "PGS000002"]).score_query() # doctest: +ELLIPSIS
        [ScoreQueryResult(pgs_id='PGS000001', ftp_url=...

        Publications and traits always return a list of score information:

        >>> CatalogQuery(accession="PGP000001").score_query() # doctest: +ELLIPSIS
        [ScoreQueryResult(pgs_id='PGS000001', ftp_url=...
        """
        match self.category:
            case CatalogCategory.SCORE:
                results = []

                for url in self.get_query_url():
                    r = httpx.get(url, timeout=5, headers=Config.API_HEADER).json()

                    if "request limit exceeded" in r.get("message", ""):
                        raise httpx.RequestError("request limit exceeded")
                    else:
                        results += r["results"]

                # return the same type as the accession input to be consistent
                match self.accession:
                    case list():
                        return ScoreQueryResult.from_query(results)
                    case str():
                        # a PGS string accession input can only ever return one result
                        try:
                            return ScoreQueryResult.from_query(results[0])
                        except IndexError:
                            raise InvalidAccessionError(
                                f"No Catalog result for accession {self.accession!r}"
                            )
                    case _:
                        raise ValueError
            case CatalogCategory.PUBLICATION:
                url = self.get_query_url()
                r = httpx.get(url, timeout=5, headers=Config.API_HEADER).json()
                try:
                    pgs_ids = [
                        score
                        for scores in list(r["associated_pgs_ids"].values())
                        for score in scores
                    ]
                except KeyError:
                    raise InvalidAccessionError(
                        f"No Catalog result for accession {self.accession!r}"
                    )
                else:
                    return CatalogQuery(accession=pgs_ids).score_query()
            case CatalogCategory.TRAIT:
                url = self.get_query_url()
                r = httpx.get(url, timeout=5, headers=Config.API_HEADER).json()
                pgs_ids = r["associated_pgs_ids"]
                if self.include_children:
                    pgs_ids.extend(r["child_associated_pgs_ids"])
                return CatalogQuery(accession=pgs_ids).score_query()


class ScoreQueryResult:
    """Class that holds score metadata with methods to extract important fields"""

    def __init__(self, *, pgs_id, ftp_url, ftp_grch37_url, ftp_grch38_url, license):
        self.pgs_id = pgs_id
        self.ftp_url = ftp_url
        self.ftp_grch37_url = ftp_grch37_url
        self.ftp_grch38_url = ftp_grch38_url
        self.license = license

    def __repr__(self):
        return (
            f"{type(self).__name__}(pgs_id={self.pgs_id!r}, ftp_url={self.ftp_url!r}, "
            f"ftp_grch37_url={self.ftp_grch37_url!r},ftp_grch38_url={self.ftp_grch38_url!r},"
            f"license={self.license!r})"
        )

    @classmethod
    def from_query(cls, result_response):
        """
        Parses PGS Catalog API JSON response

        :param result_response: PGS Catalog API JSON response
        :returns: :class:`ScoreQueryResult`

        >>> fake_response = {"id": "fake", "ftp_harmonized_scoring_files":
        ... {"GRCh37": {"positions": "fake.txt.gz"}, "GRCh38": {"positions": "fake.txt.gz"}},
        ... "license": "fake", "ftp_scoring_file": "fake.txt.gz"}
        >>> ScoreQueryResult.from_query(fake_response) # doctest: +ELLIPSIS
        ScoreQueryResult(pgs_id='fake', ftp_url='fake.txt.gz',...
        """
        try:
            pgs_id = result_response["id"]
        except TypeError:
            # assume result_response is a list of responses
            results = []
            for result in result_response:
                results.append(ScoreQueryResult.from_query(result))
            return results
        else:
            ftp_url = result_response["ftp_scoring_file"]
            harmonised_urls = result_response["ftp_harmonized_scoring_files"]

            ftp_grch37_url = harmonised_urls["GRCh37"]["positions"]
            ftp_grch38_url = harmonised_urls["GRCh38"]["positions"]
            license = result_response["license"]
            return cls(
                pgs_id=pgs_id,
                ftp_url=ftp_url,
                ftp_grch37_url=ftp_grch37_url,
                ftp_grch38_url=ftp_grch38_url,
                license=license,
            )

    def get_download_url(self, genome_build=None):
        """
        Returns scoring file download URL, with support for specifying harmonised data in a specific genome build

        >>> query = CatalogQuery(accession="PGS000001").score_query()
        >>> build = GenomeBuild.GRCh38
        >>> query.get_download_url()
        'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/PGS000001.txt.gz'
        >>> query.get_download_url(build)
        'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/Harmonized/PGS000001_hmPOS_GRCh38.txt.gz'
        """
        match build := genome_build:
            case GenomeBuild() if build == GenomeBuild.GRCh37:
                url = self.ftp_grch37_url
            case GenomeBuild() if build == GenomeBuild.GRCh38:
                url = self.ftp_grch38_url
            case None:
                url = self.ftp_url
            case _:
                raise ValueError(f"Invalid genome build {build!r}")

        logger.debug(
            f"Scoring file download URL for {self.pgs_id} with {build=}: {url}"
        )
        return url
