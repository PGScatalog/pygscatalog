""" This module contains classes that compose a ScoringFile: a file in the
PGS Catalog that contains a list of genetic variants and their effect weights.
Scoring files are used to calculate PGS for new target genomes. """
import gzip
import hashlib
import itertools
import os
import pathlib
import tempfile
import urllib

import httpx
import tenacity
from tenacity import retry

from pgscatalog.corelib import config, pgsexceptions
from pgscatalog.corelib.catalogapi import CatalogQuery, GenomeBuild, ScoreQueryResult


def _score_download_failed(retry_state):
    """Every attempt to download went wrong, let's be specific about what went wrong"""
    try:
        retry_state.outcome.result()
    except pgsexceptions.ScoreChecksumError as e:
        # the checksum kept could never be matched
        raise e
    except Exception as download_exc:
        # stuff all other exceptions into a ScoreDownloadError
        raise pgsexceptions.ScoreDownloadError(
            "Downloading score failed"
        ) from download_exc


@retry(
    stop=tenacity.stop_after_attempt(config.MAX_ATTEMPTS),
    retry=tenacity.retry_if_exception_type(
        (pgsexceptions.ScoreDownloadError, pgsexceptions.ScoreChecksumError)
    ),
    retry_error_callback=_score_download_failed,
    wait=tenacity.wait_fixed(3) + tenacity.wait_random(0, 2),
)
def _ftp_fallback(retry_state):
    """When ScoringFile.download() fails, it invokes this callback function

    Try downloading from PGS Catalog using FTP protocol instead of HTTPS.
    """
    scorefile = retry_state.args[0]
    directory = retry_state.args[1]

    ftp_url = scorefile.path.replace("https://", "ftp://")
    checksum_url = (scorefile.path + ".md5").replace("https://", "ftp://")

    fn = pathlib.Path(scorefile.path).name
    out_path = pathlib.Path(directory) / fn
    md5 = hashlib.md5()

    try:
        with (
            tempfile.NamedTemporaryFile(dir=directory, delete=False) as score_f,
            tempfile.NamedTemporaryFile(dir=directory, delete=True) as checksum_f,
        ):
            urllib.request.urlretrieve(ftp_url, score_f.name)
            urllib.request.urlretrieve(checksum_url, checksum_f.name)

            md5.update(score_f.read())

            if (checksum := md5.hexdigest()) != (
                remote := checksum_f.read().decode().split()[0]
            ):
                raise pgsexceptions.ScoreChecksumError(
                    f"Local checksum {checksum} doesn't match remote {remote}"
                )
            else:
                os.rename(score_f.name, out_path)
    except urllib.error.URLError as download_exc:
        raise pgsexceptions.ScoreDownloadError(
            f"Can't download {scorefile} over FTP"
        ) from download_exc


class ScoringFileHeader:
    """Headers are a way of storing useful metadata about the scoring file. This
    header expects a PGS Catalog header format.

    It's always best to build headers with from_path():
    >>> import pgscatalog.corelib.config
    >>> testpath = config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt.gz"
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
    """

    # slots are used here because we want a controlled vocabulary
    # random extra attributes would be bad without thinking about them
    __slots__ = (
        "pgs_id",
        "pgp_id",
        "pgs_name",
        "genome_build",
        "variants_number",
        "trait_reported",
        "trait_efo",
        "trait_mapped",
        "weight_type",
        "citation",
        "HmPOS_build",
        "HmPOS_date",
        "format_version",
        "license",
    )

    _default_license_text = (
        "PGS obtained from the Catalog should be cited appropriately, and "
        "used in accordance with any licensing restrictions set by the authors. See "
        "EBI "
        "Terms of Use (https://www.ebi.ac.uk/about/terms-of-use/) for additional "
        "details."
    )

    def __init__(
        self,
        *,
        pgs_name,
        genome_build,
        pgs_id=None,
        pgp_id=None,
        variants_number=None,
        trait_reported=None,
        trait_efo=None,
        trait_mapped=None,
        weight_type=None,
        citation=None,
        HmPOS_build=None,
        HmPOS_date=None,
        format_version=None,
        license=None,
    ):
        """kwargs are forced because this is a complicated init and from_path() is
        almost always the correct thing to do.

        Mandatory/optional fields in a header are less clear than columns. The
        Catalog provides this information automatically but scoring files from other
        places might not.

        We don't want to annoy people and force them to reformat their custom scoring
        files, but we do require some minimum information for the calculator,
        so pgs_name and genome_build are mandatory.
        """
        self.pgs_name = pgs_name
        self.genome_build = GenomeBuild.from_string(genome_build)

        if self.genome_build is None:
            # try overwriting with harmonised data
            self.genome_build = GenomeBuild.from_string(HmPOS_build)

        if self.pgs_name is None or self.genome_build is None:
            raise ValueError("pgs_name and genome_build cannot be None")

        # the rest of these fields are optional
        self.pgs_id = pgs_id
        self.pgp_id = pgp_id

        try:
            self.variants_number = int(variants_number)
        except TypeError:
            self.variants_number = None

        self.trait_reported = trait_reported
        self.trait_efo = trait_efo
        self.trait_mapped = trait_mapped
        self.weight_type = weight_type
        self.citation = citation
        self.HmPOS_build = GenomeBuild.from_string(HmPOS_build)
        self.HmPOS_date = HmPOS_date
        self.format_version = format_version
        self.license = license

        if self.license is None:
            self.license = self._default_license_text

    def __repr__(self):
        values = {x: getattr(self, x) for x in self.__slots__}
        value_strings = ", ".join([f"{key}='{value}'" for key, value in values.items()])
        return f"{type(self).__name__}({value_strings})"

    @classmethod
    def from_path(cls, path):
        raw_header: dict = read_header(path)

        if len(raw_header) == 0:
            raise ValueError(f"No header detected in scoring file {path=}")

        header_dict = {k: raw_header.get(k) for k in cls.__slots__}

        return cls(**header_dict)


class ScoringFile:
    """Represents a single scoring file.

    Can also be constructed with a ScoreQueryResult to avoid hitting the API during instantiation
    """

    def __init__(self, identifier, target_build=None, query_result=None, **kwargs):
        if query_result is None:
            self._identifier = identifier
        else:
            self._identifier = query_result

        try:
            self.header = ScoringFileHeader.from_path(self._identifier)
        except (FileNotFoundError, TypeError):
            self.include_children = kwargs.get("include_children", None)
            self._init_from_accession(self._identifier, target_build=target_build)
        else:
            # init from file path
            if target_build:
                raise NotImplementedError("Can't liftover")
            raise NotImplementedError("Local files not supported")

    def __repr__(self):
        return f"{type(self).__name__}({repr(self._identifier)})"

    def __hash__(self):
        return hash(self.pgs_id)

    def __eq__(self, other):
        if isinstance(other, ScoringFile):
            return self.pgs_id == other.pgs_id
        return False

    def _init_from_accession(self, accession, target_build):
        match self._identifier:
            case ScoreQueryResult():
                # skip hitting the API unnecessarily
                score = self._identifier
                self._identifier = self._identifier.pgs_id
            case str():
                score = CatalogQuery(
                    accession=accession, include_children=self.include_children
                ).score_query()
            case _:
                raise TypeError(f"Can't init from accession: {self._identifier!r}")

        try:
            len(score)  # was a list returned from the Catalog query?
        except TypeError:
            pass  # just a normal ScoreQueryResult, continue
        else:
            # this class can only instantiate and represent one scoring file
            raise ValueError(
                f"Can't create a ScoringFile with accession: {accession!r}. "
                "Only PGS ids are supported. Try ScoringFiles()"
            )

        self.pgs_id = score.pgs_id
        self.header = None
        self.catalog_response = score
        self.path = score.get_download_url(target_build)

    @retry(
        stop=tenacity.stop_after_attempt(config.MAX_ATTEMPTS),
        retry=tenacity.retry_if_exception_type(
            (pgsexceptions.ScoreDownloadError, pgsexceptions.ScoreChecksumError)
        ),
        retry_error_callback=_ftp_fallback,
        wait=tenacity.wait_fixed(3) + tenacity.wait_random(0, 2),
    )
    def download(self, directory, overwrite=False):
        """
        Download a ScoringFile to a specified directory with checksum validation

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

        >>> with tempfile.TemporaryDirectory() as tmp_dir:
        ...     ScoringFile("potato", GenomeBuild.GRCh38).download(tmp_dir)
        Traceback (most recent call last):
        ...
        pgscatalog.corelib.pgsexceptions.InvalidAccessionError: Invalid accession: 'potato'

        >>> with tempfile.TemporaryDirectory() as tmp_dir:
        ...     ScoringFile("PGSinvalidaccession", GenomeBuild.GRCh38).download(tmp_dir)
        Traceback (most recent call last):
        ...
        pgscatalog.corelib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGSinvalidaccession'
        """
        if config.FTP_EXCLUSIVE:
            # replace wait function to hit callback quickly
            self.download.retry.wait = tenacity.wait_none()
            raise pgsexceptions.ScoreDownloadError(
                "HTTPS downloads disabled by config.FTP_EXCLUSIVE"
            )

        try:
            fn = pathlib.Path(self.path).name
            out_path = pathlib.Path(directory) / fn

            if out_path.exists() and not overwrite:
                raise FileExistsError(f"{out_path} already exists")

            checksum_path = self.path + ".md5"
            checksum = httpx.get(checksum_path, headers=config.API_HEADER).text
            md5 = hashlib.md5()

            with tempfile.NamedTemporaryFile(dir=directory, delete=False) as f:
                with httpx.stream("GET", self.path, headers=config.API_HEADER) as r:
                    for data in r.iter_bytes():
                        f.write(data)
                        md5.update(data)

                if (calc := md5.hexdigest()) != (remote := checksum.split()[0]):
                    # will attempt to download again (see decorator)
                    raise pgsexceptions.ScoreChecksumError(
                        f"Calculated checksum {calc} doesn't match {remote}"
                    )
                else:
                    os.rename(f.name, out_path)
        except httpx.UnsupportedProtocol as protocol_exc:
            raise ValueError(
                f"Can't download a local file: {self.path!r}"
            ) from protocol_exc
        except httpx.RequestError as download_exc:
            raise pgsexceptions.ScoreDownloadError(
                "HTTPS download failed"
            ) from download_exc


class ScoringFiles:
    """This class provides methods to work with multiple ScoringFile objects.

    You can use publications or trait accessions to instantiate:
    >>> pub = ScoringFiles("PGP000001")
    >>> len(pub)
    3

    Or multiple PGS IDs:
    >>> score = ScoringFiles("PGS000001", "PGS000002")
    >>> len(score)
    2

    List input is OK too:
    >>> score = ScoringFiles(["PGS000001", "PGS000002"])
    >>> len(score)
    2

    Or any mixture of publications, traits, and scores:
    >>> score = ScoringFiles("PGP000001", "PGS000001", "PGS000002")
    >>> len(score)
    3

    Scoring files with duplicate PGS IDs (accessions) are automatically dropped.
    In the example above PGP000001 contains PGS000001, PGS000002, and PGS000003.

    Traits can have children. To include these traits, use the include_children parameter:

    >>> score_with_children = ScoringFiles("MONDO_0004975", include_children=True)
    >>> score_wo_children = ScoringFiles("MONDO_0004975", include_children=False)
    >>> len(score_with_children) > len(score_wo_children)
    True

    For example, Alzheimer's disease (MONDO_0004975) includes Late-onset Alzheier's disease (EFO_1001870) as a child trait.

    You can slice and iterate over ScoringFiles:
    >>> score[0]
    ScoringFile('PGS000001')
    >>> for x in score:
    ...     x
    ScoringFile('PGS000001')
    ScoringFile('PGS000002')
    ScoringFile('PGS000003')

    >>> ScoringFiles("PGPpotato")
    """

    def __init__(self, *args, target_build=None, **kwargs):
        # flatten args to provide a more flexible interface
        flargs = list(
            itertools.chain.from_iterable(
                arg if isinstance(arg, list) else [arg] for arg in args
            )
        )
        scorefiles = []
        pgs_batch = []
        for arg in flargs:
            match arg:
                case _ if pathlib.Path(arg).is_file():
                    raise NotImplementedError
                case str() if arg.startswith("PGP") or "_" in arg:
                    self.include_children = kwargs.get("include_children", None)
                    traitpub_query = CatalogQuery(
                        accession=arg, include_children=self.include_children
                    ).score_query()

                    # avoid unnecessary API hits by using CatalogQuery objects
                    scorefiles.extend(
                        [
                            ScoringFile(x, target_build=target_build)
                            for x in traitpub_query
                        ]
                    )
                case str() if arg.startswith("PGS"):
                    pgs_batch.append(arg)
                case str():
                    raise ValueError(f"{arg!r} is not a valid path or an accession")
                case _:
                    raise TypeError

        # batch PGS IDs to avoid overloading the API
        batched_queries = CatalogQuery(accession=pgs_batch).score_query()
        batched_scores = [
            ScoringFile(x, target_build=target_build) for x in batched_queries
        ]
        scorefiles.extend(batched_scores)

        self._elements = list(dict.fromkeys(scorefiles))

    def __repr__(self):
        args = ", ".join([repr(x.pgs_id) for x in self.elements])
        return f"{type(self).__name__}({args})"

    def __iter__(self):
        return iter(self.elements)

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, item):
        return self.elements[item]

    @property
    def elements(self):
        return self._elements

    def combine(self):
        """Combining multiple scoring files yields ScoreVariants in a consistent genome build and data format.

        This process takes care of data munging and some quality control steps."""
        raise NotImplementedError


def read_header(path: pathlib.Path):
    """Parses the header of a PGS Catalog format scorefile into a dictionary"""
    header = {}

    with auto_open(path, "rt") as f:
        header_text = generate_header_lines(f)

        for item in header_text:
            key, value = item.split("=")
            header[key[1:]] = value  # drop # character from key

    return header


def generate_header_lines(f):
    """Header lines in a PGS Catalog scoring file are structured like:
        #pgs_id=PGS000348
        #pgs_name=PRS_PrCa
    Files can be big, so we want to only read header lines and stop immediately
    """
    for line in f:
        if line.startswith("#"):
            if "=" in line:
                yield line.strip()
        else:
            # stop reading lines
            break


def auto_open(filepath, mode="rt"):
    """Automatically open a gzipped text file or an uncompressed text file"""
    with open(filepath, "rb") as test_f:
        if test_f.read(2) == b"\x1f\x8b":
            return gzip.open(filepath, mode)
        else:
            return open(filepath, mode)
