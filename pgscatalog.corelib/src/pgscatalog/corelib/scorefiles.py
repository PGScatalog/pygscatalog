""" This module contains classes that compose a ScoringFile: a file in the
PGS Catalog that contains a list of genetic variants and their effect weights.
Scoring files are used to calculate PGS for new target genomes. """

import gzip
import hashlib
import os
import pathlib
import tempfile

import httpx
import tenacity
from tenacity import retry

from pgscatalog.corelib import config
from pgscatalog.corelib.catalogapi import CatalogQuery, GenomeBuild


class ScoringFileHeader:
    """Headers are a way of storing useful metadata about the scoring file. This
    header expects a PGS Catalog header format.

    It's always best to build headers with from_path():
    >>> import pgscatalog, importlib.resources
    >>> testpath = importlib.resources.files(pgscatalog).joinpath(".").parent.parent / "tests" / "testdata" / "PGS000001_hmPOS_GRCh38.txt.gz"
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
    def __init__(self, identifier, target_build=None):
        try:
            self.header = ScoringFileHeader.from_path(identifier)
        except FileNotFoundError:
            self._init_from_accession(identifier, target_build)
        else:
            # init from file path
            if target_build:
                # todo: liftover
                raise ValueError("Don't set a target build with a path")

            self.path = identifier
            self.catalog_response = None

    def _init_from_accession(self, accession, target_build):
        score = CatalogQuery(accession=accession).score_query()
        try:
            len(score)  # was a list returned?
        except TypeError:
            pass  # just a normal ScoreQueryResult
        else:
            # a list of score results can't be used to make a ScoringFile
            raise ValueError(f"Can't create a ScoringFile with accession: {accession!r}. "
                             "Only PGS ids are supported. Try ScoringFiles()")

        self.header = None
        self.catalog_response = score
        self.path = score.get_download_url(target_build)


    @retry(stop=tenacity.stop_after_attempt(5), retry=tenacity.retry_if_exception_type(httpx.RequestError))
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
        """
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
                    raise httpx.RequestError(f"Calculated checksum {calc} doesn't match {remote}")
                else:
                    os.rename(f.name, out_path)
        except httpx.UnsupportedProtocol:
            raise ValueError(f"Can't download a local file: {self.path!r}")


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