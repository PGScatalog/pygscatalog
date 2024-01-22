""" This module contains classes that compose a ScoringFile: a file in the
PGS Catalog that contains a list of genetic variants and their effect weights.
Scoring files are used to calculate PGS for new target genomes. """
import csv
import itertools
import logging
import pathlib

from .build import GenomeBuild
from .catalogapi import ScoreQueryResult, CatalogQuery
from ._normalise import normalise
from ._download import https_download
from ._config import Config
from .pgsexceptions import ScoreFormatError
from ._read import (
    generate_header_lines,
    auto_open,
    read_rows_lazy,
    get_columns,
    detect_wide,
)

logger = logging.getLogger(__name__)


class ScoringFileHeader:
    """Headers are a way of storing useful metadata about the scoring file. This
    header expects a PGS Catalog header format.

    It's always best to build headers with from_path():
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

        if self.pgs_name is None:
            raise ValueError("pgs_name cannot be None")

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

    >>> sf = ScoringFile(Config.ROOT_DIR / "tests" / "PGS000001_hmPOS_GRCh38.txt.gz")
    >>> sf.genome_build
    GenomeBuild.GRCh38
    >>> sf.pgs_id
    'PGS000001'
    >>> for variant in sf.variants: # doctest: +ELLIPSIS
    ...     variant
    ...     break
    ScoreVariant(effect_allele='T',effect_weight='0.16220387987485377',accession='PGS000001',...

    It's important to use the .download() method when you're not working with local files,
    or many attributes and methods will fail:
    >>> sf = ScoringFile("PGS000001")
    >>> for variant in sf.variants:
    ...     variant
    ...     break
    Traceback (most recent call last):
    ...
    ValueError: Local file is missing. Did you .download()?
    """

    __slots__ = [
        "pgs_id",
        "genome_build",
        "harmonised",
        "path",
        "local_path",
        "catalog_response",
        "include_children",
        "_identifier",
        "_header",
        "_directory",
        "_is_wide",
        "_start_line",
        "_fields",
    ]

    def __init__(self, identifier, target_build=None, query_result=None, **kwargs):
        if query_result is None:
            self._identifier = identifier
        else:
            self._identifier = query_result

        try:
            self._header = ScoringFileHeader.from_path(self._identifier)
        except (FileNotFoundError, TypeError):
            self.include_children = kwargs.get("include_children", None)
            self._init_from_accession(self._identifier, target_build=target_build)
        else:
            self.local_path = pathlib.Path(self._identifier)
            self._init_from_path(target_build=target_build)

        # set up local file attributes
        try:
            start_line, fields = get_columns(self.local_path)
        except TypeError:
            self._is_wide = None
            self._start_line = None
            self._fields = None
            self._directory = None
        else:
            self._is_wide = detect_wide(fields)
            self._start_line = start_line
            self._fields = fields
            self._directory = self.local_path.parent

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
        self.catalog_response = score
        self.path = score.get_download_url(target_build)
        self.local_path = None

    def _init_from_path(self, target_build=None):
        if target_build is not None:
            raise NotImplementedError("Can't liftover coordinates yet")

        self.path = self.local_path
        self.catalog_response = None

        if (pgs_id := self._header.pgs_id) is not None:
            self.pgs_id = pgs_id
        else:
            raise ScoreFormatError("Missing pgs_id from header")

        if (build := self._header.HmPOS_build) is not None:
            self.genome_build = build
            self.harmonised = True
        else:
            self.genome_build = self._header.genome_build
            self.harmonised = False

    def _generate_variants(self):
        """Yields rows from a scoring file as ScoreVariant objects"""

        row_nr = 0

        with auto_open(self.local_path, mode="rt") as f:
            for _ in range(self._start_line + 1):
                # skip header
                next(f)

            while True:
                batch = list(itertools.islice(f, Config.BATCH_SIZE))
                if not batch:
                    break

                csv_reader = csv.reader(batch, delimiter="\t")
                yield from read_rows_lazy(
                    csv_reader=csv_reader,
                    fields=self._fields,
                    name=self.pgs_id,
                    wide=self._is_wide,
                    row_nr=row_nr,
                )
                # this is important because row_nr resets for each batch
                row_nr += len(batch)

    @property
    def variants(self):
        if self.local_path is not None:
            return self._generate_variants()
        else:
            raise ValueError("Local file is missing. Did you .download()?")

    def __repr__(self):
        return f"{type(self).__name__}({repr(self._identifier)})"

    def __hash__(self):
        return hash(self.pgs_id)

    def __eq__(self, other):
        if isinstance(other, ScoringFile):
            return self.pgs_id == other.pgs_id
        return False

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
        corelib.pgsexceptions.InvalidAccessionError: Invalid accession: 'potato'

        >>> with tempfile.TemporaryDirectory() as tmp_dir:
        ...     ScoringFile("PGSinvalidaccession", GenomeBuild.GRCh38).download(tmp_dir)
        Traceback (most recent call last):
        ...
        corelib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGSinvalidaccession'
        """
        self._directory = pathlib.Path(directory)
        fn = pathlib.Path(self.path).name
        out_path = self._directory / fn
        https_download(
            url=self.path,
            out_path=out_path,
            overwrite=overwrite,
            directory=self._directory,
        )

        # update local file attributes
        self.local_path = out_path
        start_line, fields = get_columns(self.local_path)
        self._is_wide = detect_wide(fields)
        self._start_line = start_line
        self._fields = fields

    def normalise(self, drop_missing=False, liftover=False, **kwargs):
        """Extracts key fields from a scoring file in a normalised format.

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
        available from the PGS Catalog. This data goes through a lot of validation
        and error checking.

        For example, if you set the wrong genome build, you can get incorrect
        coordinates returned and _no errors or exceptions are raised_:

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

        A LiftoverError is only raised when many converted coordinates are missing.
        """
        yield from normalise(
            scoring_file=self, drop_missing=drop_missing, liftover=liftover, **kwargs
        )

    def get_log(self, drop_missing=False, variant_log=None):
        """Get a"""
        log = {}

        for attr in self._header.__slots__:
            if (extracted_attr := getattr(self._header, attr, None)) is not None:
                log[attr] = str(extracted_attr)
            else:
                log[attr] = None

        if log["variants_number"] is None:
            # custom scoring files might not have this information
            log["variants_number"] = variant_log["n_variants"]

        if (
            variant_log is not None
            and int(log["variants_number"]) != variant_log["n_variants"]
            and not drop_missing
        ):
            logger.warning(
                f"Mismatch between header ({log['variants_number']}) and output row count ({variant_log['n_variants']}) for {self.pgs_id}"
            )
            logger.warning(
                "This can happen with older scoring files in the PGS Catalog (e.g. PGS000028)"
            )

        # multiple terms may be separated with a pipe
        if log["trait_mapped"]:
            log["trait_mapped"] = log["trait_mapped"].split("|")

        if log["trait_efo"]:
            log["trait_efo"] = log["trait_efo"].split("|")

        log["columns"] = self._fields
        log["use_harmonised"] = self.harmonised

        if variant_log is not None:
            log["sources"] = [k for k, v in variant_log.items() if k != "n_variants"]

        return {self.pgs_id: log}


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
    Traceback (most recent call last):
    ...
    corelib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGPpotato'
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
