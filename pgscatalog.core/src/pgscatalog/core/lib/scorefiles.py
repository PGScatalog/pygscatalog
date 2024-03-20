""" This module contains classes to compose and contain a ``ScoringFile``: a file
in the PGS Catalog that contains a list of genetic variants and their effect weights.
Scoring files are used to calculate PGS for new target genomes."""
import csv
import itertools
import logging
import pathlib

from xopen import xopen

try:
    import pyarrow as pa
except ImportError:
    PYARROW_AVAILABLE = False
else:
    PYARROW_AVAILABLE = True

from .scorevariant import ScoreVariant
from .genomebuild import GenomeBuild
from .catalogapi import ScoreQueryResult, CatalogQuery
from ._normalise import normalise
from ._download import https_download
from ._config import Config
from .pgsexceptions import ScoreFormatError
from ._read import (
    read_rows_lazy,
    get_columns,
    detect_wide,
    read_header,
)

logger = logging.getLogger(__name__)


class ScoringFileHeader:
    """Headers store useful metadata about a scoring file.

    This class provides convenient functions for reading and extracting information
    from the header. The header must follow PGS Catalog standards. It's always best
    to build headers with ``from_path()``:

    >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
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

    fields = (
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
        values = {x: getattr(self, x) for x in self.fields}
        value_strings = ", ".join([f"{key}='{value}'" for key, value in values.items()])
        return f"{type(self).__name__}({value_strings})"

    @classmethod
    def from_path(cls, path):
        raw_header: dict = read_header(path)

        if len(raw_header) == 0:
            raise ValueError(f"No header detected in scoring file {path=}")

        header_dict = {k: raw_header.get(k) for k in cls.fields}

        return cls(**header_dict)


class ScoringFile:
    """Represents a single scoring file in the PGS Catalog.

    :param identifier: A PGS Catalog score accession in the format ``PGS123456`` or a path to a local scoring file
    :param target_build: An optional :class:`GenomeBuild`, which represents the build you want the scoring file to align to
    :param query_result: An optional :class:`ScoreQueryResult`, if provided with an accession identifier it prevents hitting the PGS Catalog API
    :raises pgscatalog.corelib.InvalidAccessionError: If the PGS Catalog API can't find the provided accession
    :raises pgscatalog.corelib.ScoreFormatError: If you try to iterate over a ``ScoringFile`` without a local path (before downloading it)

    You can make ``ScoringFiles`` with a path to a scoring file:

    >>> sf = ScoringFile(Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz")
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
    core.lib.pgsexceptions.ScoreFormatError: Local file is missing. Did you .download()?

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
    core.lib.pgsexceptions.InvalidAccessionError: Invalid accession: 'potato'

    The same exception is raised if you provide a well formatted identifier that doesn't exist:

    >>> with tempfile.TemporaryDirectory() as tmp_dir:
    ...     ScoringFile("PGS000000", GenomeBuild.GRCh38).download(tmp_dir)
    Traceback (most recent call last):
    ...
    core.lib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGS000000'
    """

    def __init__(self, identifier, target_build=None, query_result=None, **kwargs):
        self._target_build = target_build

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
            # remote file hasn't been downloaded yet
            self.is_wide = None
            self._start_line = None
            self._fields = None
            self._directory = None
        else:
            self.is_wide = detect_wide(fields)
            self._start_line = start_line
            self._fields = fields
            self._directory = self.local_path.parent

    def _init_from_accession(self, accession, target_build):
        logger.debug("Instantiating ScoringFile from accession")

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
            raise ScoreFormatError(
                f"Can't create a ScoringFile with accession: {accession!r}. "
                "Only PGS ids are supported. Try ScoringFiles()"
            )

        self.pgs_id = score.pgs_id
        self.catalog_response = score
        self.path = score.get_download_url(target_build)
        self.local_path = None

    def _init_from_path(self, target_build=None):
        logger.debug(f"Instantiating ScoringFile from {self.local_path=}")

        if target_build is not None:
            raise ValueError(
                "target_build must be None for local files. "
                "Use .normalise(target_build=...) if you want to liftover genomic coordinates"
            )

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

        with xopen(self.local_path, mode="rt") as f:
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
                    wide=self.is_wide,
                    row_nr=row_nr,
                )
                # this is important because row_nr resets for each batch
                row_nr += len(batch)

    @property
    def variants(self):
        """A generator that yields rows from the scoring file as ``ScoreVariants``,
        if a local file is available (i.e. after downloading). Always available for
        class instances that have a valid local path."""
        if self.local_path is not None:
            return self._generate_variants()
        else:
            raise ScoreFormatError("Local file is missing. Did you .download()?")

    @property
    def target_build(self):
        """The ``GenomeBuild`` you want a ``ScoringFile`` to align to. Useful when using PGS
        Catalog accessions to instantiate this class."""
        return self._target_build

    def __repr__(self):
        if self.local_path is not None:
            return f"{type(self).__name__}({repr(str(self.local_path))}, target_build={repr(self.target_build)})"
        else:
            return f"{type(self).__name__}({repr(self.pgs_id)}, target_build={repr(self.target_build)})"

    def __hash__(self):
        return hash(self.pgs_id)

    def __eq__(self, other):
        if isinstance(other, ScoringFile):
            return self.pgs_id == other.pgs_id
        return False

    def download(self, directory, overwrite=False):
        """
        Download a ScoringFile to a specified directory with checksum validation

        :param directory: Directory to write file to
        :param overwrite: Overwrite existing file if present

        :raises pgscatalog.corelib.ScoreDownloadError: If there's an unrecoverable problem downloading the file
        :raises pgscatalog.corelib.ScoreChecksumError: If md5 validation consistently fails

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

        """
        self._directory = pathlib.Path(directory)
        fn = pathlib.Path(self.path).name
        out_path = self._directory / fn

        logger.debug(f"Downloading {self.path} to {out_path}")

        https_download(
            url=self.path,
            out_path=out_path,
            overwrite=overwrite,
            directory=self._directory,
        )

        # update local file attributes
        self.local_path = out_path
        start_line, fields = get_columns(self.local_path)
        self.is_wide = detect_wide(fields)
        self._start_line = start_line
        self._fields = fields

    def normalise(
        self, liftover=False, drop_missing=False, chain_dir=None, target_build=None
    ):
        """Extracts key fields from a scoring file in a normalised format.

        Takes care of quality control.

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
        >>> variants = ScoringFile(testpath).normalise()
        >>> for x in variants: # doctest: +ELLIPSIS
        ...     x
        ...     break
        ScoreVariant(effect_allele='T',effect_weight='0.16220387987485377',...

        Supports lifting over scoring files from GRCh37 to GRCh38:

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh37.txt"
        >>> chaindir = Config.ROOT_DIR / "tests" / "data" / "chain"
        >>> sf = ScoringFile(testpath)
        >>> sf.harmonised = False  # lying, or liftover will be skipped
        >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh38)
        >>> for x in variants:
        ...     (x.rsID, x.chr_name, x.chr_position)
        ...     break
        ('rs78540526', '11', 69516650)

        Example of lifting down (GRCh38 to GRCh37):

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt"
        >>> chaindir = Config.ROOT_DIR / "tests" / "data" / "chain"
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

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt"
        >>> chaindir = Config.ROOT_DIR / "tests" / "data" / "chain"
        >>> sf = ScoringFile(testpath)
        >>> sf.harmonised = False  # lying, or liftover will be skipped
        >>> sf.genome_build = GenomeBuild.GRCh37  # wrong build ! it's GRCh38
        >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh38)
        >>> for x in variants:
        ...     (x.rsID, x.chr_name, x.chr_position)
        ...     break
        ('rs78540526', '11', 69701882)

        A :class:`LiftoverError` is only raised when many converted coordinates are missing.

        Normalising converts the is_dominant and is_recessive optional fields in
        scoring files into an EffectType:

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000802_hmPOS_GRCh37.txt"
        >>> variants = ScoringFile(testpath).normalise()
        >>> for i, x in enumerate(variants): # doctest: +ELLIPSIS
        ...     (x.is_dominant, x.is_recessive, x.effect_type)
        ...     if i == 2:
        ...         break
        (True, False, EffectType.DOMINANT)
        (False, True, EffectType.RECESSIVE)
        (True, False, EffectType.DOMINANT)
        """
        yield from normalise(
            scoring_file=self,
            drop_missing=drop_missing,
            liftover=liftover,
            chain_dir=chain_dir,
            target_build=target_build,
        )

    def get_log(self, drop_missing=False, variant_log=None):
        """Create a JSON log from a ScoringFile's header and variant rows."""

        logger.debug(f"Creating JSON log for {self!r}")

        log = {}

        for attr in self._header.fields:
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
    """This container class provides methods to work with multiple ScoringFile objects.

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
    core.lib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGPpotato'

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
    """

    def __init__(self, *args, target_build=None, **kwargs):
        self.target_build = target_build
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
                case ScoringFile() if arg.target_build == target_build:
                    logger.info("ScoringFile build matches target build")
                    scorefiles.append(arg)
                case ScoringFile() if arg.target_build != target_build:
                    raise ValueError(
                        f"{arg.target_build=} doesn't match {target_build=}"
                    )
                case _ if pathlib.Path(arg).is_file() and target_build is None:
                    logger.info(f"Local path: {arg}, no target build is OK")
                    scorefiles.append(ScoringFile(arg))
                case _ if pathlib.Path(arg).is_file() and target_build is not None:
                    logger.critical(f"{arg} is a local file and {target_build=}")
                    raise ValueError(
                        "Can't load local scoring file when target_build is set"
                        "Try .normalise() method to do liftover, or load harmonised scoring files from PGS Catalog"
                    )
                case str() if arg.startswith("PGP") or "_" in arg:
                    logger.info(
                        "Term associated with multiple scores detected (PGP or trait)"
                    )
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
                    logger.info("PGS ID detected")
                    pgs_batch.append(arg)
                case str():
                    raise ValueError(f"{arg!r} is not a valid path or an accession")
                case _:
                    raise TypeError

        # batch PGS IDs to avoid overloading the API
        batched_queries = CatalogQuery(accession=pgs_batch).score_query()
        logger.debug(f"Batching queries to PGS Catalog API: {pgs_batch}")
        batched_scores = [
            ScoringFile(x, target_build=target_build) for x in batched_queries
        ]
        scorefiles.extend(batched_scores)

        self._elements = list(dict.fromkeys(scorefiles))

    def __repr__(self):
        ids = []
        for x in self.elements:
            if x.local_path is not None:
                ids.append(str(x.local_path))
            else:
                ids.append(x.pgs_id)

        args = ", ".join([repr(x) for x in ids])
        return f"{type(self).__name__}({args}, target_build={repr(self.target_build)})"

    def __iter__(self):
        return iter(self.elements)

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, item):
        return self.elements[item]

    def __contains__(self, item):
        for element in self.elements:
            if element == item:
                return True
        return False

    def __add__(self, other):
        if isinstance(other, type(self)):
            if self.target_build == other.target_build:
                new_elements = self._elements + other.elements
                return ScoringFiles(new_elements, target_build=self.target_build)
            else:
                return NotImplemented
        else:
            return NotImplemented

    def __mul__(self, other):
        """Intentionally not implemented. Cannot contain duplicate elements."""
        return NotImplemented

    @property
    def elements(self):
        """Returns a list of :class:`ScoringFile` objects contained inside :class:`ScoringFiles`"""
        return self._elements


def _read_normalised_rows(path):
    with xopen(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield ScoreVariant(**row)


class NormalisedScoringFile:
    """This class represents a ScoringFile that's been normalised to have a consistent format

    Its main purpose is to provide a convenient way to iterate over variants

    >>> normpath = Config.ROOT_DIR / "tests" / "data" / "combined.txt.gz"
    >>> test = NormalisedScoringFile(normpath)
    >>> test  # doctest: +ELLIPSIS
    NormalisedScoringFile('.../combined.txt.gz')

    >>> for i in test.variants:  # doctest: +ELLIPSIS
    ...     i
    ...     break
    ScoreVariant(effect_allele='T',effect_weight='0.2104229869048294',...

    >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
    >>> test = NormalisedScoringFile(ScoringFile(testpath))
    >>> test  # doctest: +ELLIPSIS
    NormalisedScoringFile(ScoringFile('.../PGS000001_hmPOS_GRCh38.txt.gz', ...))


    >>> for i in test.variants:  # doctest: +ELLIPSIS
    ...     i
    ...     break
    ScoreVariant(effect_allele='T',effect_weight='0.16220387987485377',...

    >>> for x in test.to_pa_recordbatch():  # doctest: +ELLIPSIS
    ...     x
    ...     break
    pyarrow.RecordBatch
    chr_name: string
    ...
    row_nr: [0,1,2,3,4,5,6,7,8,9,...,67,68,69,70,71,72,73,74,75,76]
    """

    def __init__(self, path):
        try:
            with xopen(path):
                pass
        except TypeError:
            self.path = False
        else:
            self.path = True
        finally:
            # either a ScoringFile or a path to a combined file
            self._scoringfile = path

    def __iter__(self):
        yield from self.variants

    @property
    def variants(self):
        if self.path:
            # get a fresh generator from the file
            self._variants = _read_normalised_rows(self._scoringfile)
        else:
            # get a fresh generator from the normalise() method
            self._variants = self._scoringfile.normalise()

        return self._variants

    def __repr__(self):
        if self.path:
            x = f"{repr(str(self._scoringfile))}"
        else:
            x = f"{repr(self._scoringfile)}"

        return f"{type(self).__name__}({x})"

    def pa_schema(self):
        """Return a pyarrow schema used when writing object to files"""
        if not PYARROW_AVAILABLE:
            raise ImportError("pyarrow is not available")

        return pa.schema(
            [
                pa.field("chr_name", pa.string()),
                pa.field("chr_position", pa.uint64()),
                pa.field("effect_allele", pa.string()),
                pa.field("other_allele", pa.string()),
                pa.field("effect_weight", pa.string()),
                pa.field("effect_type", pa.string()),
                pa.field("is_duplicated", pa.bool_()),
                pa.field("accession", pa.string()),
                pa.field("row_nr", pa.uint64()),
            ]
        )

    def to_pa_recordbatch(self):
        """Yields an iterator of pyarrow RecordBatches"""
        if not PYARROW_AVAILABLE:
            raise ImportError("pyarrow is not available")

        # accessing the property returns a fresh generator
        # don't want an infinite loop
        variants = self.variants

        while True:
            batch = list(itertools.islice(variants, Config.TARGET_BATCH_SIZE))
            if not batch:
                break

            _chrom, _pos, _ea, _oa, _ew, _et, _dup, _acc, _rownr = zip(
                *(
                    (
                        x.chr_name,
                        x.chr_position,
                        str(x.effect_allele),
                        x.other_allele,
                        x.effect_weight,
                        str(x.effect_type),
                        x.is_duplicated,
                        x.accession,
                        x.row_nr,
                    )
                    for x in batch
                ),
                strict=True,
            )

            # convert truthy strings to bools
            _dup = (True if x == "True" else False for x in _dup)

            yield pa.RecordBatch.from_arrays(
                [_chrom, _pos, _ea, _oa, _ew, _et, _dup, _acc, _rownr],
                schema=self.pa_schema(),
            )
