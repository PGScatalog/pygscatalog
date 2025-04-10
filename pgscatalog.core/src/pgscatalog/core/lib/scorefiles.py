"""This module contains classes to compose and contain a ``ScoringFile``: a file
in the PGS Catalog that contains a list of genetic variants and their effect weights.
Scoring files are used to calculate PGS for new target genomes."""

import csv
import itertools
import logging
import pathlib
from contextlib import contextmanager
from typing import Iterator, Generator

from pydantic import ValidationError
from xopen import xopen

from pgscatalog.core.lib import GenomeBuild
from pgscatalog.core.lib import models
from pgscatalog.core.lib.catalogapi import ScoreQueryResult, CatalogQuery
from pgscatalog.core.lib._normalise import normalise
from pgscatalog.core.lib._download import https_download
from pgscatalog.core.lib.models import ScoreVariant
from pgscatalog.core.lib.pgsexceptions import ScoreFormatError

logger = logging.getLogger(__name__)


class ScoringFile:
    """Represents a single scoring file in the PGS Catalog.

    :param identifier: A PGS Catalog score accession in the format ``PGS123456`` or a path to a local scoring file
    :param target_build: An optional :class:`GenomeBuild`, which represents the build you want the scoring file to align to
    :param query_result: An optional :class:`ScoreQueryResult`, if provided with an accession identifier it prevents hitting the PGS Catalog API
    :raises pgscatalog.corelib.InvalidAccessionError: If the PGS Catalog API can't find the provided accession
    :raises pgscatalog.corelib.ScoreFormatError: If you try to iterate over a ``ScoringFile`` without a local path (before downloading it)

    You can make ``ScoringFiles`` with a path to a scoring file with minimal metadata:

    >>> from pgscatalog.core.lib.genomebuild import GenomeBuild
    >>> from pgscatalog.core.lib._config import Config
    >>> sf = ScoringFile(Config.ROOT_DIR / "tests" / "data" / "custom.txt")
    >>> sf # doctest: +ELLIPSIS
    ScoringFile('.../custom.txt', target_build=None)
    >>> sf.header
    ScoreHeader(pgs_id='test', pgs_name='test', trait_reported='test trait', genome_build=GenomeBuild.GRCh37)
    >>> sf.is_harmonised
    False

    Also supports PGS Catalog header metadata:

    >>> sf = ScoringFile(Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz")
    >>> sf # doctest: +ELLIPSIS
    ScoringFile('.../PGS000001_hmPOS_GRCh38.txt.gz', target_build=None)
    >>> sf.header
    CatalogScoreHeader(pgs_id='PGS000001', pgs_name='PRS77_BC', trait_reported='Breast cancer', genome_build=None, format_version=<ScoreFormatVersion.v2: '2.0'>, trait_mapped=['breast carcinoma'], trait_efo=['EFO_0000305'], variants_number=77, weight_type=None, pgp_id='PGP000001', citation='Mavaddat N et al. J Natl Cancer Inst (2015). doi:10.1093/jnci/djv036', HmPOS_build=GenomeBuild.GRCh38, HmPOS_date=datetime.date(2022, 7, 29), HmPOS_match_pos='{"True": null, "False": null}', HmPOS_match_chr='{"True": null, "False": null}')

    Looking at the header above, the original submission lacked a genome build but has been harmonised:

    >>> sf.is_harmonised
    True

    >>> sf.genome_build
    GenomeBuild.GRCh38

    >>> sf.pgs_id
    'PGS000001'

    >>> for variant in sf.variants: # doctest: +ELLIPSIS
    ...     variant
    ...     break
    ScoreVariant(rsID='rs78540526', chr_name='11', chr_position=None, effect_allele=Allele(allele='T', is_snp=True)...

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
    FileNotFoundError: self.local_path=None: did you remember to .download()?

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
    pgscatalog.core.lib.pgsexceptions.InvalidAccessionError: Invalid accession: 'potato'

    The same exception is raised if you provide a well formatted identifier that doesn't exist:

    >>> with tempfile.TemporaryDirectory() as tmp_dir:
    ...     ScoringFile("PGS000000", GenomeBuild.GRCh38).download(tmp_dir)
    Traceback (most recent call last):
    ...
    pgscatalog.core.lib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGS000000'
    """

    def __init__(self, identifier, target_build=None, query_result=None, **kwargs):
        self._target_build = target_build
        if query_result is None:
            self._identifier = identifier
        else:
            self._identifier = query_result

        try:
            # let's try parsing a PGS Catalog header
            self._header = models.CatalogScoreHeader.from_path(self._identifier)
            logger.info(f"{identifier}: Valid PGS Catalog header")
            self.local_path = pathlib.Path(self._identifier)
            self._init_from_path(target_build=target_build)
        except ValidationError:
            logger.warning("PGS Catalog header not detected in scoring file")
            # that didn't work, let's try parsing a basic score header
            self._header = models.ScoreHeader.from_path(self._identifier)
            logger.info(f"{identifier}: Valid simple score header")
            self.local_path = pathlib.Path(self._identifier)
            self._init_from_path(target_build=target_build)
        except (FileNotFoundError, TypeError):
            # was it an accession?
            self.include_children = kwargs.get("include_children", None)
            self._init_from_accession(self._identifier, target_build=target_build)

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

        self._pgs_id = score.pgs_id
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
        self._pgs_id = self._header.pgs_id

    @property
    def is_wide(self) -> bool:
        # TODO: temporarily disabled because wide scoring files are so far untested and undocumented
        return False

    @property
    def pgs_id(self) -> str:
        return self._pgs_id

    @property
    def is_harmonised(self) -> bool:
        return self._header.is_harmonised

    @property
    def genome_build(self) -> GenomeBuild:
        if self.is_harmonised:
            return self._header.HmPOS_build
        else:
            return self._header.genome_build

    @property
    def header(self):
        return self._header

    @property
    def variants(self) -> Generator[ScoreVariant, None, None]:
        """A generator that yields rows from the scoring file as ``ScoreVariants``,
        if a local file is available (i.e. after downloading). Always available for
        class instances that have a valid local path.

        >>> from pgscatalog.core.lib._config import Config
        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000802_hmPOS_GRCh37.txt"
        >>> sf = ScoringFile(testpath)
        >>> for variant in sf.variants:
        ...     variant
        ...     break
        ScoreVariant(rsID='rs10936599', chr_name='3', chr_position=170974795...
        """
        return self.read_variants()

    @property
    def target_build(self) -> GenomeBuild:
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
        >>> from pgscatalog.core.lib import GenomeBuild

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
        directory = pathlib.Path(directory)

        if not directory.is_dir():
            raise NotADirectoryError(directory)

        out_path = directory / pathlib.Path(self.path).name

        logger.debug(f"Downloading {self.path} to {out_path}")

        https_download(
            url=self.path,
            out_path=out_path,
            overwrite=overwrite,
            directory=directory,
        )

        # update local file attributes
        self.local_path = out_path

    def normalise(
        self, liftover=False, drop_missing=False, chain_dir=None, target_build=None
    ):
        """Extracts key fields from a scoring file in a normalised format.

        Takes care of quality control.

        >>> from pgscatalog.core.lib import GenomeBuild
        >>> from pgscatalog.core.lib._config import Config
        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000001_hmPOS_GRCh38.txt.gz"
        >>> variants = ScoringFile(testpath).normalise()
        >>> for x in variants: # doctest: +ELLIPSIS
        ...     x
        ...     break
        ScoreVariant(rsID='rs78540526', chr_name='11', chr_position=69516650, effect_allele=Allele(allele='T', is_snp=True), ...

        Supports lifting over scoring files from GRCh37 to GRCh38:

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "lift_to_grch38.txt"
        >>> chaindir = Config.ROOT_DIR / "tests" / "data" / "chain"
        >>> sf = ScoringFile(testpath)
        >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh38)
        >>> for x in variants:
        ...     (x.rsID, x.chr_name, x.chr_position)
        ...     break
        ('rs78540526', '11', 69516650)

        Example of lifting down (GRCh38 to GRCh37):

        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "lift_to_grch37.txt"
        >>> chaindir = Config.ROOT_DIR / "tests" / "data" / "chain"
        >>> sf = ScoringFile(testpath)
        >>> variants = sf.normalise(liftover=True, chain_dir=chaindir, target_build=GenomeBuild.GRCh37)
        >>> for x in variants:
        ...     (x.rsID, x.chr_name, x.chr_position)
        ...     break
        ('rs78540526', '11', 69331418)

        Liftover support is only really useful for custom scoring files that aren't
        in the PGS Catalog. It's always best to use harmonised data when it's
        available from the PGS Catalog. Harmonised data goes through a lot of validation
        and error checking.

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

    @contextmanager
    def read(self) -> Iterator[csv.DictReader]:
        """A simple method of reading variants from a scoring file.

        Returns a csv.DictReader, so each row is a variant in a dictionary.

        No data validation is done. Combine the returned dictionaries with the pydantic models if you want to do that (CatalogScoreVariants).

        This method must be called with a context manager:

        >>> from pgscatalog.core.lib._config import Config
        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000802_hmPOS_GRCh37.txt"
        >>> sf = ScoringFile(testpath)
        >>> with sf.read() as reader:
        ...     for variant in reader:
        ...         variant
        ...         break
        {'rsID': 'rs10936599', 'chr_name': '3', 'chr_position': '170974795', 'effect_allele': 'T', 'other_allele': 'C', 'effect_weight': '0.123', 'allelefrequency_effect': '0.377', 'is_dominant': 'True', 'is_recessive': 'False', 'locus_name': 'MYNN', 'hm_source': 'ENSEMBL', 'hm_rsID': 'rs10936599', 'hm_chr': '3', 'hm_pos': '169492101', 'hm_inferOtherAllele': ''}

        Calling this method directly isn't helpful:

        >>> sf.read()
        <contextlib._GeneratorContextManager object ...>

        Only local scoring files can be read (download them first):

        >>> sf = ScoringFile("PGS001229")
        >>> with sf.read() as f:
        ...     pass
        Traceback (most recent call last):
        ...
        FileNotFoundError: self.local_path=None: did you remember to .download()?
        """
        if self.local_path is None:
            raise FileNotFoundError(
                f"{self.local_path=}: did you remember to .download()?"
            )

        with xopen(self.local_path) as f:
            yield csv.DictReader(
                (line for line in f if not line.startswith("#")), delimiter="\t"
            )

    def read_variants(self) -> Generator[ScoreVariant, None, None]:
        """Yields rows from a scoring file as ScoreVariants

        ScoreVariants are pydantic models with data validation (PGS Catalog standards)

        >>> from pgscatalog.core.lib._config import Config
        >>> testpath = Config.ROOT_DIR / "tests" / "data" / "PGS000802_hmPOS_GRCh37.txt"
        >>> sf = ScoringFile(testpath)
        >>> variants = sf.read_variants()
        >>> for i, variant in enumerate(variants):
        ...     variant
        ...     if i == 2:
        ...         break
        ScoreVariant(rsID='rs10936599', chr_name='3', ...
        ScoreVariant(rsID='rs6061231', chr_name='20', ...
        ScoreVariant(rsID='rs10774214', chr_name='12', ...
        """
        with self.read() as reader:
            for row_nr, variant in enumerate(reader):
                yield ScoreVariant(
                    **variant, **{"accession": self.pgs_id, "row_nr": row_nr}
                )


class ScoringFiles:
    """This container class provides methods to work with multiple ScoringFile objects.

    You can use publications or trait accessions to instantiate:

    >>> from pgscatalog.core.lib import GenomeBuild
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
    pgscatalog.core.lib.pgsexceptions.InvalidAccessionError: No Catalog result for accession 'PGPpotato'

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
