"""The cache: the portable interface between stages, and the shape checks on it.

The cache is the zarr store written by LOAD and read by every later stage. Its
hierarchy is::

    /{sampleset}/{chromosome}/{filename}/genotypes
                                        /meta

Stages deliberately fan out at different levels of it - MATCH per chromosome,
SCORE per shard - so this module is the one place that knows how those levels are
arranged. Callers ask it which shards a staged cache holds, which chromosome, which
single shard, whose samples those shards hold and what variant metadata they carry;
nobody outside it slices a path, counts a depth or walks the tree a second time
(ADR 0003).

The shape invariants live here too, because they are checks on the cache. They
read zarr attributes and variant metadata only: nothing here touches the 3D
genotype arrays.
"""

from __future__ import annotations

import itertools
import pathlib
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Any, cast

import zarr
import zarr.storage

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    import duckdb

    from pgscatalog.calc.lib.types import Pathish

# what an absent variant means when LOAD declared nothing (ADR 0004)
_DEFAULT_ABSENCE_SEMANTICS = "unobserved"


def leaf_group_path(sampleset: str, chromosome: str | None, filename: str) -> str:
    """Where a shard's genotypes land in the cache hierarchy.

    Built here and nowhere else, because LOAD writes this group and samplesheet
    validation has to check the string LOAD will actually write to. Two functions
    building the rule separately agree only by coincidence, and if they drifted the
    guard would pass while two shards of one chromosome sharing a basename collided
    at merge time, one leaf silently overwriting the other (ADR 0003).

    Parameters
    ----------
    sampleset : str
        The sampleset the shard belongs to.
    chromosome : str | None
        The chromosome the input was declared for, or None for unsplit whole-genome
        input, which occupies the sentinel chromosome `ALL` so the depth is the same
        whatever the input shape.
    filename : str
        The target genome file's basename. Directories are stripped, which is what
        makes the uniqueness validation necessary.
    """
    return f"{sampleset}/{chromosome or 'ALL'}/{filename}"


def meta_group_path(leaf_path: str) -> str:
    """Where a leaf keeps its 1D variant metadata.

    Built here and nowhere else, because LOAD writes this group and every later
    stage reads it.

    Parameters
    ----------
    leaf_path : str
        The leaf group's path relative to the store root.
    """
    return f"{leaf_path}/meta"


@dataclass(frozen=True)
class Shard:
    """One leaf of the cache hierarchy: the cached genotypes of one target genome.

    Parameters
    ----------
    path : str
        The leaf group's path relative to the store root, which is how the shard
        is identified downstream (the `zarr_group` column).
    chromosome : str
        The chromosome group this leaf sits under, or `ALL` for unsplit input.
        Taken from the shard's position in the tree, never from its path.
    attrs : dict[str, Any]
        The leaf group's zarr attributes, as LOAD wrote them.
    """

    path: str
    chromosome: str
    attrs: dict[str, Any]

    @property
    def samples(self) -> list[str]:
        """The sample IDs this shard holds"""
        return cast("list[str]", self.attrs["samples"])

    @property
    def absence_semantics(self) -> str:
        """What an absent variant means here, as declared per sampleset (ADR 0004)"""
        return cast(
            "str", self.attrs.get("absence_semantics", _DEFAULT_ABSENCE_SEMANTICS)
        )


class GenotypeCache:
    """The cache a stage has been staged with: a zarr store on a local filesystem.

    Parameters
    ----------
    path : Pathish
        Path to the zarr store (`genotypes.zarr` in the stage's directory).
    sampleset : str | None
        The sampleset name. MATCH_COMBINE and SCORE_COMBINE have no store to read
        it from, so they pass it in.
    """

    def __init__(self, path: Pathish, sampleset: str | None = None) -> None:
        self._path = pathlib.Path(path)
        self._declared_sampleset = sampleset

    # ponytail: `root` still hands out a raw zarr group because SCORE reads the 3D
    # genotype arrays through it (`make_dosage_arrays`), which is an array read
    # rather than a walk of the hierarchy. It goes when dosage takes a shard instead
    # of the root group, which splitting the stages left alone.
    @property
    def root(self) -> zarr.Group:
        """The store's root group, opened read only"""
        store = zarr.storage.LocalStore(self._path)
        return zarr.open_group(store=store, mode="r")

    @property
    def _sampleset_group(self) -> zarr.Group:
        """The group holding this sampleset's chromosomes"""
        return cast("zarr.Group", self.root[self.sampleset])

    @cached_property
    def sampleset(self) -> str:
        """The sampleset this cache holds.

        One per cache: shards of one sampleset are slices of one callset, so a
        second sampleset in the store means two cohorts were staged together.
        """
        if self._declared_sampleset is not None:
            return self._declared_sampleset

        samplesets = list(self.root.keys())

        if len(samplesets) > 1:
            raise NotImplementedError(
                f"One than one sampleset not supported, found {samplesets}"
            )
        if len(samplesets) == 0:
            raise FileNotFoundError(
                "zarr archive looks badly structured, can't find sampleset group"
            )
        return samplesets[0]

    @cached_property
    def shards(self) -> list[Shard]:
        """Every shard this cache holds.

        One traversal, cached: a stage asks this several times and the store does
        not change under it.

        The chromosome is the group a leaf sits under rather than a slice of its
        path, and leaves are found by carrying sample IDs rather than by being at
        a known depth - they are the only nodes that do (ADR 0003).
        """
        return [
            Shard(path=node.path, chromosome=chromosome, attrs=dict(node.attrs))
            for chromosome, chromosome_group in self._sampleset_group.groups()
            for _, node in chromosome_group.members(max_depth=None)
            if "samples" in node.attrs
        ]

    @property
    def chromosome(self) -> str:
        """The single chromosome this cache holds, which is what a MATCH task is
        staged with (ADR 0003)"""
        return _one_chromosome(self.shards)

    @property
    def shard(self) -> Shard:
        """The single leaf group this SCORE task holds.

        weight_mat_row_nr indexes this task's dosage array, so staging two shards
        at once would silently misalign it.
        """
        if len(self.shards) != 1:
            paths = sorted(shard.path for shard in self.shards)
            raise ValueError(
                f"SCORE runs per shard but this cache holds {paths}. "
                "Stage one shard per task."
            )
        return self.shards[0]

    @cached_property
    def sample_ids(self) -> list[str]:
        """The sample IDs this cache's shards hold, in genotype array order.

        Every shard has to agree, which is why this is one list rather than a union.
        A SCORE task is staged with one shard, so there is nothing left to disagree
        with: the assumption stopped being load-bearing rather than being fixed, and
        the check is kept as it was (ADR 0003).
        """
        samples: list[str] | None = None

        for shard in self.shards:
            if samples is None:
                # initialise the sample list
                samples = shard.samples
            elif samples != shard.samples:
                # compare sample names
                raise NotImplementedError(
                    "Different samples across target genome files"
                )

        if samples is None:
            raise TypeError("Couldn't get samples from zarr attribute")

        return samples

    def variant_metadata(self) -> Iterator[tuple[str, dict[str, Any]]]:
        """Each shard's 1D variant metadata, paired with the shard's path.

        Read one shard at a time, because these arrays get big - especially once
        they're loaded into a dataframe - and MATCH only needs one shard's worth in
        memory at a time.
        """
        root = self.root

        for shard in self.shards:
            meta = cast("zarr.Group", root[meta_group_path(shard.path)])
            # loads the numpy arrays into memory now
            yield shard.path, {name: array[:] for name, array in meta.arrays()}


def _one_chromosome(shards: Sequence[Shard]) -> str:
    """The single chromosome a set of shards belongs to.

    MATCH is staged with every shard of one chromosome, so a second chromosome is
    a fan-out wiring mistake rather than a bigger job (ADR 0003).
    """
    chromosomes = {shard.chromosome for shard in shards}
    if len(chromosomes) != 1:
        raise ValueError(
            f"MATCH runs per chromosome but this cache holds {sorted(chromosomes)}. "
            "Stage one chromosome's shards per task."
        )
    return chromosomes.pop()


def check_cache_shape(
    conn: duckdb.DuckDBPyConnection, shards: Sequence[Shard], sampleset: str
) -> str:
    """Check the shards under this chromosome can be scored independently.

    Reads zarr attributes and variant metadata only - nothing touches the 3D
    genotype arrays. Returns the chromosome.

    Parameters
    ----------
    conn : duckdb.DuckDBPyConnection
        A connection to the score database, which holds the `targetvariants` table
        loaded from this chromosome's shards.
    shards : Sequence[Shard]
        The shards the MATCH task was staged with.
    sampleset : str
        The sampleset they belong to, named if there are no shards at all.
    """
    if not shards:
        raise ValueError(f"No shards found under {sampleset!r}")

    chromosome = _one_chromosome(shards)

    # the invariants are about a set of leaves keyed by path, which is also what
    # `targetvariants.zarr_group` holds, so they work on the mapping
    attrs = {shard.path: shard.attrs for shard in shards}

    check_absence_semantics(attrs)
    check_samples_disjoint(attrs)
    check_shards_are_rectangular(conn, attrs)

    return chromosome


def check_absence_semantics(shards: dict[str, dict[str, Any]]) -> None:
    """Absence semantics is a property of the callset, so every shard of it must
    declare the same thing, and only `unobserved` is implemented (ADR 0004)"""
    declared = {
        path: attrs.get("absence_semantics", _DEFAULT_ABSENCE_SEMANTICS)
        for path, attrs in shards.items()
    }

    if len(set(declared.values())) > 1:
        raise ValueError(
            f"Shards disagree about absence semantics: {declared}. Shards of one "
            "sampleset are slices of one callset, so they cannot mean different "
            "things by an absent variant."
        )

    if (semantics := next(iter(declared.values()))) != _DEFAULT_ABSENCE_SEMANTICS:
        raise NotImplementedError(
            f"absence_semantics={semantics!r} is not implemented yet. Absent "
            "variants are always treated as never assessed in this release."
        )


def check_samples_disjoint(shards: dict[str, dict[str, Any]]) -> None:
    """Shards under one chromosome must not share samples.

    SCORE_COMBINE sums partial scores grouped by sample ID, which is only valid
    while every (sample, variant) pair appears in exactly one shard. A shared
    sample double-counts its dosage: no error, no NaN, just an inflated score.
    """
    # ponytail: O(n^2) over shards-per-chromosome, which is a handful of sample
    # batches. If that ever stops being true, count sample occurrences in one pass.
    samples = {path: set(attrs["samples"]) for path, attrs in shards.items()}

    for (left, left_samples), (right, right_samples) in itertools.combinations(
        samples.items(), 2
    ):
        if shared := left_samples & right_samples:
            examples = ", ".join(sorted(shared)[:3])
            raise ValueError(
                f"Shards {left!r} and {right!r} share {len(shared)} samples "
                f"(e.g. {examples}). Shards under one chromosome must hold disjoint "
                "samples: the chromosome is the variant axis, so splitting variants "
                "within one is not supported. Split by chromosome instead."
            )


def check_shards_are_rectangular(
    conn: duckdb.DuckDBPyConnection, shards: dict[str, dict[str, Any]]
) -> None:
    """Shards under one chromosome must agree on which variants they contain.

    Defined over rows with real data: the cache back-fills queried-but-absent
    positions as sentinel rows with a null ref, so raw position sets are identical
    across shards by construction and would never look ragged.
    """
    if len(shards) < 2:
        return

    # targetvariants is keyed on (variant_id, zarr_group), so a variant can appear at
    # most once per shard and counting its rows counts the shards holding it. Anything
    # short of every shard is a variant one shard assessed and another didn't.
    ragged = conn.execute(
        """
        SELECT variant_id, list(zarr_group ORDER BY zarr_group) AS present_in
        FROM targetvariants
        WHERE ref IS NOT NULL AND zarr_group IN $shards
        GROUP BY variant_id
        HAVING COUNT(*) <> $n_shards
        ORDER BY variant_id
        LIMIT 3
        """,
        {"shards": list(shards), "n_shards": len(shards)},
    ).fetchall()

    if not ragged:
        return

    # only on the error path, so the happy path pays for one aggregate and no more
    counts = conn.execute(
        """
        SELECT zarr_group, COUNT(*)
        FROM targetvariants
        WHERE ref IS NOT NULL AND zarr_group IN $shards
        GROUP BY zarr_group
        ORDER BY zarr_group
        """,
        {"shards": list(shards)},
    ).fetchall()

    summary = ", ".join(f"{path} ({n} variants)" for path, n in counts)
    examples = ", ".join(
        f"{variant_id} (only in {', '.join(present_in)})"
        for variant_id, present_in in ragged
    )
    raise ValueError(
        "Shards under one chromosome contain different variants: "
        f"{summary}. For example: {examples}. Every shard of a chromosome must "
        "interrogate the same variants, because a variant missing from one shard and "
        "present in another cannot be told apart from a variant that was never "
        "assessed there. Re-export the target genomes without dropping sites that are "
        "monomorphic within a batch, or split the input by chromosome rather than by "
        "sample."
    )
